import xlsxwriter
import copy

from mutant import Mutant, PointMutant, InsertionMutant, DeletionMutant
from codontable import Translate
from barcodeset import FormatBarcode
from mutant import SpacedAa

class WorksheetWriter:
	def __init__(self,worksheet):
		self.worksheet=worksheet
		self.row=0
	def write(self,*args):
		try:
			if type(args[0])==type([]):
				self.worksheet.write_row(self.row,0,args[0])
			elif type(args[0])==type(()):
				self.worksheet.write_row(self.row,0,args[0])
			else:
				self.worksheet.write_row(self.row,0,args)
			self.row+=1
		except TypeError:
			print args
			raise

class MutantSet:
	"""
	MutantSet contains the base sequence and parameters for a mutant
	screen.

	Nucleotide 1 of each mutation range must be the first 
	nucleotide of a codon.
	base sequence ----------------------------------------------------
	synth sequence        --------------------------------------
	mutated seq                  ---------------------------
	5' overhang seq       -------
	3' overhang seq                                         ----
	The synth sequence may not be in frame because the 5' overhang may
	not be a multiple of 3 in length, but the mutate sequence must
	be in frame. The mutate sequence is represented by a Mutant
	object, as are all the mutants derived from it.
	
	Revised 7/13/2017
	Scaling down screen due to lab technical issues. Now using few
	bar codes and re-using them.
	Specifying in frame insertions (ifi's) and in-frame deletions (ifd's).
	If ifi's or ifd's not specified then generate all of them, otherwise
	only use those that are specified.
	7/27/2017 - added ms_barcode_choices and ifd_barcode_choices parameters.
	These let you specify a list of barcodes for ms and ifi/ifd mutants
	respectively.
	"""
	def __init__(self,screen,fragment,basesequence,synth_from,synth_to,mut_from,mut_to,first_aa_coord,num_barcodes=0,ms=[],ifi=[],ifd=[],ms_barcode_choices=[],id_barcode_choices=[]):
		# Name of the MutantSet.
		self.screen=screen
		self.fragment=fragment
		# Sequence synthesis coords, 1-based, applies to entire 
		# base sequence.
		self.synth_from=synth_from
		self.synth_to=synth_to
		# Sequence mutagenesis coords, 1-based, applies to entire
		# base sequence.
		self.mut_from=mut_from
		self.mut_to=mut_to
		# The cooridinate in the WT protein of the first amino acid
		# that is getting mutated. This is used for the GVCF notation
		# of the protein mutations.
		self.first_aa_coord=first_aa_coord

		# Store the number of bar codes per mutation to generate.
		self.num_barcodes=num_barcodes

		# Store the lists of missense, in-frame insertions (ifi's) and 
		# in-frame deletions (ifd's) to be generated. If these
		# lists are empty then generate all of them.
		self.ms=ms[:]
		self.ifi = ifi[:]
		self.ifd = ifd[:]

		# Sequence being synthesized.
		self.synth_nt_sequence=basesequence[self.synth_from-1:self.synth_to]

		# Overhang sequences before and after coding regions.
		self.overhang_seq_5 = basesequence[self.synth_from-1:self.mut_from-1]
		self.overhang_seq_3 = basesequence[self.mut_to:self.synth_to]

		# Sequence being mutated.
		self.coding_sequence = Mutant(self.screen,self.fragment,basesequence[self.mut_from-1:self.mut_to],self.overhang_seq_5,self.overhang_seq_3,self.first_aa_coord)

		# The overhang_seq_5 + coding_sequence.nt_sequence + overhang_seq_3 
		# should be identical to the synthesized sequence.
		assert len(self.overhang_seq_5) + len(self.coding_sequence.coding_sequence) + len(self.overhang_seq_3) == len(self.synth_nt_sequence), "Length of fragments (%d, %d, %d) != total synthesized length (%d)" % ( len(self.overhang_seq_5), len(self.coding_sequence.coding_sequence) , len(self.overhang_seq_3), len(self.synth_nt_sequence))

		# Mutagenesis coords, 0-based, applies to synthesized
		# sequence.
		self.mut_from_0=self.mut_from - self.synth_from
		self.mut_to_0=self.mut_from_0 + (self.mut_to - self.mut_from)

		# List of mutants to be generated.
		self.mutants=[]
	
	def __str__(self):
		return "Screen: %s, fragment: %s, synth_from: %d, synth_to: %d, mut_from: %d, mut_to: %d, synth_nt_sequence %s...%s, mut_nt_sequence %s...%s, aa_sequence %s, mut_from_0: %d, mut_to_0: %d, first_aa_coord: %d" % ( self.screen, self.fragment, self.synth_from, self.synth_to, self.mut_from, self.mut_to, self.synth_nt_sequence[0:5], self.synth_nt_sequence[-5:], self.coding_sequence.coding_sequence[0:5], self.coding_sequence.coding_sequence[-5:], self.coding_sequence.aa_sequence, self.mut_from_0, self.mut_to_0,self.first_aa_coord)
	
	def Aa(self,i):
		"""
		Returns amino acid i (zero-based) from the amino acid sequence.
		"""
		return self.coding_sequence.aa_sequence[i]

	def GenerateMutants(self):
		"""
		Populates list of all Mutants in this MutantSet. These represent
		the missense substitutions possible via single nucleotide 
		substitutions, not all single nucleotide substitutions.
		"""
		print "Generating mutants for mutant set %s %s." % (self.screen, self.fragment)
		self.mutants=[self.coding_sequence]
		self.mutants += self.coding_sequence.PointMutants(self.ms)
		self.mutants += self.coding_sequence.InsertionMutants(self.ifi)
		self.mutants += self.coding_sequence.DeletionMutants(self.ifd)
		print "Generated %d mutants." % len(self.mutants)
	
	def FilterMutants(self):
		"""
		Applies filters to mutants, discarding any that don't
		pass filters.
		"""
		print "Applying filters to mutant set %s %s." % (self.screen, self.fragment)
		passed_filters=0
		for m in self.mutants:
			m.ApplyFilters()
			if m.passes_filters:
				passed_filters+=1
		print "%d mutants pass filters." % passed_filters
	
	def CreateBarcodes(self):
		"""
		Creates bar codes for each Mutant object. Replaces
		this mutantset's list of mutants with a new and
		expanded list of bar coded mutants.
		"""
		print "Creating bar codes for mutant set %s %s." % (self.screen, self.fragment)
		orig_mutant_list = self.mutants[:]
		self.mutants=[]
		# Create bar coded versions of mutants for each mutant in the list.
		for mutant in orig_mutant_list:
			self.mutants.append(mutant)
			if mutant.passes_filters:
				self.mutants += mutant.GenerateBarcodedMutants()
			else:
				print "NOT generating bar codes for %s mutant %s." % ( mutant.type, mutant.name )

		# Then, for each bar coded mutant, create a copy of the WT
		# coding sequence that has the same bar code.
		print "Creating barcoded copies of WT sequence."
		barcode_wt = []
		barcode_problems=0
		for mutant in self.mutants:
			if mutant.barcode:
				barcode = mutant.AdjustedBarcode()
				wt=copy.deepcopy(self.coding_sequence)
				print "Applying barcode %s from %s to WT." % ( FormatBarcode(barcode), mutant.name )
				try:
					wt.ApplyBarcode(barcode)
					wt.ApplyFilters()
					barcode_wt.append(wt)
				except AssertionError:
					barcode_problems += 1
					print "AssertionError:"
					print "aa sequence changed after bar code applied."
					print "%-5s: %s" % ('mut', SpacedAa(mutant.aa_sequence))
					print "%-5s: %s" % ('mut', mutant.coding_sequence)
					print "%-5s: %s" % ('orig', SpacedAa(self.coding_sequence.aa_sequence))
					print "%-5s: %s" % ('orig', self.coding_sequence.coding_sequence)
					print "%-5s: %s" % ('alt', SpacedAa(Translate(wt.coding_sequence)))
					print "%-5s: %s" % ('alt', wt.coding_sequence)
					raise
		print "Generated %d barcoded WT sequences, %d problems." % (len(barcode_wt), barcode_problems)
		self.mutants += barcode_wt

		print "Done - mutant set contains %d mutants." % len(self.mutants)
	
	def Write(self,version):
		"""
		Write() generates an excel file for the MutantSet.
		version is the version number of the package.
		"""
		filename=self.screen+'_'+self.fragment+".xlsx"
		workbook=xlsxwriter.Workbook(filename,{'strings_to_numbers':  True,})

		# Summary information.
		summary = WorksheetWriter(workbook.add_worksheet("Summary"))
		summary.write("screen",self.screen)
		summary.write("fragment",self.fragment)
		summary.write("version",version)
		summary.write("synth_from",self.synth_from)
		summary.write("synth_to",self.synth_to)
		summary.write("mut_from",self.mut_from)
		summary.write("mut_to",self.mut_to)
		summary.write("synth_nt_sequence",self.synth_nt_sequence)
		summary.write("mut_nt_sequence",self.coding_sequence.coding_sequence)
		summary.write("aa_sequence",self.coding_sequence.aa_sequence)
		summary.write("mut_from_0",self.mut_from_0)
		summary.write("mut_to_0",self.mut_to_0)
		
		pages = [
			('Point Mutations', 'point'),
			('In-frame Insertions', 'insertion'),
			('In-frame Deletions', 'deletion'),
			('WT Controls', 'WT'),
		]
		for (worksheet_name,mutant_type) in pages:
			page=WorksheetWriter(workbook.add_worksheet(worksheet_name))
			page.worksheet.freeze_panes(1, 0)
			page.write(Mutant.TableHeader().split())
			for m in filter(lambda x:x.type==mutant_type,self.mutants):
				if m.barcode or not m.passes_filters:
					page.write(m.AsTuple())
		workbook.close()

	def AsFasta(self):
		"""
		AsFasta() writes the mutant set as a fasta file.
		"""
		filename = self.screen + '_' + self.fragment + ".fasta"
		print "Writing fasta file for mutant set %s %s to %s." % ( self.screen, self.fragment, filename )
		ofs=open(filename,"w")
		#ofs.write(self.__str__())
		#ofs.write("\n")
		for m in self.mutants:
			#if m.barcode or not m.passes_filters:
			#	ofs.write(m.AsFasta())
			ofs.write(m.AsFasta())
		ofs.close()
	
