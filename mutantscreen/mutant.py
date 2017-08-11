import re
import copy

from barcodeset import BarcodeSet, FormatBarcode, ParseBarcode
from codontable import codons, reverse_codons, Translate

def FastaSequence(s):
	"""
	Returns sequence with newlines inserted every 60 bases for 
	displaying nucleotide sequences as fasta.
	"""
	retval = ''
	while s:
		retval += s[0:60] + "\n"
		s=s[60:]
	return retval
		
def Deletion(seq,position,length=3):
	"""
	Returns sequence seq with a deletion of length 'length' at
	position 'position'.
	"""
	return seq[0:position] + seq[position+length:]

def Insert(s1,s2,pos):
	"""
	Inserts s2 into s1 and position pos.
	"""
	return s1[0:pos] + s2 + s1[pos:]

def Substitute(s1,s2,pos):
	"""
	Substitutes s2 for character at position pos in s1.
	"""
	return s1[0:pos] + s2 + s1[pos+1:]

def BooleanAsIntString(value):
	return str( int( value ))

def MutateSequence(sequence,nucleotide,location):
	return sequence[0:location]+nucleotide+sequence[location+1:]

def GCContent(sequence):
	"""
	Returns the GC content of a sequence.
	"""
	return (sequence.count('C')+sequence.count('G'))/float(len(sequence))

def HighGC(sequence):
	"""
	Retuns True if sequence has GC content > 75%.
	"""
	return GCContent(sequence) > 0.75

def LowGC(sequence):
	"""
	Retuns True if sequence has GC content < 25%.
	"""
	return GCContent(sequence) < 0.25

def SpacedAa(aa_sequence):
	"""
	Takes an amino acid sequence and returns it s.t. there are two 
	spaces after each amino acid letter. This lets us print a.a. sequences
	along side their coding n.t equivalent.
	"""
	return ''.join(map(lambda x:x+'  ',list(aa_sequence)))

class Mutant:
	"""
	Mutant is a nucleotide sequence AND the amino acid sequence
	that it encodes - this represents the nucleotide / protein pair.
	A requirement of this is that the length of the nucleotide sequence
	is a multiple of 3, and is non-zero.
	"""
	# Indexes into the filter_results list:
	HasNoMonoSixmers=0
	HasNoGFourMers=1
	HasOKGcContent=2

	def __init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent=None):
		assert len(coding_sequence)>0 and len(coding_sequence)%3==0
		self.screen = screen
		self.fragment = fragment
		self.coding_sequence=coding_sequence
		self.overhang_seq_5 = overhang_seq_5
		self.overhang_seq_3 = overhang_seq_3
		self.complete_sequence = overhang_seq_5 + self.coding_sequence + overhang_seq_3
		self.first_aa_coord=first_aa_coord
		self.aa_sequence=Translate(self.coding_sequence)
		self.mutated_codon=None
		# Filter results - this is list of boolean values.
		self.filter_results=[False]*3
		# Pass_filters: whether mutant passes all filters.
		# Defaults to True until mutant fails a filter.
		self.passes_filters = True
		#self.barcode = '_'
		self.barcode = ''

		# Fields to describe the mutation that this Mutant
		# possesses. All these fields will be None if the Mutant
		# represents the wild type.
		self.mut_aa = ''
		self.ref_aa = ''
		self.ref_codon = ''
		self.mut_codon = ''
		# mut_codon_coords is a list of the nucleotide positions
		# that are in the mutated codon. This gets populated in 
		# InitMutation(). Meaning:
		# point: coords of exon containing point mutation
		# deletion: - none
		# insertion: coords of inserted exon
		self.mut_codon_coords=[]
		self.parent = parent
		self.mut_location = None
		# mutant name is amino acid change.
		self.name = "WT"
		self.type = 'WT'
		self.ref_nucleotide = ''
		self.mut_nucleotide = ''

	def AdjustedBarcode(self):
		"""
		Returns the parsed bar code. If this mutant contains an
		insertion or deletion then the coordinates following the
		mutation are adjusted. The intent is that this bar code
		can be applied to the WT sequence.
		"""
		return ParseBarcode(self.barcode)
		
	def FastaHeader(self):
		return '>%s %s %s %s %s' % (self.screen, self.fragment, self.type, self.name, self.barcode )
		
	def ProteinMutation(self):
		"""
		Returns impact of mutation on protein. For wild type sequence
		the impact is 'None'. 
		For point mutations:
			native aa, aa position, mutant aa, e.g. "V45M". 
		For deletions:
			native aa, aa position, 'X', e.g. "V45X". 
		For insertions:
			
		"""
		return self.name


	def IsEfficientTranslator(self):
		"""
		Checks if this mutant's codon has a frequency as good or
		better than the wild type codon.
		"""
		ref_codon_frequency = codons[self.ref_codon][1]
		mut_codon_frequency = codons[self.mut_codon][1]

		return mut_codon_frequency >= ref_codon_frequency
	
	def RNASecondaryStructureOK(self):
		"""
		"""
		pass
	
	def PositionNearMutation(self,position):
		"""
		Checks if the given position in the nucleotide sequence
		is in the vicinity of the mutation. This test is used to
		avoid using bar code nucleotides close to the mutation.
		"""
		return position in self.mut_codon_coords

	def GenerateBarcodedMutants(self):
		"""
		New version of this function uses the suggested dinucleotide
		bar codes and the num_barcodes value to determine how many
		new mutants to create.
		"""
		newmutants=[]
		print "Generating barcodes for %s mutant %s." % (self.type, self.name)
		if self.type=='point':
			barcode_choices=self.ms_barcode_choices
		else:
			barcode_choices=self.id_barcode_choices

		for barcode in barcode_choices:
			# Make a copy of self.
			bcm=copy.deepcopy(self)
			try:
				bcm.ApplyBarcode(bs.next())
			except IndexError:
				print "Warning - ran out of bar codes from pool of %d for mutant %s. Created %d bar coded mutants." % ( bs.Size(), self.name, len(newmutants) )
				print Mutant.TableHeader()
				print self.AsTabDelim()
				print self
				break
			bcm.ApplyFilters()
			if bcm.passes_filters:
				newmutants.append(bcm)
			else:
				print "Warning: bar coded mutant failed filters."
			for 
		return newmutants

	def OldGenerateBarcodedMutants(self):
		"""
		Generates bar codes for this Mutant. Returns a list of new
		Mutant objects that are barcoded versions of self.
		"""
		# Make a list of empty lists as long as the sequence.
		print "Generating barcodes for %s mutant %s." % (self.type, self.name)
		barcode_sites=[]
		while(len(barcode_sites)<len(self.coding_sequence)):
			barcode_sites.append([])

		native_aa_sequence = Translate(self.coding_sequence)
		for position in range(0,len(self.coding_sequence)):
			if self.PositionNearMutation(position):
				# Skipping sites that are close to the
				# mutation. This is necessary so the
				# bar codes can be applied to the WT sequence.
				continue
			codon_start = position/3*3
			for candidate in ['A','C','G','T']:
				s2=Substitute(self.coding_sequence,candidate,position)
				if native_aa_sequence == Translate(s2):
					barcode_sites[position].append(candidate)
		# Translate barcode sites into (position,[nucleotide]) tuples.
		sites=[]
		for position in range(0,len(barcode_sites)):
			if len(barcode_sites[position])>1:
				sites.append( (position,barcode_sites[position]) )
		bs = BarcodeSet(sites,self.mut_codon_coords)
		newmutants=[]
		while len(newmutants) < 5:
			bcm=copy.deepcopy(self)
			try:
				bcm.ApplyBarcode(bs.next())
			except IndexError:
				print "Warning - ran out of bar codes from pool of %d for mutant %s. Created %d bar coded mutants." % ( bs.Size(), self.name, len(newmutants) )
				print Mutant.TableHeader()
				print self.AsTabDelim()
				print self
				break
			bcm.ApplyFilters()
			if bcm.passes_filters:
				newmutants.append(bcm)
			else:
				print "Warning: bar coded mutant failed filters."
		return newmutants
	
	def ApplyFilters(self):
		self.filter_results[Mutant.HasNoMonoSixmers]=self.NoMonoSixMers()
		self.filter_results[Mutant.HasNoGFourMers]=self.NoGFourMers()
		self.filter_results[Mutant.HasOKGcContent]=self.GCContentNotExtreme()
		self.passes_filters = all(self.filter_results)
	
	def ApplyBarcode(self,barcode):
		"""
		This applies the bar code to self's coding sequence.
		The bar code is a list of (position,nucleotide) tuples
		that don't change the amino acid translation.
		"""
		for ( position, nucleotide ) in barcode:
			self.coding_sequence = Substitute(self.coding_sequence, nucleotide, position )
		self.barcode = FormatBarcode(barcode)
		#self.aa_sequence = Translate(self.coding_sequence)
		assert self.aa_sequence == Translate(self.coding_sequence)

	def AsFasta(self):
		"""
		AsFasta returns mutant object in Fasta format.
		"""
		return self.FastaHeader() + '\n' + FastaSequence(self.overhang_seq_5 + self.coding_sequence + self.overhang_seq_3)

	def AsTuple(self):
		"""
		AsTuple returns the values of the mutant object as a tuple.
		"""
		return ( self.screen,
			self.fragment,
			self.type,
			self.name,
			self.barcode,
			self.aa_sequence,
			'%s' % self.mut_location,
			self.ref_nucleotide,
			self.mut_nucleotide,
			self.ref_codon,
			str(self.mut_codon),
			str(self.ref_aa),
			self.mut_aa,
			self.passes_filters
			) + tuple(self.filter_results)

	def AsTabDelim(self):
		"""
		AsTabDelim generates a single row of tab-delimited output
		for a single Mutant object.
		"""
		values = [ self.type,
			self.name,
			self.barcode,
			self.aa_sequence,
			'%s' % self.mut_location,
			self.ref_nucleotide,
			self.mut_nucleotide,
			self.ref_codon,
			str(self.mut_codon),
			str(self.ref_aa),
			self.mut_aa,
			BooleanAsIntString(self.passes_filters)
			#'%d' % int(self.passes_filters)
			]
		#values += map(str,self.filter_results)
		values += map(BooleanAsIntString,self.filter_results)

		try:
			output = '\t'.join(values) + "\n"
		except TypeError:
			print values
			raise
		return output

	def DeletionMutants(self,ifds=[]):
		"""
		Creates a list of new Mutant objects that contain every
		possible 3-base deletion. These are all in-frame, since they
		are all 3 bases long. Any deletion that creates a stop
		is discarded.

		7/13/2017 - revised to take a list of in-frame deletions.
		If this list is empty, then generate all in-frame deletions,
		otherwise only generate the ones specified in the list. These
		are in p. notation.
		"""
		# List of mutants to be returned.
		mutants=[]

		# for position in 1 ...
		for position in range(0,len(self.coding_sequence)-3):
			# create mutant.
			d = DeletionMutant(
				self.screen,
				self.fragment,
				Deletion(self.coding_sequence,position),
				self.overhang_seq_5,
				self.overhang_seq_3,
				self.first_aa_coord,
				self)
			# Exclude mutants that have a stop codon.
			if '*' in d.aa_sequence:
				continue

			d.InitMutation(position)

			# Exclude mutants that are not in ifd list, if a list
			# of ifds was given.
			if ifds and d.name not in ifds:
				continue

			mutants.append(d)
		print "Created %d in-frame deletions." % len(mutants)
		return mutants

	def InsertionMutants(self,ifis=[]):
		"""
		Creates a list of new Mutant objects each of
		which includes a new alanine at each position. For
		example given the aa sequence MSDFWYY this returns
		Mutant objects for:
			AMSDFWYY
			MASDFWYY
			MSADFWYY
			MSDAFWYY
			MSDFAWYY
			MSDFWAYY
			MSDFWYAY
			MSDFWYYA
		7/13/2017 - revised to take a list of in-frame insertions.
		If this list is empty, then generate all in-frame insertion,
		otherwise only generate the ones specified in the list. These
		are in p. notation.
		"""
		# List of mutants to be returned.
		mutants=[]

		# Get a list of codons coding for alanine, in order by
		# usage frequency.
		alanine_codons=map(lambda x:x[0],reverse_codons['A'])

		# For each in-frame position (and the +1 adds the end of the
		# coding sequence)...
		for position in range(0,len(self.coding_sequence)+1,3):
			for codon in alanine_codons:
				#print "Inserting codon %s at position %d" % (codon,position)
				new_mutant=InsertionMutant(
					self.screen,
					self.fragment,
					Insert(self.coding_sequence,codon,position),
					self.overhang_seq_5,
					self.overhang_seq_3,
					self.first_aa_coord,
					self)
				
				new_mutant.InitMutation(position)
				if '*' not in new_mutant.aa_sequence:
					mutants.append(new_mutant)
					break
		print "Created %d insertion mutants." % len(mutants)
		return mutants

	def PointMutants(self,ms=[]):
		"""
		Creates a list of new Mutant objects each of which
		has a single nucleotide change that results in an amino acid 
		change.

		7/27/2017 - revised to take a list of missense mutations.
		If this list is empty, then generate all point mutations,
		otherwise only generate the ones specified in the list. These
		are in p. notation.
		"""
		mutants=[]
		# For each nucleotide location ...
		for n in range(0,len(self.coding_sequence)):
			# ... try the 3 possible nucleotide changes there ...
			for mutant_nucleotide in ['A','C','G','T']:
				if self.coding_sequence[n] == mutant_nucleotide:
					continue
				new_nt_sequence = self.coding_sequence[0:n] + mutant_nucleotide + self.coding_sequence[n+1:]
				# ... and create the resulting Mutant.
				new_cs = PointMutant(self.screen,self.fragment,new_nt_sequence,self.overhang_seq_5,self.overhang_seq_3,self.first_aa_coord,self)
				# If the new Mutant has the same
				# amino acid sequence than self or it 
				# contains a stop codon then skip it.
				if new_cs.aa_sequence == self.aa_sequence or '*' in new_cs.aa_sequence:
					continue
				# Initialized the mutation info in the
				# new coding sequence.
				new_cs.InitMutation(n)

				# If there is a list of desired point mutations
				# and this mutant is not in it, then skip it.
				if ms and new_cs.name not in ms:
					continue

				# Add the new sequence to the list.
				mutants.append(new_cs)

		# Return the list of all mutants found.
		print "Created %d point mutants." % len(mutants)
		return mutants
	
	def NoGFourMers(self):
		"""
		Checks if coding sequence contains a GGGG sequence. If so,
		fails this test.
		"""
		if self.complete_sequence.find('G'*4)>-1:
			return False
		return True

	def NoMonoSixMers(self):
		"""
		Checks if sequence contains homopolymer > 5 bases long.
		If so, fails this test.
		"""
		for nucleotide in ['A','C','G','T']:
			if self.complete_sequence.find(nucleotide*6)>-1:
				return False
		return True
	
	def GCContentNotExtreme(self):
		"""
		Checks if nucleotide sequence has very high or very low
		GC content.
		"""
		if HighGC(self.complete_sequence) or LowGC(self.complete_sequence):
			return False
		return True
	
	def NotSplicing(self):
		"""
		Returns True if sequence does not contain a splice donor/
		acceptor site.
		"""
		pattern='GT.*AG'
		if re.search(pattern,self.complete_sequence):
			return False
		return True
	
	def __str__(self):
		spacer = ' '*len(self.overhang_seq_5)
		return "  coding_sequence: %s\n      aa_sequence: %s\ncomplete_sequence: %s\n" % ( spacer + self.coding_sequence, spacer + SpacedAa(self.aa_sequence), self.complete_sequence )

	@staticmethod
	def TableHeader():
		"""
		Returns the header of the mutant table as a tab-delimited
		list. Make sure the filter name order is the same as the
		order in which they are called!
		"""
		return '\t'.join([
			"Screen",
			"Fragment",
			"Type",
			"ProteinVariant",
			"Barcode",
			"Protein",
			"Location",
			"Ref_nuc",
			"Mut_nuc",
			"Ref_codon",
			"Mut_codon",
			"Ref_aa",
			"Mut_aa",
			"PassesTests",
			"NoMonoSixMers",
			"NoGFourMers",
			#"TransEff",
			"GCContentNotExtreme" ]) + "\n"

class PointMutant(Mutant):
	"""
	This is a Mutant that contains a point mutation.
	"""
	def __init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent):
		Mutant.__init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent)
		self.type = "point"

	def InitMutation(self,location):
		"""
		Initializes fields that store mutation info.
		"""
		# Find the start of the mutated codon.
		mutant_codon_location=location/3*3
		self.mut_codon_coords=range(mutant_codon_location,mutant_codon_location+3)

		# Reference and mutant codons.
		self.ref_codon = self.parent.coding_sequence[mutant_codon_location:mutant_codon_location+3]
		self.mut_codon = self.coding_sequence[mutant_codon_location:mutant_codon_location+3]

		# Reference and mutant amino acid.
		self.mut_aa = codons[self.mut_codon][0]
		try:
			self.ref_aa = codons[self.ref_codon][0]
		except KeyError:
			# Check if this error from an empty codon at the
			# end of the sequence. This is due to the alanine
			# scan mutants.
			if self.ref_codon == '':
				pass
			else:
				raise

		self.mut_location = location
		self.mut_aa_location = self.mut_location/3+self.first_aa_coord
		try:
			self.ref_nucleotide = self.parent.coding_sequence[location]
		except IndexError:
			# This may be due to the insertional mutants.
			if location >= len(self.parent.coding_sequence):
				self.ref_nucleotide='-'
			else:
				raise
		self.mut_nucleotide = self.coding_sequence[location]
		self.name = "p.%s%d%s" % ( self.ref_aa, self.mut_aa_location, self.mut_aa)

class InsertionMutant(Mutant):
	"""
	This is a Mutant that contains an insertional mutation.
	"""
	def __init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent):
		Mutant.__init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent)
		self.type = "insertion"

	def InitMutation(self,location):
		"""
		Initializes fields that store mutation info.
		"""
		# Find the start of the mutated codon.
		mutant_codon_location=location/3*3
		self.mut_codon_coords=range(mutant_codon_location,mutant_codon_location+3)

		# Reference and mutant codons.
		self.ref_codon = self.parent.coding_sequence[mutant_codon_location:mutant_codon_location+3]
		self.mut_codon = self.coding_sequence[mutant_codon_location:mutant_codon_location+3]

		# Reference and mutant amino acid.
		self.mut_aa = codons[self.mut_codon][0]
		try:
			self.ref_aa = codons[self.ref_codon][0]
		except KeyError:
			# Check if this error from an empty codon at the
			# end of the sequence. This is due to the alanine
			# scan mutants.
			if self.ref_codon == '':
				pass
			else:
				raise

		self.mut_location = location
		self.mut_aa_location = self.mut_location/3+self.first_aa_coord
		# Need to store name as p.Xn_Yn+1insA notation to be consistent
		# with ifi list.
		#self.name = "p.%s%dinsA" % ( self.ref_aa, self.mut_aa_location )
		self.name = "p.%s%dinsA" % ( self.ref_aa, self.mut_aa_location )
		try:
			self.ref_nucleotide = self.parent.coding_sequence[location]
		except IndexError:
			# This may be due to the insertional mutants.
			if location >= len(self.parent.coding_sequence):
				self.ref_nucleotide='-'
			else:
				raise
		self.mut_nucleotide = self.coding_sequence[location]

	def IsEfficientTranslator(self):
		"""
		Checks if this mutant's codon has a frequency as good or
		better than the wild type codon.
		"""
		try:
			ref_codon_frequency = codons[self.ref_codon][1]
		except KeyError:
			return True
		mut_codon_frequency = codons[self.mut_codon][1]

		return mut_codon_frequency >= ref_codon_frequency

	def AdjustedBarcode(self):
		"""
		Returns the parsed bar code. If this mutant contains an
		insertion or deletion then the coordinates following the
		mutation are adjusted. The intent is that this bar code
		can be applied to the WT sequence.
		"""
		barcode_list = ParseBarcode(self.barcode)
		# Adjust the locations of nucleotides following the 
		# insertion. Need to subtract 3 from these locations
		# so they map to the corresponding location in the WT.
		adjusted_barcode = []
		for (location, nucleotide) in barcode_list:
			if location >= self.mut_codon_coords[0]:
				adjusted_barcode.append((location-3,nucleotide))
			else:
				adjusted_barcode.append((location,nucleotide))
		return adjusted_barcode

class DeletionMutant(Mutant):
	"""
	This is a Mutant that contains a deletional mutation.
	"""
	def __init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent):
		Mutant.__init__(self,screen,fragment,coding_sequence,overhang_seq_5,overhang_seq_3,first_aa_coord,parent)
		self.type = "deletion"

	def IsEfficientTranslator(self):
		return True

	def InitMutation(self,location):
		self.mut_location = location
		self.ref_nucleotide = self.parent.coding_sequence[location]
		self.mut_nucleotide = self.coding_sequence[location]
		mutant_codon_location=location/3*3

		# Reference and mutant codons.
		self.ref_codon = self.parent.coding_sequence[mutant_codon_location:mutant_codon_location+3]
		self.mut_codon = self.coding_sequence[mutant_codon_location:mutant_codon_location+3]

		self.ref_aa = codons[self.ref_codon][0]
		self.mut_aa_location = self.mut_location/3+self.first_aa_coord
		mut_aa_index=self.mut_location/3

		# Name of the deletion is tricky. See 
		# http://www.hgvs.org/mutnomen/recs-prot.html#del .
		# If frame 1 deletion, name like p.Q3del
		# If not frame 1, name like p.Q3_L4delinsE, ie the deletion
		# takes out Q3 and L4, and leaves an E in its place.
		self.mut_aa = self.aa_sequence[mut_aa_index]
		# If frame 1...
		if mutant_codon_location == location:
			self.name = 'p.%s%ddel' % ( self.ref_aa, self.mut_aa_location)
			# No need to exclude any site from bar code creation.
			self.mut_codon_coords=[]
		else:
			next_aa = self.parent.aa_sequence[mut_aa_index+1]
			self.name = 'p.%s%d_%s%ddelins%s' % ( self.ref_aa, self.mut_aa_location, next_aa, self.mut_aa_location + 1, self.mut_aa )
			# Exclude a wide area around the deletion for barcode
			# creation. Just to be safe.
			self.mut_codon_coords=range(mutant_codon_location,mutant_codon_location+6)

	def AdjustedBarcode(self):
		"""
		Returns the parsed bar code. If this mutant contains an
		insertion or deletion then the coordinates following the
		mutation are adjusted. The intent is that this bar code
		can be applied to the WT sequence.
		"""
		barcode_list = ParseBarcode(self.barcode)
		# Adjust the locations of nucleotides following the 
		# deletion. Need to add 3 to these locations
		# so they map to the corresponding location in the WT.
		adjusted_barcode = []
		for (location, nucleotide) in barcode_list:
			if location >= self.mut_location:
				adjusted_barcode.append((location+3,nucleotide))
			else:
				adjusted_barcode.append((location,nucleotide))
		return adjusted_barcode

if __name__ == "__main__":
	# Create a Mutant, print it, and its point mutants.
	#orig = Mutant('CCAGGT','AAAAA','GGGGG')
	#print "Original:"
	#print orig
	#print "Mutants:"
	#for mutant in orig.PointMutants():
	#	print mutant
	print FormatBarcode([(2,'A'),(3,'C'),(17,'G')])
