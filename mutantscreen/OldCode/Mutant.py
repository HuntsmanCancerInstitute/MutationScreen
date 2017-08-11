from CodonTable import codons, reverse_codons

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

class Mutant:
	"""
	Mutant represents a specific nucleotide change to a single
	location in the parent MutantSet's sequence.

	The Mutant class contains a number of filter methods. These should
	return True if the mutant is good, False if the mutant should be
	discarded.
	"""

	def __init__(self,mutantset,location,mutant_nucleotide,ref_nucleotide):
		# Parent mutant set.
		self.mutantset=mutantset

		# Position in the synthesized sequence, 0-based.
		self.location=location

		# Mutant nucleotide
		self.mutant_nucleotide=mutant_nucleotide

		# Reference (original) nucleotide.
		self.ref_nucleotide=ref_nucleotide

		# Store the mutated synthesized sequence, since we need
		# this for several tests.
		self.mutated_sequence = MutateSequence(self.mutantset.sequence,mutant_nucleotide,location)
		
		# Codon_start_location - 0-based coordinate of start
		# of codon containing mutant.
		self.codon_start_location = (self.location/3)*3

		# List of bar codes for this mutant.
		self.barcodes = []

		# Pass_filters: whether mutant passes all filters.
		# Defaults to True until mutant fails a filter.
		self.passes_filters = True

		self.ref_codon = self.mutantset.sequence[self.codon_start_location:self.codon_start_location+3]
		self.mut_codon = self.mutated_sequence[self.codon_start_location:self.codon_start_location+3]
		self.ref_aa = codons[self.ref_codon][0]
		self.mut_aa = codons[self.mut_codon][0]

		self.filter_results=[]
	
	def NotSilent(self):
		"""
		Returns True if this is a silent mutation.
		"""
		return self.ref_aa != self.mut_aa

	def NotSplicing(self):
		"""
		Returns True if Mutant does not contain a splice donor/
		acceptor site.
		"""
		pattern='GT.*AG'
		if re.search(pattern,self.mutated_sequence):
			return False
		return True
	
	def NoGFourMers(self):
		"""
		Checks if mutation creates a GGGG sequence. If so,
		mutant fails this test.
		"""
		if self.mutated_sequence.find('G'*4)>-1:
			return False
		return True

	def NoSixMers(self):
		"""
		Checks if mutation creates homopolymer > 5 bases in
		the sequence. If so, mutant fails this test.
		"""
		for nucleotide in ['A','C','G','T']:
			if self.mutated_sequence.find(nucleotide*6)>-1:
				return False
		return True
	
	def GCContentNotExtreme(self):
		if HighGC(self.mutated_sequence) or LowGC(self.mutated_sequence):
			return False
		return True
		
	def IsEfficientTranslator(self):
		"""
		Checks if this mutant's codon has the greatest frequency
		of any codon the this mutant's amino acid.
		"""
		best_mutant_codon = reverse_codons[self.mut_aa][0][0]
		return best_mutant_codon == self.mut_codon
	
	def RNASecondaryStructureOK(self):
		"""
		"""
		pass
	
	def GenerateBarcodes(self):
		"""
		Generates bar codes for this Mutant.
		"""
		pass
	
	def ApplyFilters(self):
		self.filter_results.append(self.NotSilent())
		self.filter_results.append(self.NoSixMers())
		self.filter_results.append(self.NoGFourMers())
		self.filter_results.append(self.IsEfficientTranslator())
		self.filter_results.append(self.GCContentNotExtreme())
		self.passes_filters = all(self.filter_results)

	def AsTabDelim(self):
		"""
		ToTable generates a single row of tab-delimited output
		for a single Mutant object.
		"""
		values = [ '%d' % self.location,
			self.ref_nucleotide,
			self.mutant_nucleotide,
			self.ref_codon,
			self.mut_codon,
			self.ref_aa,
			self.mut_aa,
			'%s' % self.passes_filters
			]
		values += map(str,self.filter_results)

		return '\t'.join(values) + "\n"
	
	@staticmethod
	def TableHeader():
		"""
		Returns the header of the mutant table as a tab-delimited
		list. Make sure the filter name order is the same as the
		order in which they are called!
		"""
		return '\t'.join([
			"Location",
			"Ref_nuc",
			"Mut_nuc",
			"Ref_codon",
			"Mut_codon",
			"Ref_aa",
			"Mut_aa",
			"PassesTests",
			"NotSilent",
			"NoSixMers",
			"NoGFourMers",
			"TransEff",
			"GCContentNotExtreme" ]) + "\n"

