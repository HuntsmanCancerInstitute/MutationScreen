def MissenseSubstitute(aa_sequence,new_aa,location):
	return aa_sequence[0:location]+new_aa+aa_sequence[location+1:]

amino_acids=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y',]

class Protein:
	"""
	Protein represents an amino acid sequence.
	"""
	def __init__(self,aa_sequence):
		self.aa_sequence=aa_sequence
	
	def GenerateMissenseSubstitutions(self):
		"""
		Returns list of all possible missense substitutions from the
		original sequence. The number of these should be 
		19*len(orig_sequence).
		"""
		sequences=[]
		for position in range(0,len(self.aa_sequence)):
			for aa in amino_acids:
				if aa != self.aa_sequence[position]:
					sequences.append(MissenseSubstitute(self.aa_sequence,aa,position))
		sequences.sort()
		return sequences
	
	def IdMissenseSubstitution(self,other_aa_sequence):
		"""
		Finds first missense substitution.
		"""
		l=zip(list(self.aa_sequence),list(other_aa_sequence))
		for location in range(0,len(l)):
			if l[location][0] != l[location][1]:
				return (location,l[location][1])
		return (None,None)

if __name__ == "__main__":
	p=Protein("ASDF")
	print p.GenerateMissenseSubstitutions()
