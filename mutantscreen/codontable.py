codons= {
	'TTT':('F',0.46),
	'TTC':('F',0.54),
	'TTA':('L',0.08),
	'TTG':('L',0.13),
	'CTT':('L',0.13),
	'CTC':('L',0.20),
	'CTA':('L',0.07),
	'CTG':('L',0.40),
	'ATT':('I',0.36),
	'ATC':('I',0.47),
	'ATA':('I',0.17),
	'ATG':('M',1.00),
	'GTT':('V',0.18),
	'GTC':('V',0.24),
	'GTA':('V',0.12),
	'GTG':('V',0.46),
	'TCT':('S',0.19),
	'TCC':('S',0.22),
	'TCA':('S',0.15),
	'TCG':('S',0.05),
	'CCT':('P',0.29),
	'CCC':('P',0.32),
	'CCA':('P',0.28),
	'CCG':('P',0.11),
	'ACT':('T',0.25),
	'ACC':('T',0.36),
	'ACA':('T',0.28),
	'ACG':('T',0.11),
	'GCT':('A',0.27),
	'GCC':('A',0.40),
	'GCA':('A',0.23),
	'GCG':('A',0.11),
	'TAT':('Y',0.44),
	'TAC':('Y',0.56),
	'TAA':('*',0.30),
	'TAG':('*',0.24),
	'CAT':('H',0.42),
	'CAC':('H',0.58),
	'CAA':('Q',0.27),
	'CAG':('Q',0.73),
	'AAT':('N',0.47),
	'AAC':('N',0.53),
	'AAA':('K',0.43),
	'AAG':('K',0.57),
	'GAT':('D',0.46),
	'GAC':('D',0.54),
	'GAA':('E',0.42),
	'GAG':('E',0.58),
	'TGT':('C',0.46),
	'TGC':('C',0.54),
	'TGA':('*',0.47),
	'TGG':('W',1.00),
	'CGT':('R',0.08),
	'CGC':('R',0.18),
	'CGA':('R',0.11),
	'CGG':('R',0.20),
	'AGT':('S',0.15),
	'AGC':('S',0.24),
	'AGA':('R',0.21),
	'AGG':('R',0.21),
	'GGT':('G',0.16),
	'GGC':('G',0.34),
	'GGA':('G',0.25),
	'GGG':('G',0.25),
}

# reverse_codons is a dictionary indexed by amino acid of lists of codons
# and codon frequency (in human) for that amino acid. Codons for each
# amino acid are in descending frequency order.
reverse_codons={
	'*':[ ('TGA',0.47), ('TAA',0.30), ('TAG',0.24), ],
	'A':[ ('GCC',0.40), ('GCT',0.27), ('GCA',0.23), ('GCG',0.11), ],
	'C':[ ('TGC',0.54), ('TGT',0.46), ],
	'D':[ ('GAC',0.54), ('GAT',0.46), ],
	'E':[ ('GAG',0.58), ('GAA',0.42), ],
	'F':[ ('TTC',0.54), ('TTT',0.46), ],
	'G':[ ('GGC',0.34), ('GGA',0.25), ('GGG',0.25), ('GGT',0.16), ],
	'H':[ ('CAC',0.58), ('CAT',0.42), ],
	'I':[ ('ATC',0.47), ('ATT',0.36), ('ATA',0.17), ],
	'K':[ ('AAG',0.57), ('AAA',0.43), ],
	'L':[ ('CTG',0.40), ('CTC',0.20), ('CTT',0.13), ('TTG',0.13), ('TTA',0.08), ('CTA',0.07), ],
	'M':[ ('ATG',1.00), ],
	'N':[ ('AAC',0.53), ('AAT',0.47), ],
	'P':[ ('CCC',0.32), ('CCT',0.29), ('CCA',0.28), ('CCG',0.11), ],
	'Q':[ ('CAG',0.73), ('CAA',0.27), ],
	'R':[ ('AGA',0.21), ('AGG',0.21), ('CGG',0.20), ('CGC',0.18), ('CGA',0.11), ('CGT',0.08), ],
	'S':[ ('AGC',0.24), ('TCC',0.22), ('TCT',0.19), ('AGT',0.15), ('TCA',0.15), ('TCG',0.05), ],
	'T':[ ('ACC',0.36), ('ACA',0.28), ('ACT',0.25), ('ACG',0.11), ],
	'V':[ ('GTG',0.46), ('GTC',0.24), ('GTT',0.18), ('GTA',0.12), ],
	'W':[ ('TGG',1.00), ],
	'Y':[ ('TAC',0.56), ('TAT',0.44), ],
}

def MakeCodons(s,codon_list=[]):
	"""
	MakeCodons transforms nucleotide string into a list of codons. 
	If string length is not a multiple of 3 then any extra nucleotides 
	are ignored.
	"""
	if len(s)<3:
		# We're done.
		return codon_list
	else:
		return MakeCodons(s[3:],codon_list+[s[0:3]])

def Translate(nt_sequence):
	"""
	Translates a nucleotide sequence into amino acids.
	"""
	return ''.join(map(lambda x:codons[x][0], MakeCodons(nt_sequence)))

def MakePossibleSNSCodons(codon):
	"""
	MakePossibleSNSCodons(codon) returns a list of the codons resulting
	from all possible single nucleotide substitutions to a given codon.
	"""
	l=[]
	s=set(list("ACGT"))	# Set containing all 4 nucleotides.
	# First base changes
	for n in s.difference(set(codon[0])):
		l.append(n+codon[1:])
	# Second base changes
	for n in s.difference(set(codon[1])):
		l.append(codon[0]+n+codon[2])
	# Third base changes
	for n in s.difference(set(codon[2])):
		l.append(codon[0:2]+n)

	return l
		
def PossibleMissenseSubstitutions(codon):
	"""
	Given a codon, return the codons that result from single nucleotide
	substitutions as long as they create a missense substitution.
	Return these as a dictionary of codon lists indexed by amino acid.
	"""
	d={}
	# This is what the codon codes for.
	orig_aa = codons[codon][0]

	# Generate all the possible codons that are single-nucleotide 
	# substitutions. Screen those codons to make sure they result
	# in a missense substitution. If they pass, store them in a
	# dictionary.
	for candidate_codon in MakePossibleSNSCodons(codon):
		new_aa = codons[candidate_codon][0]
		if new_aa != orig_aa and new_aa != '*':
			try:
				d[new_aa].append(candidate_codon)
			except KeyError:
				d[new_aa] = [candidate_codon]
	return d

if __name__ == "__main__":
	codon="GGG"
	possible_codons=MakePossibleSNSCodons(codon)
	print "%d possible SNS codons: %s" % ( len(possible_codons), possible_codons )

	d=PossibleMissenseSubstitutions(codon)
	print "Possible amino acid changes:", d.keys()
	print d
