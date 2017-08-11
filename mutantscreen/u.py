from CodonTable import Translate
from Mutant import Substitute

s="GATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAG"
print s

# Make a list of empty lists as long as the sequence.
barcode_sites=[]
while(len(barcode_sites)<len(s)):
	barcode_sites.append([])
#print barcode_sites
nucleotides=['A','C','G','T']
for position in range(0,len(s)):
	codon_start = position/3*3
	for candidate in nucleotides:
		#if candidate == s[position]:
		#	continue
		s2=Substitute(s,candidate,position)
		#if Translate(s[codon_start:codon_start+3]) == Translate(s2[codon_start:codon_start+3]):
		if Translate(s) == Translate(s2):
			barcode_sites[position].append(candidate)
		#print position
		#print candidate
		#print s
		#print s2
		#print Translate(s)
		#print Translate(s2)
		#print barcode_sites
# Translate barcode sites into (position,[nucleotide]) tuples.
sites=[]
for position in range(0,len(barcode_sites)):
	if len(barcode_sites[position])>1:
		sites.append( (position,barcode_sites[position]) )
print sites

# Count locations with 0, 1, 2, or 3 bases.
d={0:0,1:0,2:0,3:0,4:0}
for location in barcode_sites:
	d[len(location)]+=1
print d
