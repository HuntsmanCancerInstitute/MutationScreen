# Code to id missense substitution in a protein sequence.

s1="ASDFASDFASDFASDF"
s2="ASDFASDGASDFASDF"
s3="MSDFASDFASDFASDF"
s4="ASDFASDFASDFASDW"

def FindMSSub(s1,s2):
	"""
	Finds first missense substitution.
	"""
	l=zip(list(s1),list(s2))
	for location in range(0,len(l)):
		if l[location][0] != l[location][1]:
			return (location,l[location][1])
	return (None,None)

#print FindMSSub(s1,s1)
#print FindMSSub(s1,s2)
#print FindMSSub(s1,s3)
#print FindMSSub(s1,s4)

