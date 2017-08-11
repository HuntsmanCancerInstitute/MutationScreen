import math
import random

def ParseBarcode(barcode):
	"""
	This is the inverse function on FormatBarcode, ie it takes a
	barcode string (e.g "2C_5G_17T") and returns a list of
	(location,nucleotide) tuples, such as [(2,'C'),(5,'G'),(17,'T')].
	"""
	# This lambda expression takes a string like '14G' and returns
	# a tuple of (14,'G').
	return map(lambda x:(int(x[:-1]),x[-1]),barcode.split('_'))

def FormatBarcode(barcode):
	"""
	barcode is a list of (location,nucleotide) tuples. FormatBarcode
	returns the barcode as an underscore-delimited string of location
	nucleotide pairs. For example 2C_5G_17T.
	"""
	return '_'.join(map(lambda x:'%d%s' % x, barcode))

def Combine(l1,l2):
	retval=[]
	for item1 in l1:
		for item2 in l2:
			retval.append(item1+item2)
	return retval

def MakeBarcodes(l):
	"""
	"""
	(location,choices)=l[0]
	barcodes=map(lambda x:[(location,x)],choices)
	if len(l)==1:
		return barcodes
	else:
		return Combine(barcodes,MakeBarcodes(l[1:]))

class BarcodeSet:
	"""
	BarcodeSet holds all the bar codes possible given the variant
	locations and the possible nucleotides at those locations.

	This implementation looks for sites where any nucleotide is 
	allowed, and generates combinations of them.

	A BarcodeSet is an iterable object, so one could call it as:
		for barcode in BarcodeSet(sites):
			do something ...
	"""
	max_barcode_length = 8
	barcodes_distributed = set()

	def __init__(self,sites,excluded_locations=[]):
		#print "Creating BarcodeSet using sites", sites

		# Find locations with the most options.
		l = map(lambda x:(len(x[1]),x[0],x[1]), sites)
		l.sort()
		l.reverse()
		#print l
		numoptions=l[0][0]
		#print "Choosing sites with %d options." % numoptions
		chosensites=[]
		for site in l:
			if site[0] == numoptions and site[1] not in excluded_locations:
				chosensites.append((site[1],site[2]))
			if len(chosensites) >= BarcodeSet.max_barcode_length:
				break
		chosensites.sort()
		#print "These sites have %d options:" % numoptions
		#print chosensites

		# Generate bar codes.
		#print "Found %d sites with %d options each." % (len(chosensites),numoptions)
		#print "Generating %d bar codes." % int(math.pow(numoptions,len(chosensites)))
		self.barcodes = MakeBarcodes(chosensites)
		self.set_size = len(self.barcodes)
		random.shuffle(self.barcodes)
		#print "Created %d bar codes." % len(self.set_size)
		
	def __iter__(self):
		return self

	def next(self):
		"""
		next() returns another bar code from the BarcodeSet. This
		method checks the barcodes_distributed set to make sure
		no two bar codes are handed out during the same run.
		"""
		next_barcode = self.barcodes.pop()
		while FormatBarcode(next_barcode) in BarcodeSet.barcodes_distributed:
			next_barcode = self.barcodes.pop()
		BarcodeSet.barcodes_distributed.add(FormatBarcode(next_barcode))
		return next_barcode
	
	def Size(self):
		return self.set_size


if __name__ == "__main__":
	sites=[(2, ['C', 'T']), (3, ['C', 'T']), (5, ['A', 'G']), (8, ['A', 'C', 'G', 'T']), (11, ['A', 'C', 'G', 'T']), (14, ['A', 'C', 'G', 'T']), (17, ['A', 'C', 'G', 'T']), (20, ['A', 'C', 'G', 'T']), (23, ['A', 'G']), (26, ['A', 'G']), (29, ['A', 'C', 'G', 'T']), (32, ['A', 'G']), (35, ['C', 'T']), (38, ['A', 'C', 'G', 'T']), (41, ['A', 'C', 'T']), (44, ['C', 'T']), (47, ['A', 'C', 'G', 'T']), (53, ['A', 'G'])]
	#s1=MakeBarcodes([(2,['A','C'])])
	#s2=MakeBarcodes([(3,['G','T'])])
	#print s1
	#print s2
	#print Combine(s1,s2)
	bs=BarcodeSet(sites)

	# Get next 5 random bar codes.
	for i in range(0,5):
		print bs.next()
