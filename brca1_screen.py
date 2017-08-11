#!/usr/bin/env python
import sys
import mutantscreen
from mutantscreen.mutantset import MutantSet

def main():
	sys.stdout.write("%s version %s\n" % ( sys.argv[0], mutantscreen.__version__))

	brca1_sequence="""gttgacattgattattgactagttattaatagtaatcaattacggggtcattagttcatagcccatatatggagttccgcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatggtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagctctctggctaactagagaacccactgcttactggcttatcgaaattaatacgactcactatagggagacccaagctggctagccaccatggatttatctgctcttcgcgttgaagaagtacaaaatgtcattaatgctatgcagaaaatcttagagtgtcccatctgtctggagttgatcaaggaacctgtctccacaaagtgtgaccacatattttgcaaattttgcatgctgaaacttctcaaccagaagaaagggccttcacagtgtcctttatgtaagaatgatataaccaaaaggagcctacaagaaagtacgagatttagtcaacttgttgaagagctattgaaaatcatttgtgcttttcagcttgacacaggtttggagtatgcaaacagctataattttgcaaaaaaggaaaataactctcctgaacatctaaaagatgaagtttctatcatccaaagtatgggctacagaaaccgtgccaaaagacttctacagagtgaacccgaaaatccttccttgcaggaaaccagtctcagtgtccaactctctaaccttggaactgtgagaactctgaggacaaagcagcggatacaacctcaaaagacgtctgtctacattgaattgggatctggaggttctggcggtggatctggtcgggctgactacaaagaccatgacggtgattataaagatcatgacatcgactacaaggatgacgatgacaagtctccagcgattccgtcgacaccacctactccctctccaggcggatccaagctactgtcttctatcgaacaagcatgcgatatttgccgacttaaaaagctcaagtgctccaaagaaaaaccgaagtgcgccaagtgtctgaagaacaactgggagtgtcgctactctcccaaaaccaaaaggtctccgctgactagggcacatctgacagaagtggaatcaaggctagaaagactggaacagctatttctactgatttttcctcgagaagaccttgacatgattttgaaaatggattctttacaggatataaaagcattgttaacaggattatttgtacaagataatgtgaataaagatgccgtcacagatagattggcttcagtggagactgatatgcctctaacattgagacagcatagaataagtgcgacatcatcatcggaagagagtagtaacaaaggtcaaagacagttgactgtatgataagcggccgctcgagtctagagggcccgtttaaacccgctgatcagcctcgactgtgccttctagttgccagccatctgttgtttgcccctcccccgtgccttccttgaccctggaaggtgccactcccactgtcctttcctaataaaatgaggaaattgcatcgcattgtctgagtaggtgtcattctattctggggggtggggtggggcaggacagcaagggggaggattgggaagacaatagcaggcatgctggggatgcggtgggctctatgg""".upper()
	
	fragments = [ MutantSet( "BRCA1",
			"Fragment1",
			brca1_sequence,
			synth_from=658,
			synth_to=799,		# Synthesizing 142 nt
			mut_from=677,
			mut_to=775,		# Mutating 99 nt, 33 aa
			first_aa_coord=2),
		MutantSet( "BRCA1",
			"Fragment2",
			brca1_sequence,
			synth_from=753,
			synth_to=900,		# Synthesizing 148 nt
			mut_from=776,
			mut_to=874,
			first_aa_coord=35),		# Mutating 99 nt, 33 aa
		MutantSet( "BRCA1",
			"Fragment3",
			brca1_sequence,
			synth_from=846,
			synth_to=999,		# Synthesizing 154 nt.
			mut_from=875,
			mut_to=973,
			first_aa_coord=68), ]		# Mutating 99 nt, 33 aa

	# Nucleotide 1 of each mutation range must be the first 
	# nucleotide of a codon.
	# base sequence ...----------------------------------------------...
	# synth sequence        --------------------------------------
	# mutate seq                   ---------------------------
	# overhang seq          -------
	# The synth sequence may not be in frame because the overhang may
	# not be a multiple of 3 in length, but the mutate sequence must
	# be in frame.
	
	for fragment in fragments:

		print fragment

		# Generate mutants.
		fragment.GenerateMutants()

		# Filter mutants.
		fragment.FilterMutants()

		# Generate bar codes.
		fragment.CreateBarcodes()

		# Write them to a file.
		fragment.Write(mutantscreen.__version__)

		# Write them to a Fasta file.
		fragment.AsFasta()
		print

		# stop after the first fragment for testing.
		#break

if __name__ == "__main__":
	main()
