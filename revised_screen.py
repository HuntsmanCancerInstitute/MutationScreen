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
			synth_from=654,
			synth_to=798,		# Synthesizing 142 nt
			mut_from=677,
			mut_to=775,		# Mutating 99 nt, 33 aa
			first_aa_coord=2,
			num_barcodes=1,
			ms=[
				'p.D2Y',
				'p.L3S',
				'p.S4F',
				'p.A5D',
				'p.L6R',
				'p.R7C',
				'p.V8D',
				'p.E9V',
				'p.E10V',
				'p.V11E',
				'p.Q12L',
				'p.N13I',
				'p.V14D',
				'p.I15T',
				'p.I15N',
				'p.N16I',
				'p.N16Y',
				'p.A17D',
				'p.M18K',
				'p.Q19P',
				'p.K20I',
				'p.I21S',
				'p.L22S',
				'p.E23V',
				'p.C24F',
				'p.P25L',
				'p.I26N',
				'p.C27W',
				'p.L28Q',
				'p.E29V',
				'p.L30S',
				'p.I31N',
				'p.K32M',
				'p.E33V',
				'p.P34R',
				],
			ifi=[ 'p.M1_D2insA', 'p.Q19_K20A', 'p.I26_C27insA', ],
			ifd=[ 'p.D2del', 'p.L3del', 'p.S4del', 'p.A5del',
			'p.Q12del', 'p.Q19del', 'p.I21del', 'p.C27del', ],
			ms_barcode_choices=[('c.21C>T','c.96G>A'),],
			id_barcode_choices=[('c.42C>T','c.75C>A'),],
			),
		MutantSet( "BRCA1",
			"Fragment2",
			brca1_sequence,
			synth_from=753,
			synth_to=901,		# Synthesizing 148 nt
			mut_from=776,
			mut_to=874,		# Mutating 99 nt, 33 aa
			first_aa_coord=35,
			num_barcodes=4,
			ms=[],			# Generate all MS.
			ifi=[
				'p.S36_T37insA',
				'p.F43_C44insA',
				'p.L52_N53insA',
				'p.N53_Q54insA',
				'p.Q54_K55insA',
				'p.K55_K56insA',
				'p.K56_G57insA',
				'p.S59_Q60insA',
				'p.Q60_C61insA',
				'p.K65_N66insA',
			],
			ifd=[	
				'p.V35del',
				'p.S36del',
				'p.T37del',
				'p.K38del',
				'p.C39del',
				'p.D40del',
				'p.H41del',
				'p.I42del',
				'p.F43del',
				'p.C44del',
				'p.K45del',
				'p.F46del',
				'p.C47del',
				'p.M48del',
				'p.L49del',
				'p.K50del',
				'p.L51del',
				'p.L52del',
				'p.N53del',
				'p.Q54del',
				'p.K55del',
				'p.K56del',
				'p.G57del',
				'p.P58del',
				'p.S59del',
				'p.Q60del',
				'p.C61del',
				'p.P62del',
				'p.L63del',
				'p.C64del',
				'p.K65del',
				'p.N66del',
				'p.D67del',
			],
			ms_barcode_choices=[
				('c.105C>T','c.120C>T'),
				('c.105C>T','c.156C>T'),
				('c.105C>T','c.180G>A'),
				('c.105C>T','c.189A>G'),
				],
			id_barcode_choices=[
				('c.108C>G','c.114G>A'),
				('c.108C>G','c.141C>T'),
				('c.108C>G','c.171G>C'),
				('c.108C>G','c.174T>A'),
				],
			),		
		MutantSet( "BRCA1",
			"Fragment3",
			brca1_sequence,
			synth_from=852,
			synth_to=996,		# Synthesizing 154 nt.
			mut_from=875,
			mut_to=973,		# Mutating 99 nt, 33 aa
			first_aa_coord=68,
			num_barcodes=4,
			ms=[
				'p.I68K',
				'p.T69I',
				'p.K70I',
				'p.R71G',
				'p.S72R',
				'p.L73Q',
				'p.Q74L',
				'p.E75V',
				'p.S76I',
				'p.T77M',
				'p.R78I',
				'p.F79S',
				'p.S80C',
				'p.Q81L',
				'p.L82R',
				'p.V83D',
				'p.E84V',
				'p.E85K',
				'p.L86P',
				'p.L87S',
				'p.K88I',
				'p.I89N',
				'p.I90N',
				'p.C91W',
				'p.A92D',
				'p.F93C',
				'p.Q94L',
				'p.L95P',
				'p.D96Y',
				'p.T97I',
				'p.G98C',
				'p.L99S',
				'p.E100V',
				'p.R78G',
				'p.S80I',
				'p.E85V',
				'p.L86Q',
				'p.L86R',
				'p.C91G',
				'p.L95R',
				],
			ifi=[],
			ifd=[]), 
		]		

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
