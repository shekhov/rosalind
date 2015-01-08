""" Script for only Rosalind tasks """
import aa
import dna
import tools
from collections import Counter

def ORF (fasta_file_loc):
	""" 
	Return all open reading frames from the given DNA string.
	Open reading frame starting from start codon, and finished as soon as the stop codon appears
	http://rosalind.info/problems/orf/
	run: print (convert_array_to_string(ORF('test_files/rosalind_splc.txt'), "\n"))
	"""
	
	DNA = tools.fasta_to_sequence(tools.file_to_string(fasta_file_loc))
	dna.isNucleotide(DNA)
	
	RNA = dna.dnaToRna(DNA)
	RNA_2 = dna.dnaToRna(dna.reverse_dna(DNA))
	
	temp_result = aa.translation(RNA) + aa.translation(RNA_2)
	peptides = []
	for p in Counter(temp_result):
		peptides.append(p)
	return peptides
	
def FSH (fasta_file_loc):
	"""
		Returns longest shared motif from all given DNA sequences 
		http://rosalind.info/problems/lcsm/
		run:
	"""
	result = ''
	DNA = tools.fasta_to_collection (tools.file_to_string (fasta_file_loc))
	
	return result