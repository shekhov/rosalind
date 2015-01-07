""" Script for only Rosalind tasks """
import aa
import dna
import tools

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
	RNA_2 = dna.dnaToRna(tools.reverseSequence(DNA))
	peptides = aa.translation(RNA) +aa.translation(RNA_2)
	return peptides