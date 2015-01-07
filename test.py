import unittest
import dna
import aa
import tools


NUCLEOTIDES = "ACTGU"
TEST_FASTA_LOC = 'fasta_test.txt'
TEST_FASTA_STRING = """>Rosalind_1497
GAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATA
AGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG
AAGCGGCACCTAGTCGAGCCTCCCACCCTTATCTTGAACTGAACCAGAGTCCATGGCCGT
GGTCGTTAGGGGACTATGGAGACCGCTTCATGGGGAGGGAAAACGTTTTGCTATGTACTG
GGTGCTGTTTTCGGAGTCGTGGGAAGAGACATTTACAAACGTGGCAAAAGCTCTCATGTG
GCTCCATGTAGGCGTCCAGTTCTCGACATTTCAGTACTCAACGGTCTGCTTCGTTCCGGG
CAATCTCTTTAGACCCAGCGGGGTCTTATCATCGCGTTTGAAGACATTCCGTGGAAGGGG
GCGGTCTCGGGAGGATATCGCCTGGGATGGTCTAGCCCGGAGACGTATCTCTGTAGCTAC
AGAGATACGTCTCCGGGCTAGACCATCGCCCGTTCTAACGATACCCGCGAAGGAAATTAA
GTACGAGCTTCTGGATATAGATGTTCGTTTCTTCTTATCCCTATGACGTCTCTCATGTAA
CCGATACTGTTTCCGCTATTGATGTTTACCTGGATAAAGTTTCGTGCCTAGACTTTTAGG
GGACTCGGCTCCTGTCCCCACAGCCGGTATGAGAAGAGCTCATAAACTCGACAACATCAA
TGATTCATCTAACTGGGTATTATCAGGTCTGCGTAACCTTTGAGCAATTCGCATGTAGGT
CTATGGTCGCTCTAGTCCACACGCGCTATGCGCAACACTAGCCCTCGCGGTGTTGCCATT
GACCGCCTCCAACGCATCCTCATCTGTCCCCGCTCCTGCATTCCGGAACCCAACAAGCTT
TCAAGATGCATTTTAGCTCAGCTGCTGGGCGATAAACGTCTCTCATTGAGTG"""

class CountNucleotidesTest (unittest.TestCase):

	def testArray(self):
		""" count_nucleotides should return an array with 4 numbers """
		self.assertEqual (4, len(dna.count_nucleotides("")))
	
	def testResult(self):
		""" count_nucleotides should return right amount of nucleotides """
		sequence = "AAATTTTGGGGCCCCC"
		result = [3, 4, 4, 5]
		self.assertEqual (result, dna.count_nucleotides(sequence))

class NucleotideStringTestCase:
	# All methods working with sequences of DNA and RNA should be tested here
	test_method = None
	def testUpperCase (self):
		""" All letters should be in uppercase """
		self.assertRaises (dna.InvalidCharacterError, self.test_method, "acgc")
			
	def testNucleotides (self):
		""" count_nucleotides should work only with 5 letters: A T G C U """
		self.assertRaises (dna.InvalidSequenceError, self.test_method, 'R')

	def testT_or_U (self):
		""" The sequence should not contain both Thymine and Uracil """
		self.assertRaises (dna.InvalidSequenceError, self.test_method, NUCLEOTIDES)
	
		
		
class DNAStringCheckFailed (unittest.TestCase, NucleotideStringTestCase):
	def setUp (self):	
		self.test_method = dna.count_nucleotides
		
class DnaToRnaFailed (unittest.TestCase):
	def testT_or_U (self):
		""" The sequence should not contain both Thymine and Uracil """
		self.assertRaises (dna.InvalidSequenceError, self.test_method, NUCLEOTIDES)
		
class FastaTest (unittest.TestCase):
	def testFileExistFailed (self):
		""" File should exist """
		self.assertRaises (tools.InvalidFileLocation, tools.file_to_string('not a place'))

	def testFastaToString (self):
		""" Should return string object without any deletion in it"""
		result = tools.file_to_string(TEST_FASTA_LOC)
		self.assertEqual (TEST_FASTA_STRING, result)

	def testFastaFormatFailed (self):
		""" argument should be in the fasta format """
		self.assertRaises (tools.InvalidFastaFormatError, tools.fasta_to_sequence, 'just a string')
		
	def testFastaRightSequence (self):
		""" should get right sequence from FASTA string """
		tf = """>Rosalind_1497
					GAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATA
					AGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG"""
		res = 'GAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATAAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG'
		self.assertEqual (res, tools.fasta_to_sequence(tf))
		
	def testResultingDNAStringFailed (self):
		""" DNA string should not have any other characters then DNA nucleotides """
		dna = aa.get_ORF (TEST_FASTA_LOC)
		

		
	
		
		
if __name__  == "__main__":
	unittest.main()