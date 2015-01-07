import unittest
import dna
import aa
import tools


NUCLEOTIDES = "ACTGU"
TEST_FASTA_LOC = 'fasta_test.txt'
TEST_FASTA_STRING = ">Rosalind_1497\nGAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATA\nAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG\nAAGCGGCACCTAGTCGAGCCTCCCACCCTTATCTTGAACTGAACCAGAGTCCATGGCCGT\nGGTCGTTAGGGGACTATGGAGACCGCTTCATGGGGAGGGAAAACGTTTTGCTATGTACTG\nGGTGCTGTTTTCGGAGTCGTGGGAAGAGACATTTACAAACGTGGCAAAAGCTCTCATGTG\nGCTCCATGTAGGCGTCCAGTTCTCGACATTTCAGTACTCAACGGTCTGCTTCGTTCCGGG\nCAATCTCTTTAGACCCAGCGGGGTCTTATCATCGCGTTTGAAGACATTCCGTGGAAGGGG\nGCGGTCTCGGGAGGATATCGCCTGGGATGGTCTAGCCCGGAGACGTATCTCTGTAGCTAC\nAGAGATACGTCTCCGGGCTAGACCATCGCCCGTTCTAACGATACCCGCGAAGGAAATTAA\nGTACGAGCTTCTGGATATAGATGTTCGTTTCTTCTTATCCCTATGACGTCTCTCATGTAA\nCCGATACTGTTTCCGCTATTGATGTTTACCTGGATAAAGTTTCGTGCCTAGACTTTTAGG\nGGACTCGGCTCCTGTCCCCACAGCCGGTATGAGAAGAGCTCATAAACTCGACAACATCAA\nTGATTCATCTAACTGGGTATTATCAGGTCTGCGTAACCTTTGAGCAATTCGCATGTAGGT\nCTATGGTCGCTCTAGTCCACACGCGCTATGCGCAACACTAGCCCTCGCGGTGTTGCCATT\nGACCGCCTCCAACGCATCCTCATCTGTCCCCGCTCCTGCATTCCGGAACCCAACAAGCTT\nTCAAGATGCATTTTAGCTCAGCTGCTGGGCGATAAACGTCTCTCATTGAGTG\n"

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
		self.assertRaises (dna.InvalidSequenceError, dna.dnaToRna, NUCLEOTIDES)
		
class FastaTest (unittest.TestCase):
	# File To string Test
	def testFileExistFailed (self):
		""" File should exist """
		self.assertRaises (tools.InvalidFileLocationError, tools.file_to_string, 'not a place')

	def testFastaToString (self):
		""" Should return string object without any deletion in it"""
		result = tools.file_to_string(TEST_FASTA_LOC)
		self.assertEqual (TEST_FASTA_STRING, result)
	# --------------------------
	# FASTA to sequence Test
	def testFastaFormatFailed (self):
		""" argument should be in the fasta format """
		self.assertRaises (tools.InvalidFastaFormatError, tools.fasta_to_sequence, 'just a string')
		
	def testFastaRightSequence (self):
		""" should get right sequence from FASTA string """
		tf = ">Rosalind_1497\nGAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATA\nAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG"
		res = 'GAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATAAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG'
		self.assertEqual (res, tools.fasta_to_sequence(tf))
	# -------------------------
		

		
	
		
		
if __name__  == "__main__":
	unittest.main()