import unittest
import dna

NUCLEOTIDES = "ACTGU"

class CountNucleotidesTest (unittest.TestCase):

	def testArray(self):
		""" count_nucleotides should return an array with 4 numbers """
		self.assertEqual (4, len(dna.count_nucleotides("")))
	
	def testResult(self):
		""" count_nucleotides should return right amount of nucleotides """
		sequence = "AAATTTTGGGGCCCCC"
		result = [3, 4, 4, 5]
		self.assertEqual (result, dna.count_nucleotides(sequence))

class NucleotideStringTestCase (unittest.TestCase):
	# All methods working with sequences of DNA and RNA should be tested here
	def __init__ (self, method_to_test):
		super (NucleotideStringTestCase, self).__init__()
		self.test_method = method_to_test	

	def testUpperCase (self):
		""" All letters should be in uppercase """
		self.assertRaises (dna.InvalidCharacterError, method_to_test, "acgc")
			
	def testNucleotides (self):
		""" count_nucleotides should work only with 5 letters: A T G C U """
		self.assertRaises (dna.InvalidSequenceError, method_to_test, 'R')

	def testT_or_U (self):
		""" The sequence should not contain both Thymine and Uracil """
		self.assertRaises (dna.InvalidSequenceError, method_to_test, NUCLEOTIDES)	
	
		
		
class DNAStringCheckFailed (NucleotideStringTestCase):
	def __init__(self, method):
		super(DNAStringCheckFailed, self).__init__(dna.count_nucleotides)

		
	
		
		
if __name__  == "__main__":
	unittest.main()