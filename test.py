import unittest
import dna
import aa
import tools


NUCLEOTIDES = "ACTGU"
TEST_FASTA_LOC = 'fasta_test.txt'
TEST_FASTA_STRING = ">Rosalind_1497\nGAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATA\nAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG\nAAGCGGCACCTAGTCGAGCCTCCCACCCTTATCTTGAACTGAACCAGAGTCCATGGCCGT\nGGTCGTTAGGGGACTATGGAGACCGCTTCATGGGGAGGGAAAACGTTTTGCTATGTACTG\nGGTGCTGTTTTCGGAGTCGTGGGAAGAGACATTTACAAACGTGGCAAAAGCTCTCATGTG\nGCTCCATGTAGGCGTCCAGTTCTCGACATTTCAGTACTCAACGGTCTGCTTCGTTCCGGG\nCAATCTCTTTAGACCCAGCGGGGTCTTATCATCGCGTTTGAAGACATTCCGTGGAAGGGG\nGCGGTCTCGGGAGGATATCGCCTGGGATGGTCTAGCCCGGAGACGTATCTCTGTAGCTAC\nAGAGATACGTCTCCGGGCTAGACCATCGCCCGTTCTAACGATACCCGCGAAGGAAATTAA\nGTACGAGCTTCTGGATATAGATGTTCGTTTCTTCTTATCCCTATGACGTCTCTCATGTAA\nCCGATACTGTTTCCGCTATTGATGTTTACCTGGATAAAGTTTCGTGCCTAGACTTTTAGG\nGGACTCGGCTCCTGTCCCCACAGCCGGTATGAGAAGAGCTCATAAACTCGACAACATCAA\nTGATTCATCTAACTGGGTATTATCAGGTCTGCGTAACCTTTGAGCAATTCGCATGTAGGT\nCTATGGTCGCTCTAGTCCACACGCGCTATGCGCAACACTAGCCCTCGCGGTGTTGCCATT\nGACCGCCTCCAACGCATCCTCATCTGTCCCCGCTCCTGCATTCCGGAACCCAACAAGCTT\nTCAAGATGCATTTTAGCTCAGCTGCTGGGCGATAAACGTCTCTCATTGAGTG\n"

class NucleotideStringTestCase:
	# All methods working with sequences of DNA and RNA should be tested here
	#test_method = None
	def testUpperCase (self):
		""" All letters should be in uppercase """
		self.assertRaises (dna.InvalidCharacterError, self.test_method, "acgc")
			
	def testNucleotides (self):
		""" should work only with 5 letters: A T G C U """
		self.assertRaises (dna.InvalidSequenceError, self.test_method, 'RFK')

	def testT_or_U (self):
		""" The sequence should not contain both Thymine and Uracil """
		self.assertRaises (dna.InvalidSequenceError, self.test_method, NUCLEOTIDES)
		self.assertRaises (dna.InvalidSequenceError, self.test_method, "TTTCCCUUU")

class CountNucleotidesTest (unittest.TestCase, NucleotideStringTestCase):
	def setUp (self):	
		self.test_method = dna.count_nucleotides
		
	def testArray(self):
		""" count_nucleotides should return an array with 4 numbers """
		self.assertEqual (4, len(dna.count_nucleotides("ATGC")))
	
	def testResult(self):
		""" count_nucleotides should return right amount of nucleotides """
		sequence = "AAATTTTGGGGCCCCC"
		result = [3, 4, 4, 5]
		self.assertEqual (result, dna.count_nucleotides(sequence))
		
class DnaToRnaFailed (unittest.TestCase):
	def testT_or_U (self):
		""" The sequence should not contain both Thymine and Uracil """
		self.assertRaises (dna.InvalidSequenceError, dna.dnaToRna, NUCLEOTIDES)
		
class NucleoPatternCase (unittest.TestCase, NucleotideStringTestCase):
	def setUp (self):
		self.test_method = dna.isNucleotide
		
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

class ReturnPeptideTest (unittest.TestCase, NucleotideStringTestCase):
	""" Checked rna to peptide transformation """
	def setUp (self):
		self.test_method = aa.return_peptide
		
	def testGivenRNA(self):
		"""Should not contain T"""
		self.assertRaises (dna.InvalidSequenceError, aa.return_peptide, 'ACGT')
	
	def testReturningValue (self):
		""" Should return known sequences """
		p = ["MAKR", "MDNL"]
		rna = ["AUGGCAAAAAGA", "AUGGACAACCUU"]
		result = []
		for n in range (len(rna)):
			self.assertEqual (p[n], aa.return_peptide(rna[n]))
		
	def testStopCodonDetect (self):
		""" Should stop when meet stop codon """
		p = "M"
		rna = "AUGUGAGCA"
		self.assertEqual (p, aa.return_peptide(rna))
		
	def testNotByThreeDevided (self):
		""" Should work fine with sequences /3 != 0 """
		self.assertEqual ("M", aa.return_peptide("AUGAA"))
		
class TranslationTest (unittest.TestCase, NucleotideStringTestCase):
	def setUp (self):
		self.test_method = aa.translation
		
	def testStartCodonDetect (self):
		""" Should skip all codons before start codon """
		p = "MA"
		rna = "GCC UUC AUG GCA UGA".replace(" ", "")
		result = aa.translation(rna)[0]
		self.assertEqual (p, result)

	def testEmptyResult (self):
		"""If empty argument, or no start codon was present return empty array"""
		self.assertEqual(aa.translation(""), [])
		self.assertEqual(aa.translation("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), [])

	def testNoStopCodonGiven (self):
		"""Should raise error if no codon was present, when start codon was """
		self.assertRaises (aa.NoStopCodonError, aa.translation, "AUGAAAGGUAUG")

	def testGivenRNA(self):
		"""Should not contain T"""
		self.assertRaises (dna.InvalidSequenceError, aa.transcription, 'ACGT')
		
class OpenFrameDetectionTest (unittest.TestCase):

	def testNoStartCodon (self):
		""" Should return empty array """
		self.assertEqual (aa.get_orf(''), [])
		self.assertEqual (aa.get_orf ('UGAGCAGACAACCUUCAAAAA'), [])
		
	def testNoStopCodon (self):
		""" Should not work with no stop codon if start codon was in the sequence """
		self.assertRaises (aa.NoStopCodonError, aa.get_orf, "AUGAAAGGUAUG")
		
	def testRightAnswer (self):
		""" should return known sequence """
		rna = 'AAA AUG AAA GGU AUG UGA'.replace(" ", "")
		result = [[3,15], [12, 15]]
		self.assertEqual (result, aa.get_orf(rna))

class FindStopCodonTest (unittest.TestCase):

	def testNoStopCodon (self):
		""" Should return false if no stop codon found """
		self.assertEqual(aa.find_next_stop_codon("AAAAAAAAA"), False)

	def testStopCodonFound (self):
		""" Should return right position of stop codon"""
		pos = 3
		self.assertEqual (aa.find_next_stop_codon("AAAUAG"), pos)
		self.assertEqual(aa.find_next_stop_codon("AAGAAAUAGGUU"), 6)

	def testStopCodonNotInRightPlace (self):
		""" When codon found not in the +3 position from start return False """
		self.assertEqual(aa.find_next_stop_codon("AAAAUAGAAA"), False)
		
	
		
		
if __name__  == "__main__":
	unittest.main()