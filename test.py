import unittest
import dna
import aa
import tools
import tasks
from collections import Counter


NUCLEOTIDES = "ACTGU"
TEST_FASTA_LOC = 'fasta_test.txt'
TEST_FASTA_STRING = ">Rosalind_1497\nGAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATA\nAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCG\nAAGCGGCACCTAGTCGAGCCTCCCACCCTTATCTTGAACTGAACCAGAGTCCATGGCCGT\nGGTCGTTAGGGGACTATGGAGACCGCTTCATGGGGAGGGAAAACGTTTTGCTATGTACTG\nGGTGCTGTTTTCGGAGTCGTGGGAAGAGACATTTACAAACGTGGCAAAAGCTCTCATGTG\nGCTCCATGTAGGCGTCCAGTTCTCGACATTTCAGTACTCAACGGTCTGCTTCGTTCCGGG\nCAATCTCTTTAGACCCAGCGGGGTCTTATCATCGCGTTTGAAGACATTCCGTGGAAGGGG\nGCGGTCTCGGGAGGATATCGCCTGGGATGGTCTAGCCCGGAGACGTATCTCTGTAGCTAC\nAGAGATACGTCTCCGGGCTAGACCATCGCCCGTTCTAACGATACCCGCGAAGGAAATTAA\nGTACGAGCTTCTGGATATAGATGTTCGTTTCTTCTTATCCCTATGACGTCTCTCATGTAA\nCCGATACTGTTTCCGCTATTGATGTTTACCTGGATAAAGTTTCGTGCCTAGACTTTTAGG\nGGACTCGGCTCCTGTCCCCACAGCCGGTATGAGAAGAGCTCATAAACTCGACAACATCAA\nTGATTCATCTAACTGGGTATTATCAGGTCTGCGTAACCTTTGAGCAATTCGCATGTAGGT\nCTATGGTCGCTCTAGTCCACACGCGCTATGCGCAACACTAGCCCTCGCGGTGTTGCCATT\nGACCGCCTCCAACGCATCCTCATCTGTCCCCGCTCCTGCATTCCGGAACCCAACAAGCTT\nTCAAGATGCATTTTAGCTCAGCTGCTGGGCGATAAACGTCTCTCATTGAGTG\n"
TEST_FASTA_STRING_R = "GAGGGATACAAGTGACAGGCGTAAAGTTGTTATAGCTAGTGACGATCTGTCTTTTTAATAAGCCAGAAGGGGTCTCTTAAAGCGCAAGTATTCAGCGGTTTTTCCATTGGCTGCACCTCGAAGCGGCACCTAGTCGAGCCTCCCACCCTTATCTTGAACTGAACCAGAGTCCATGGCCGTGGTCGTTAGGGGACTATGGAGACCGCTTCATGGGGAGGGAAAACGTTTTGCTATGTACTGGGTGCTGTTTTCGGAGTCGTGGGAAGAGACATTTACAAACGTGGCAAAAGCTCTCATGTGGCTCCATGTAGGCGTCCAGTTCTCGACATTTCAGTACTCAACGGTCTGCTTCGTTCCGGGCAATCTCTTTAGACCCAGCGGGGTCTTATCATCGCGTTTGAAGACATTCCGTGGAAGGGGGCGGTCTCGGGAGGATATCGCCTGGGATGGTCTAGCCCGGAGACGTATCTCTGTAGCTACAGAGATACGTCTCCGGGCTAGACCATCGCCCGTTCTAACGATACCCGCGAAGGAAATTAAGTACGAGCTTCTGGATATAGATGTTCGTTTCTTCTTATCCCTATGACGTCTCTCATGTAACCGATACTGTTTCCGCTATTGATGTTTACCTGGATAAAGTTTCGTGCCTAGACTTTTAGGGGACTCGGCTCCTGTCCCCACAGCCGGTATGAGAAGAGCTCATAAACTCGACAACATCAATGATTCATCTAACTGGGTATTATCAGGTCTGCGTAACCTTTGAGCAATTCGCATGTAGGTCTATGGTCGCTCTAGTCCACACGCGCTATGCGCAACACTAGCCCTCGCGGTGTTGCCATTGACCGCCTCCAACGCATCCTCATCTGTCCCCGCTCCTGCATTCCGGAACCCAACAAGCTTTCAAGATGCATTTTAGCTCAGCTGCTGGGCGATAAACGTCTCTCATTGAGTG"

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

		
class NucleoPatternCase (unittest.TestCase, NucleotideStringTestCase):
	def setUp (self):
		self.test_method = dna.isNucleotide		
		
class CountNucleotidesTest (unittest.TestCase):		
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
		
		second_test = TEST_FASTA_STRING
		second_res = TEST_FASTA_STRING_R
		self.assertEqual (second_res, tools.fasta_to_sequence(second_test))
	# -------------------------
	# FASTA to collection of strings
	def testFastaFormatCollectionFailed (self):
		""" Work only with fasta format string """
		self.assertRaises (tools.InvalidFastaFormatError, tools.fasta_to_collection, 'just a string')
		
	def testFastaRightCollectionSize (self):
		""" The length of the collection should be correct """
		fasta_string = ">Rosalind_1\nGATTACA\n>Rosalind_2\nTAGACCA\n>Rosalind_3\nATACA"
		self.assertEqual (len(tools.fasta_to_collection(fasta_string)), 3)
		
	def testFastaCollectionNames (self):
		""" keys of the dictionary should be correct """
		fasta_string = ">Rosalind_1\nGATTACA\n>Rosalind_2\nTAGACCA\n>Rosalind_3\nATACA"
		names = ["Rosalind_1", "Rosalind_2", "Rosalind_3"]
		r = tools.fasta_to_collection(fasta_string)
		for n in names:
			self.assertIn (n, r)
			
	def testFastaCollectionValues (self):
		""" values in the dictionary should be correct """
		fasta_string = ">Rosalind_1\nGATTACA\n>Rosalind_2\nTAGACCA\n>Rosalind_3\nATACA"
		names = ["Rosalind_1", "Rosalind_2", "Rosalind_3"]
		values = ['GATTACA', 'TAGACCA', 'ATACA']
		r = tools.fasta_to_collection(fasta_string)
		for v in values:
			self.assertIn (v, r.values())
	# --------------------------

class ReturnPeptideTest (unittest.TestCase):
	""" Checked rna to peptide transformation """		
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
		
class TranslationTest (unittest.TestCase):		
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
		self.assertRaises (dna.InvalidSequenceError, aa.translation, 'ACGT')
		
class OpenFrameDetectionTest (unittest.TestCase):
	def testFindAllStartCodons (self):
		""" Should return just all M in sequence"""
		result = [172, 195, 209, 232, 295, 305, 446, 560, 582, 594, 621, 688, 719, 772, 782, 807,
 905]
		data = dna.dnaToRna(TEST_FASTA_STRING_R)
		test = dna.find_all_start_codons (data)
		self.assertEqual (test, result)
		self.assertEqual (len(test), len(result))
	
	def testNoStartCodon (self):
		""" Should return empty array """
		self.assertEqual (dna.get_orf(''), [])
		self.assertEqual (dna.get_orf ('UGAGCAGACAACCUUCAAAAA'), [])
		
	def testNoStopCodon (self):
		""" Should not work with no stop codon if start codon was in the sequence """
		self.assertRaises (aa.NoStopCodonError, dna.get_orf, "AUGAAAGGUAUG")
		
	def testRightAnswer (self):
		""" should return known sequence """
		rna = 'AAA AUG AAA GGU AUG UGA'.replace(" ", "")
		result = [[3,15], [12, 15]]
		self.assertEqual (result, dna.get_orf(rna))

class FindStopCodonTest (unittest.TestCase):
	def testNoStopCodon (self):
		""" Should return false if no stop codon found """
		self.assertEqual(dna.find_next_stop_codon("AAAAAAAAA"), False)

	def testStopCodonFound (self):
		""" Should return right position of stop codon"""
		pos = 3
		self.assertEqual (dna.find_next_stop_codon("AAAUAG"), pos)
		self.assertEqual(dna.find_next_stop_codon("AAGAAAUAGGUU"), 6)

	def testStopCodonNotInRightPlace (self):
		""" When codon found not in the +3 position from start return False """
		self.assertEqual(dna.find_next_stop_codon("AAAAUAGAAA"), False)

class ORF_TestCase (unittest.TestCase):
	def testFileExistFailed (self):
		""" File should exist """
		self.assertRaises (tools.InvalidFileLocationError, tasks.ORF, 'not a place')

	def testKnownResult (self):
		"""Example from Rosalind"""
		known_r = ['MLLGSFRLIPKETLIQVAGSSPCNLS', 'M', 'MGMTPRLGLESLLE', 'MTPRLGLESLLE']
		loc = 'test_files/orf_test.txt'
		
		result = tasks.ORF(loc)
		self.assertEqual (len(result), len(known_r))
		for r in known_r:
			self.assertIn (r, result)	
			
		hard_loc = 'test_files/rosalind_orf_test.txt'	
		hard_dna = tools.fasta_to_sequence(tools.file_to_string(hard_loc))
		hard_r = ['MAYR', 'MISSLNHGLSVSYCRGWHYMYAAPLPHSVHIPIKMGAH', 'MGLPYPSLLGN', 'MPHPYLIVFTSLSRWAHTKAPAPRQFQLPKRLGYGSPIREPVPIHGPSCSHACHRNDR', 'M', 'MV', 'MHV', 'MYAAPLPHSVHIPIKMGAH', 'MSRSWTQICISMPGLGIHP', 'MGA', 'MLMQIWVQERLMILISNPFSTNW', 'MPGLGIHP', 'MLHTLCPYRVQSITNWSKKDCLLRS', 'MIRSCQTPLEPL', 'MRVNVENVKRSRSIRIGSTSWSSGTDLLVRVVLVAFDRNVS', 'MDRGWGPAPEWDSRTPASWVTETDVVPVP', 'MNAKTGHLFAAPCAKLDVMGA', 'MPATAVAHR', 'MHWNAATYSDVHRYRKTNYTFAFYKWSLNGLR', 'MPRPGIYLRLPAQNSTSWGPRIGMILVVAMRVNVENVKRSRSIRIGSTSWSSGTDLLVRVVLVAFDRNVS', 'MTSSFAQGAANRCPVLAFIPS', 'MRAGWTVDGDRLPNGTPVPQPLG', 'MVLEWIALN', 'MCIDTVRPITRLRFINGP', 'MGRTIRPFIATGNRLISNGVCGTL', 'MNAKTGHADADLGPRTTHDLNKQSFFDQLVID', 'MGTRCVT', 'MGTGSRMGLPYPSLLGN','MHASRMDRGWGPAPEWDSRTPASWVTETDVVPVP', 'MHVTETIDEAANLV', 'MATTKIIPILGPHDVEFCAGSRK','MQIWVQERLMILISNPFSTNW', 'MGAH', 'MA', 'MILISNPFSTNW', 'MSQKRSMRQPI', 'MFDWYAWILLRLCPLVIISTR', 'MRQPI', 'MPPPIVMCIDTVRPITRLRFINGP', 'MR', 'MILVVAMRVNVENVKRSRSIRIGSTSWSSGTDLLVRVVLVAFDRNVS', 'MPRPGMLMQIWVQERLMILISNPFSTNW', 'MDCVKLDWLPHRSFL']
			
		hard_result = tasks.ORF(hard_loc)
		self.assertEqual (len(hard_result), len (hard_r))
		for r in hard_r:
			self.assertIn (r, hard_result)

	def testMultipleResult (self):
		""" If result was already made, do not duplicate it"""
		loc = 'test_files/orf_test.txt'
		res = tasks.ORF(loc)
		this_result = Counter(res)
		for r in this_result:
			self.assertEqual(1, this_result[r])

class SharedMotifCase (unittest.TestCase):
	def testKnownResult (self):
		""" Should work with Sample Dataset from Rosalind """
		# at first, and with big dataset afterwords
		pass
		
	def testFileExistFailed (self):
		""" File should exist """
		self.assertRaises (tools.InvalidFileLocationError, tasks.FSH, 'not a place')	
		
if __name__  == "__main__":
	unittest.main()