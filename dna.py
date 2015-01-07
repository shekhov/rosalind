""" DNA module. Library to work with DNA sequence """
import re

class DNA_Error (Exception): pass
class InvalidCharacterError (DNA_Error): pass
class InvalidSequenceError (DNA_Error): pass

NUCLEOTIDES = "ATGCU"
N_PATTERN = '[^ATUGC]'

def isNucleotide (sequence): pass
	#if not sequence.isupper(): raise InvalidCharacterError
	#if 'T' in sequence and 'U' in sequence: raise InvalidSequenceError
	#if re.search (N_PATTERN, sequence):	raise InvalidSequenceError



def count_nucleotides (sequence):
	""" Returns an array with quantity of nucleotides inside in the following order (A T G C) """
	result = [0, 0, 0, 0] # A T G C
	isNucleotide(sequence)
	
	for n in sequence:
		if n == 'A': result[0]+=1
		elif n == 'T' or n == 'U': result[1]+=1
		elif n == 'G': result[2]+=1
		elif n == 'C': result[3]+=1
			
	return result
	
def dnaToRna (dna):
	""" Replace Thymines to Uracils """
	isNucleotide(dna)	
	return dna.replace('T', 'U')