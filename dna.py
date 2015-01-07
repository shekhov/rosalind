""" DNA module. Library to work with DNA sequence """

class DNA_Error (Exception): pass
class InvalidCharacterError (DNA_Error): pass
class InvalidSequenceError (DNA_Error): pass

NUCLEOTIDES = "ATGCU"

def count_nucleotides (sequence):
	""" Returns an array with quantity of nucleotides inside in the following order (A T G C) """
	result = [0, 0, 0, 0] # A T G C
	
	if 'T' and 'U' in sequence: raise InvalidSequenceError	
	for n in sequence:
		if n.upper() not in NUCLEOTIDES: raise InvalidSequenceError
		if n.lower() == n: raise InvalidCharacterError
		if n == 'A': result[0]+=1
		elif n == 'T' or n == 'U': result[1]+=1
		elif n == 'G': result[2]+=1
		elif n == 'C': result[3]+=1
			
	return result
	
def dnaToRna (dna):
	""" Replace Thymines to Uracils """
	if 'T' and 'U' in dna: raise InvalidSequenceError		
	return dna.replace('T', 'U')