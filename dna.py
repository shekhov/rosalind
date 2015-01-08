""" DNA module. Library to work with DNA sequence """
import re
import copy
import aa
from aa import RNA_code
from aa import START_CODON

class DNA_Error (Exception): pass
class InvalidCharacterError (DNA_Error): pass
class InvalidSequenceError (DNA_Error): pass

NUCLEOTIDES = "ATGCU"
N_PATTERN = '[^ATUGC]'

def isNucleotide (sequence): #pass
	if not sequence.isupper(): raise InvalidCharacterError
	if 'T' in sequence and 'U' in sequence: raise InvalidSequenceError
	if re.search (N_PATTERN, sequence):	raise InvalidSequenceError

def reverse_dna (sequence):
	""" Return complimentared sequence and turn it the way it will be read by peptides """
	# Can be done by one line
	return(sequence[::-1].translate(str.maketrans('ACGT', 'TGCA')))
	
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
	
def find_all_start_codons (RNA):
	""" Return an array with indexes of AGU sequence"""
	result = []
	for n in range (len(RNA)):
		if RNA[n:n+3] == START_CODON:
			result.append (n)
	return (result)
	
def find_next_stop_codon (RNA):
	""" Codon (+3) dependant """
	id = 0
	while id < len(RNA):
		tRNA = RNA[id:id+3]
		if not tRNA in RNA_code: return False
		if RNA_code[tRNA] == False:
			return id
		id+= 3
	return False

def get_orf (RNA):
	""" Return array with indexes of open frames """
	result = []
	
	start_codons = find_all_start_codons (RNA)
	for codon in start_codons:
		next_stop = find_next_stop_codon (RNA[codon:])
		if not next_stop:
			if len(result) == 0: raise aa.NoStopCodonError #Exit of no stop-codon in the whole sequence
			else: break #Exit if no more codons around
			
		result.append ([codon, codon + next_stop])
	return result
	
def getSetOfMotifs (sequence):
	""" Return the list of all posible motifs from given sequence """
	motifs = set()
	for n in range (len(sequence)): 
		l = n+1 # lenght of the motif
		for t in range (len(sequence)-1):
			m = sequence[t: t+l]
			motifs.add(m)
	return motifs
	
def findSimilarMotif (array):
	""" Return dictionary with only motifs that are present in all strings """
	result = {}
	# Initializative filling 
	start_set = getSetOfMotifs (array[0])
	for each in start_set:
		result[each] = True
		
	for s in array[1:]:
		this_set = getSetOfMotifs(s)
		for motif in copy.copy(result):
			if motif not in this_set:
				del result[motif]
				
	return result
	