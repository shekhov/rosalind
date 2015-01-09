""" DNA module. Library to work with DNA sequence """
import re
import copy
import aa
from aa import RNA_code
from aa import START_CODON
import tools

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
	
def getSetOfMotifs (sequence, minL=None, maxL=None):
	""" Return the list of all posible motifs from given sequence """
	if maxL is None: maxL = len(sequence)-1
	if minL is None: minL = 1
	
	thisMax = 0	
	motifs = set()
	
	for n in range (len(sequence)):		
		if maxL > n: thisMax = n
		else: thisMax = maxL
		
		for l in range (minL, thisMax):
			m = sequence [n-l:n]
			motifs.add(m)

	return motifs
	
def findSimilarMotif (array, percCut=None):
	""" Return dictionary with only motifs that are present in all strings """
	# percCut is responsible for what percent of minimum motif length we should cut off 
	if percCut is None: percCut = 0
		
	result = {}
	# Initializative filling 
	start_set = getSetOfMotifs (array[0])
	for each in start_set:
		result[each] = True
	
	maxMotifL = 0		
	for s in array[1:]:
		this_set = getSetOfMotifs(s, minL = int(maxMotifL*percCut*0.01))
		old_res = copy.copy(result)

		for motif in old_res:
			if motif not in this_set:
				del result[motif]
				
		maxMotifL = len(tools.getLongestKey(result))		
		print ("%i -> %i maxL: %i" % (len(old_res), len(result), maxMotifL))	
				
	return result
	