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
	if maxL is None: maxL = len(sequence)
	if minL is None: minL = 1
		
	#print ("Search from %i to %i lenght" % (minL, maxL))
	
	thisMax = 0	
	motifs = set()
	
	for n in range (minL-1, len (sequence)):
		if maxL >= n: thisMax = n+1
		else: thisMax = maxL
		
		for l in range (minL-1, maxL): # Length 
			m = sequence [n-l: n+1]
			if (len(m))==0:continue
			#print ("n = %i, max = %i, l=%i, m = %s" % (n, thisMax, l, m))
			motifs.add(m)
	
	# for n in range (minL-1, len(sequence)):		
		# if maxL > n: thisMax = n+1
		# else: thisMax = maxL
		
		# for l in range (minL, thisMax):
			# m = sequence [n-l:n]
			# motifs.add(m)

	return motifs
	
def findSimilarMotif (array):
	""" Return dictionary with only motifs that are present in all strings """
	init_motif = getSetOfMotifs (array[0])
	
	dic = {}
	for motif in init_motif:
		#result[motif] = True
		l = len(motif)
		if l in dic: dic[l][motif] = True
		else: dic[l] = {}; dic[l][motif] = True
		
	
	for s in array[1:]: # For all given strings
		c = copy.deepcopy (dic) # Make copy for allowing righting in original dictionary
		found = False
		while (True): # Goes through dictionary until we meet element that presents in next string
			this_max = max(dic) # Take the maximum index from deleted dictionary
			next_motif = getSetOfMotifs (s, this_max, this_max) # Search for this index in the string
			if len (next_motif) == 0: # If nothing in the next string with this length, go again with next lenght
				del (dic[this_max])
				continue
			
			for el in c[this_max]: # In the array search for equal element in motifs
				if el not in next_motif: # If in the next motif there is no such element
					del (dic[this_max][el]) # delete this element from dictionary for the next simplifying
					if (len (dic[this_max]) == 0): del(dic[this_max])
				else: 
					found = True# IF MATCH, then job for this iteration is done, and we can go for another sequence
			if found: break	
			if len (dic) == 0: break
	return dic

def getLongestMotifs (dic):
	""" arg: dic is a dictionary with key=length of elements in the values """
	maxKey = max(dic)
	keys = list (dic[maxKey].keys())
	return keys	
	