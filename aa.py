""" Script for working with amino acids and protein sequences """

import tools
import dna

class AminoAcidError (Exception): pass
class NoStopCodonError (AminoAcidError): pass
class InvalidCodonError (AminoAcidError): pass


START_CODON = "AUG"

STOP_CODONS = ['UAG', 'UGA', 'UAA']

RNA_code = {
    'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
    'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
    'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
    'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
    'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
    'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
    'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'TAA': False,   'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'TAG': False,   'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
    'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'TGA': False,   'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}
AA_M = {'D': 115, 'E': 129, 'F': 147, 'G': 57, 'A': 71, 'C': 103, 'L': 113, 'M': 131, 'N': 114, 'H': 137, 'I': 113, 'K': 128, 'T': 101, 'V': 99, 'W': 186, 'P': 97, 'Q': 128, 'R': 156, 'S': 87, 'Y': 163}

AA_EXACT_MASS = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}

AA_REV_M = {128: 'KQ', 97: 'P', 99: 'V', 131: 'M', 101: 'T', 103: 'C', 129: 'E', 137: 'H', 71: 'A', 113: 'IL', 114: 'N', 115: 'D', 147: 'F', 87: 'S', 57: 'G', 163: 'Y', 156: 'R', 186: 'W'}

MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

def get_aa (tRNA):
	""" Return aminoacid from RNA_code """
	pass

def return_peptide (RNA):
	""" Return correspondance peptide based on given RNA """
	# Should work with already made sequence without stop codons. Only one peptide per call

	dna.isNucleotide(RNA)
	if 'T' in RNA: raise dna.InvalidSequenceError	
	peptide = ""
	id = 0
	while id < len(RNA):
		tRNA = RNA[id:id+3]
		if not tRNA in RNA_code: break
		AA = RNA_code[tRNA]
		if AA == False: # In the case if was mistake in the previous step
			break
		peptide += AA
		id += 3
	
	return peptide

def find_next_stop_codon (RNA):
	""" Codon (+3) dependant """
	dna.isNucleotide(RNA)
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
	start_codon = RNA.find(START_CODON)
	if start_codon == -1: return result

	frame = [start_codon]
	id = start_codon
	while id < len(RNA):
		stop = find_next_stop_codon(RNA[id:])
		if not stop:
			if len(result) == 0: raise NoStopCodonError #Exit of no stop-codon in the whole sequence
			else: break #Exit with already one result
		frame.append(stop + id)
		result.append(frame)

		#Next iteration
		start_codon = RNA[id+1:].find(START_CODON)
		if start_codon == -1: break #Exit if no more codons around

		frame = [start_codon+id+1]
		id = frame[0]
	return result



def translation (RNA):
	""" transform rna to corresponding peptides, start-stop codon dependant """
	result = []
	if len(RNA) == 0: return result
	dna.isNucleotide(RNA)
	if "T" in RNA: raise dna.InvalidSequenceError

	orf = get_orf(RNA)

	for frame in orf:
		peptide = return_peptide(RNA[frame[0]:frame[1]+3])
		result.append(peptide)
	return result

