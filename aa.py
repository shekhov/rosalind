""" Script for working with amino acids and protein sequences """

import tools

RNA_code = {'UCC': 'S', 'AUG': 'M', 'AUA': 'I', 'CAA': 'Q', 'AUC': 'I', 'GUG': 'V', 'GAG': 'E', 'UAG': False, 'GUC': 'V', 'UCA': 'S', 'GUA': 'V', 'AUU': 'I', 'UGC': 'C', 'UCU': 'S', 'UGU': 'C', 'UAU': 'Y', 'UCG': 'S', 'GUU': 'V', 'GCU': 'A', 'UUC': 'F', 'ACA': 'T', 'AGC': 'S', 'GAA': 'E', 'AGG': 'R', 'GCG': 'A', 'GCA': 'A', 'GCC': 'A', 'GGA': 'G', 'GGC': 'G', 'ACC': 'T', 'GGG': 'G', 'UUA': 'L', 'CAU': 'H', 'CCU': 'P', 'GGU': 'G', 'UUG': 'L', 'AAA': 'K', 'UAA': False, 'CGG': 'R', 'CGA': 'R', 'CGC': 'R', 'ACU': 'T', 'CAG': 'Q', 'ACG': 'T', 'CCC': 'P', 'CAC': 'H', 'UAC': 'Y', 'CCG': 'P', 'CGU': 'R', 'AAC': 'N', 'AAU': 'N', 'CCA': 'P', 'UGA': False, 'CUU': 'L', 'AGU': 'S', 'CUC': 'L', 'GAC': 'D', 'CUA': 'L', 'CUG': 'L', 'GAU': 'D', 'UGG': 'W', 'AAG': 'K', 'AGA': 'R', 'UUU': 'F'}

AA_M = {'D': 115, 'E': 129, 'F': 147, 'G': 57, 'A': 71, 'C': 103, 'L': 113, 'M': 131, 'N': 114, 'H': 137, 'I': 113, 'K': 128, 'T': 101, 'V': 99, 'W': 186, 'P': 97, 'Q': 128, 'R': 156, 'S': 87, 'Y': 163}

AA_EXACT_MASS = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}

AA_REV_M = {128: 'KQ', 97: 'P', 99: 'V', 131: 'M', 101: 'T', 103: 'C', 129: 'E', 137: 'H', 71: 'A', 113: 'IL', 114: 'N', 115: 'D', 147: 'F', 87: 'S', 57: 'G', 163: 'Y', 156: 'R', 186: 'W'}

MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

def get_aa (tRNA):
	""" Return aminoacid from RNA_code """
	pass

def return_peptide (RNA):
	""" Return correspondance peptide based on given RNA """
	pass
