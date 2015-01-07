""" Script with helpful tools """
import os.path

class ToolsError (Exception): pass
class InvalidFileLocationError (ToolsError): pass
class InvalidFastaFormatError (ToolsError): pass

def file_to_string (location):
	""" return strings from the file """
	if not os.path.isfile (location): raise InvalidFileLocationError
	file = open (location, 'r')
	result = ""
	for line in file.readlines():
		result = result + line
	file.close()
	return result

def fasta_to_sequence (fasta_string):
	""" :return: string of dna without fasta formating """
	#string should have > characters with following name. Check only > character
	if '>' not in fasta_string: raise InvalidFastaFormatError
	# Cut names of fasta format and merge the rest of the sequence in one string
	
	result_s = ""
	n = fasta_string.find('>')
	while (n < len (fasta_string)):
		end_name = fasta_string[n:].find("\n")
		next_name = fasta_string[end_name:].find ('>')
		if next_name == -1: n = len(fasta_string)
		else: n = next_name
		result_s += fasta_string[end_name:n].replace("\n", "")
		
	return result_s
		
