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

def fasta_to_collection (fasta_string):
	""" Return dictionary with keys=identificators of fasta string, value=correspondent sequence"""
	if '>' not in fasta_string: raise InvalidFastaFormatError
	result = {}
	# Find all enterence of >name\n
	name_start = []
	name_end = []
	for i in range (len(fasta_string)):
		if fasta_string[i] == '>': 
			name_start.append (i+1)
			this_end = fasta_string[i:].find('\n')
			if this_end == -1: raise InvalidFastaFormatError
			name_end.append (i+this_end)
	name_start.append (len(fasta_string)+1) # Need for not throwing the error in the next step		
	# Build dictionary
	for n in range (len(name_end)):
		name = fasta_string[name_start[n]:name_end[n]]
		value = fasta_string[name_end[n]+1:name_start[n+1]-1]
		result[name] = value.replace("\n", "")
	return result
	
def convert_array_to_string (array, sep):
	""" Convert given array to the string, where elements of array separates by sep """
	return sep.join(str(element) for element in array)
