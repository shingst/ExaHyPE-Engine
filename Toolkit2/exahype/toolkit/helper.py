import fnmatch
import os

def find(parent_directory,file_glob):
	"""
	Searches `parent_directory` recursively for files
	matching `file_glob`.

	Returns a list of paths to files 
	matching `file_glob`

	NOTE: Function works for Python 2.2 and above.
	Starting with Python 3.5, more elegant solutions are available.
	"""
	matches = []
	for root, dirnames, filenames in os.walk(parent_directory):
		for filename in fnmatch.filter(filenames, file_glob):
			matches.append(os.path.join(parent_directory, filename))
	return matches

class BadSpecificationFile(Exception):
	pass
	
class TemplateNotFound(Exception):
	_templateFile = None
	def __init__(self,templateFile):
		self._templateFile=templateFile
	pass
	
def parse_variables(solver,field):
	"""
	Parse the 'variables' parameter and similarly
	structured parameters ('material_parameters','global_observables')
	"""
	if field in solver:
		if type(solver[field]) is int:
			return [{"name" : "Q", "multiplicity" : solver[field]}]
		else:
			return solver[field]
	else:
		return [{"name" : "Q", "multiplicity" : 0}]
		
def variables_to_str(solver,field):
	"""
	Represent a variables field as string.
	"""
	if field in solver:
		if type(solver[field]) is int:
			return str(solver[field])
		else:
			return ", ".join([item["name"]+" : "+str(item["multiplicity"]) for item in solver[field]])
	else:
		return "0"

def count_variables(variables_as_list):
	"""
	Sum up the multiplicities.
	"""
	number = 0;
	for variable in variables_as_list:
		number += variable["multiplicity"]
	return number
