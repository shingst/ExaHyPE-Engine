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
