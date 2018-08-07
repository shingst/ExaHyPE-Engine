
class BadSpecificationFile(Exception):
	pass
	
class TemplateNotFound(Exception):
	_templateFile = None
	def __init__(self,templateFile):
		self._templateFile=templateFile
	pass
