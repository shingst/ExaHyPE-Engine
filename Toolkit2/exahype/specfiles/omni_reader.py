#!/usr/bin/env python3

"""
An "omni-reader" for trying to read different kind of file formats.
"""

from collections import OrderedDict
import validate # local

# We work with two major exceptions: ParserError, a generalized
# decoder error and ImportError, if a library was not available
class ParserError(RuntimeError): pass

readers = OrderedDict()
def register(func, format_name):
	def register(the_func):
		readers[format_name] = the_func
		return the_func
	return register

@register("json")
def read_json(string):
	import json
	try:
		return json.loads(string)
	except json.JSONDecodeError as e:
		raise ParserError(e)

@register("hjson")
def read_hjson(string):
	"""
	Read human readable JSON, https://hjson.org/
	Python implementation: https://github.com/hjson/hjson-py
	"""
	import hjson # pip install hjson
	try:
		hijson.loads(string)
	except hijson.HjsonDecodeError as e:
		raise ParserError(e)

@register("yaml")
def read_yaml(stream):
	import yaml
	try:
		return yaml.load(stream)
	except yaml.YAMLError as e:
		raise ParserError(e)
	
@register("mexa")
def read_mexa(file_handle):
	import mexa # local directory, but requires networkx
	try:
		mf = mexa.mexafile.from_filehandle(file_handle)
		return mf.tree()
	except ValueError as e:
		raise ParserError(e)

def read_omni(string, validate=False):
	"""
	Tries to read anything.
	Validates if required
	"""
	for format_name, format_func in readers.items():
		try:
			structure = format_func(string)
		except ImportError as e:
			# library is not installed
			pass
		except ParserError as e:
			# this is most likely not a file of that type
			pass
	
		if validate:
			return validate.validate(structure)

		return structure
