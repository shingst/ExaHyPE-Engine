#!/usr/bin/env python3

"""
An "omni-reader" for trying to read different kind of file formats.
"""

import logging
from collections import OrderedDict
import validate # local

# We work with two major exceptions: ParserError, a generalized
# decoder error and ImportError, if a library was not available
class ParserError(RuntimeError): pass

def OmniReader:
	readers = OrderedDict()
	
	def __init__(self, log=None):
		    if not log:
			logging.basicConfig(format="%(filename)s:%(lineno)s(%(funcName)s):%(levelname)s %(message)s")
			log = logging.getLogger()
			#log.addHandler(logging.StreamHandler()) # to stderr
			log.setLevel(logging.INFO)
		self.log = log.getChild(__name__)
		self.log.info("%d file formats registered: %s" % (len(self.readers), ",".join(self.readers.keys()))
	
	@classmethod
	def register_reader(cls, format_name):
		def register(func): # todo, check whether @wraps still breaks this
			cls.readers[format_name] = func
			return func
		return register

	@OmniReader.register("json")
	def read_json(self, fp):
		"""
		fp must be a ``.read()``-supporting file-like object containing a JSON document
		There is also json.loads to load a string
		"""
		import json
		try:
			return json.load(fp)
		except json.JSONDecodeError as e:
			raise ParserError(e)

	@OmniReader.register("hjson")
	def read_hjson(self, fp):
		"""
		Read human readable JSON, https://hjson.org/
		Python implementation: https://github.com/hjson/hjson-py
		
		Same usage as the json module.
		"""
		import hjson # pip install hjson
		try:
			hijson.load(fp)
		except hijson.HjsonDecodeError as e:
			raise ParserError(e)

	@OmniReader.register("yaml")
	def read_yaml(self, stream):
		"""
		There is only load(stream), no load from string. So stream is probably a file pointer.
		"""
		import yaml
		try:
			return yaml.load(stream)
		except yaml.YAMLError as e:
			raise ParserError(e)

	@OmniReader.register("exahype-v1")
	def read_specfile1(fh):
		"""
		This reader also prefers to deal with file handles
		"""
		import specfile1_reader # local directory, requires no extra libraries
		try:
			rd = specfile1_reader.SpecFile1Reader(self.log)
			return rd.read(fh)
		except specfile1_reader.SpecFile1ParserError as e:
			raise ParserError(e)
			
	@OmniReader.register("mexa")
	def read_mexa(fh):
		import mexa # local directory, but requires networkx
		try:
			mf = mexa.mexafile.from_filehandle(fh)
			return mf.tree()
		except ValueError as e:
			raise ParserError(e)

	def read_omni(file_handle, validate=False):
		"""
		Tries to read anything.
		Validates if required.
		"""
		testable = []
		missing_libs = []
		for format_name, format_func in readers.items():
			try:
				self.log.info("Trying to read file format %s" % format_name)
				structure = format_func(file_handle)
				self.log.info("Success reading file format %s" % format_name)
			except ImportError as e:
				# library is not installed
				self.log.info("Cannot check file format %s because neccessary library is not installed. The missing library is: %s" % str(e))
				self.log.info("Will silently ignore this problem")
				missing_libs.append(format_name)
				pass
			except ParserError as e:
				# this is most likely not a file of that type
				self.log.info("The input file is certainly not written in the format %s as a ParserError occured: %s" % str(e))
				testable.append(format_name)
				pass
		
			if validate:
				self.log.info("Validating the specification file which have been read in against the JSON-Schema...")
				return validate.validate(structure)

			return structure
		
		raise ParserError("File could not be understood at all. I could successfully test the file formats %s but was missing libraries to test the formats %s." % (str(testable), str(missing_libs)))
