#
# A token-based parser for the original Specification file format.
# It works a bit like the C++ parser used to work
#

import re, numbers
from collections import namedtuple

Token = namedtuple("Token", ["name", "pos"])

class Parser:
	def __init__(self, string):
		# First tokenize the whole file, keep the seperators in the stream
		# The regex includes: (line endings|general seperators|multiline comments|single line comments|the "const" keyword)
		seps = r"((?:\r?\n\r?)+|[\{\}:=,\s+]|\*/|/\*|//|#|const)"
		#seps = r"[\{\}:=,\s+]"
		self.tokensWithWitespace = re.split(seps, string)
		
		self.stringTokens = []
		inSingleLineComment = False
		inMultiLineComment = False
		for token in self.tokensWithWitespace:
			if token == "//" or token == "#":
				inSingleLineComment = True
				continue
			if inSingleLineComment and re.match(r"^(?:\r?\n\r?)+$", token):
				inSingleLineComment = False
				continue
			if token == "/*":
				inMultiLineComment = True
				continue
			if inMultiLineComment and token == "*/":
				inMultiLineComment = False
				continue
			if token and not re.match(seps, token):
				# all non-empty strings and non-seperators are real tokens
				self.stringTokens.append(token)
				
		# convert tokenstream to native python types
		for token in self.stringTokens:
			if token == "yes":
				self.tokens.append(True)
				continue
			if token == "no":
				self.tokens.append(False)
				continue
			
			try:
				self.tokens.append(int(token))
				continue
			except ValueError:
				pass

			try:
				self.tokens.append(float(token))
				continue
			except ValueError:
				pass
			
			# probably just a string
			self.tokens.append(token)

class Mapper:
	def __init__(self, tokens):
		self.tokens = tokens
	
	def slice(self, start, stop=None):
		return Mapper(self.tokens[start:] if stop is None else self.tokens[start:stop])
	
	def token_after(self, token, additional_tokens_to_skip=0):
		pos = self.tokens.index(token) # raises ValueError, could return tokenNotFound
		tokenpos = pos + additional_tokens_to_skip
		return Token(self.tokens[tokenpos], tokenpos)
	
	def all_numeric_tokens_after(self, token):
		start = self.tokens.index(token) # raises ValueError, could return tokenNotFound
		ret = []
		stop = pos
		for token in self.tokens[pos:]:
			if instanceof(token, numbers.Number):
				ret.append(token)
				stop += 1
			else:	break
		return Token(ret, range(start, stop))
	
	def within_section(self, sectionname):
		beginSectionName, start = self.token_after(sectionname)
		
		endSectionName = None
		stop = None
		while endSectionName != beginSectionName:
			endSectionName, stop = self.slice(start,stop).token_after("end")
		
		return Token(beginSectionName, range(start, stop))
	
	def set_if_token_exists(self, token):
		try:
			pos = self.tokens.index(token)
			return Token(token, pos)
		except ValueError:
			return Token(None, None)
	
	def map_token(self, **mapping):
		"""
		Maps the first token in the token string according to mapping.
		"""
		pos = 0
		token = self.token[pos]
		print("Mapping token %s according to mapping %s" % (token, str(mapping)))
		if token in mapping:
			return Token(mapping[token], pos)
		elif:
			raise ValueError("Cannot map %s because it is not in the mapping %s" % (token, str(mapping)))
		
	def filter(self, *tokens):
		"""
		Returns a list of tokens which are in the tokenstream.
		"""
		ret = [ token for token in tokens if token in self.tokens ]
		pos = [ self.tokens.index(token) for token in ret ]
		return Token(ret, pos)
	
	def call(self, funcname, something):
		# makeup an argument list
		if isinstance(something, list):
			args = something
			kwargs = {}
		elif isinstance(something, dict):
			args = []
			kwargs = something
		else:
			args = [something]
			kwargs = {}
		method = getattr(self, funcname)
		return method(*args, **kwargs)
	
	def lookup(self, instructions_list):
		mapper = self
		for instruction in instructions_list:
			start = None
			stop = None
			if not isinstance(instruction, dict):
				raise ValueError("Instructions must be dictionaries.")
			call = lambda mapper, funcname: mapper.call(funcname, instruction[funcname])
			if "token_after" in instruction:
				result, start = call(mapper, "token_after")
			elif "all_numeric_tokens_after" in instruction:
				result, (start,stop) = call(mapper, "all_numeric_tokens_after")
			elif "within_section" in instruction:
				result, (start,stop) = call(mapper, "within_section")
			elif "set_if_token_exists" in instruction:
				result, start = call(mapper, "set_if_token_exists")
			elif "map_token" in instruction:
				result, start = call(mapper, "map_token")
			elif "filter_token" in instruction:
				result, positions = call(mapper, "filter_token")
			mapper = self.slice(start,stop)
		return result

	def mapSchema(self, schema, instruction_stack=[]):
		"""
		Traverse the JSON-Schema file and lookup the positions of the fields
		
		This works recursively.
		"""
		if "old_format_prepend" in schema:
			instruction_stack = [ schema["old_format_prepend"] ] + instruction_stack
		
		if "old_format" in schema:
			self.lookup(instruction_stack + [schema["old_format"]])
		
		if "properties" in schema:
			for key, subschema in schema["properties"].items():
				return self.mapSchema(subschema, instruction_stack)
		


	
	# replace - by _ in the keywords

test = """
exahype-project  Euler

  peano-kernel-path const          = ./Peano
  exahype-path const               = ./ExaHyPE
  output-directory const           = ./ApplicationExamples/EulerFlow/EulerFlow
  architecture const               = noarch

  computational-domain
    dimension const          = 3/*
    width                    = 15.0,*/ 15.0, 15.0
    offset                   = 0.0, //0.0, 0.0
    end-time                 =// 0.1
  end computational-domain

  shared-memory
    identifier             #  = dummy
    configure                = {}
    cores                    = 2
    properties-file          = sharedmemory.properties
  end shared-memory

  distributed-memory
    identifier               = static_load_balancing
    configure                = {hotspot,fair,ranks-per-node:4}
    buffer-size              = 64
"""

p = Parser(test)
print("Full tokens:")
print(p.fullTokenStream)
print("Simple tokens:")
print(p.tokens)
