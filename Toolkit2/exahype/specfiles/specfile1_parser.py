#
# A token-based parser for the original Specification file format.
# It works a bit like the C++ parser used to work
#

from __future__ import print_function
import re, numbers, logging, inspect, json
from collections import namedtuple

from validate import schema # local namepace

logging.basicConfig(format="%(filename)s:%(lineno)s(%(funcName)s):%(levelname)s %(message)s")
log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler()) # to stderr
log.setLevel(logging.DEBUG)


Token = namedtuple("Token", ["name", "pos"])

TokenPos = namedtuple("TokenPos", ["token","start","stop"])

def list_all_methods(object_or_class):
	"List all methods/callables from an object or class"
	return dict(inspect.getmembers(object_or_class, predicate=inspect.isroutine)).keys()

class TokenNotFound(Exception):
	pass

class Parser:
	"""
	This is a simple tokenizer for ExaHyPE Specification Files (v1, the "original" ones).
	It understands single line comments and multi line comments and it acts pretty close
	to the C++ parser for the original specfile format.
	
	Usage is like:
	
	> p = Parser()
	> print(p.tokens) # is a list of native tokens
	
	"""
	
	# (line_endings|general_seperators|comments|const_keyword)
	# general_seperators: probably ":" should not be part (vtk::foo::bar cf. algorithmic-factor:2), or : should be but not ::
	seperators = r"((?:\r?\n\r?)+|[\{\}:=,\s+]|\*/|/\*|//|#|const)"

	def __init__(self, string=None):
		if string:
			self.tokensWithWitespace = self.tokenize(string)
			self.stringTokens = self.removeComments(self.tokensWithWitespace)
			self.tokens = self.nativeTokens(self.stringTokens)

	@classmethod
	def tokenize(cls, string):
		# First tokenize the whole file, keep the seperators in the stream
		return  [ token for token in re.split(cls.seperators, string) if token != "" ] # don't include empty tokens that carry no information
	
	@classmethod
	def removeComments(cls, tokensWithWitespace):
		# Returns tokens without comments and whitespace, but as a list of strings
		stringTokens = []
		inSingleLineComment = False
		inMultiLineComment = False
		for token in tokensWithWitespace:
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
			if token and not re.match(cls.seperators, token):
				# all non-empty strings and non-seperators are real tokens
				stringTokens.append(token)
		return stringTokens
	
	@staticmethod
	def nativeTokens(stringTokens):
		tokens = []
		# convert tokenstream to native python types
		for token in stringTokens:
			if token == "yes":
				tokens.append(True)
				continue
			if token == "no":
				tokens.append(False)
				continue
			
			try:
				tokens.append(int(token))
				continue
			except ValueError:
				pass

			try:
				tokens.append(float(token))
				continue
			except ValueError:
				pass
			
			# probably just a string
			tokens.append(token)
		return tokens

class ParserView:
	"""
	This is view on the token stream of the Parser.
	"""

	def __init__(self, tokens=None):
		if isinstance(tokens, list):
			self.tokens = tokens
		elif isinstance(tokens, Parser):
			self.tokens = tokens.tokens
		else:
			raise ValueError("Expect list of tokens or Parser instance")
	
	def slice(self, start, stop=None):
		if isinstance(start, list):
			if not isinstance(stop, list):
				stop = [ stop for x in range(len(start)) ] # bring stop to same shape as start
			return [ self.slice(a,b) for a,b in zip(start,stop)] # threaden this function, return a list of ParserViews
		else:
			return ParserView(self.tokens[slice(start,stop)])
	
	def token_after(self, token, additional_tokens_to_skip=0):
		try:
			pos = self.tokens.index(token)
		except ValueError as e:
			# value not in list
			raise TokenNotFound("Token %s not in tokenlist %s" % (token, self.tokens))
		tokenpos = pos + additional_tokens_to_skip
		return TokenPos(self.tokens[tokenpos], tokenpos, None)
	
	def all_numeric_tokens_after(self, token):
		try:
			start = self.tokens.index(token)
		except ValueError as e:
			raise TokenNotFound("Token %s not in tokenlist %s" % (token, self.tokens))
		ret = []
		stop = pos
		for token in self.tokens[pos:]:
			if instanceof(token, numbers.Number):
				ret.append(token)
				stop += 1
			else:	break
		return TokenPos(ret, start, stop)
	
	def within_section(self, sectionname):
		beginSectionName, start, stop = self.token_after(sectionname)
		
		# search for the "end" which belongs to the start. This assumes sections to look like
		#  start sectionname
		#    .... start someting-else ... end something-else ...
		#  end sectionname
		endSectionName = None
		while endSectionName != beginSectionName:
			endSectionName, stop = self.slice(start,stop).token_after("end")
		
		return TokenPos(beginSectionName, start, stop)
	
	def set_if_token_exists(self, token):
		"""
		Returns the token if it exists, otherwise None.
		The harmless version of token_after().
		"""		
		try:
			pos = self.tokens.index(token)
			return TokenPos(token, pos, None)
		except ValueError:
			return TokenPos(None, None, None)
	
	def map_token(self, **mapping):
		"""
		Maps the first token in the token string according to mapping.
		"""
		pos = 0
		token = self.token[pos]
		print("Mapping token %s according to mapping %s" % (token, str(mapping)))
		if token in mapping:
			return TokenPos(mapping[token], pos, None)
		else:
			raise ValueError("Cannot map %s because it is not in the mapping %s" % (token, str(mapping)))
		
	def filter(self, *tokens):
		"""
		Returns a list of tokens which are in the tokenstream.
		"""
		list_of_tokenpos = map(self.set_if_token_exists, tokens)
		return zip(*list_of_tokenpos) # gives three tuples with a list each


class Mapper:
	def __init__(self, parser):
		self.parserView = ParserView(parser)
	
	@staticmethod
	def call(parserView, funcname, something):
		"""
		makeup an argument list from something, where something typically comes from JSON
		and is of any tye (scalar, list, dict).
		"""
		if isinstance(something, list):
			args = something
			kwargs = {}
		elif isinstance(something, dict):
			args = []
			kwargs = something
		else:
			args = [something]
			kwargs = {}
		method = getattr(parserView, funcname)
		log.info("Calling %s(%s,%s) on tokenstream=%s" % (funcname, str(args), str(kwargs), str(parserView.tokens)))
		return method(*args, **kwargs)
	
	@staticmethod
	def lookup(parserViewInstance, instructions_list, schema):
		"""
		Instruction_list will be worked off,
		Schema only to pass, just in case it is needed.
		"""
		parserView = parserViewInstance
		for instruction in instructions_list:
			if callable(instruction):
				# This is a hacky solution where the Mapper code himself injects callable
				# closures. Such data cannot come from JSON.
				result, parserview = instruction()
			
			start = None
			stop = None
			if not isinstance(instruction, dict):
				raise ValueError("Instructions must be dictionaries.")

			for funcname in sorted(instruction.keys()):
				# Note: If you want several instructions to be called, especially in order,
				# put them as an Instruction_list, not inside a single dict.
				if funcname in list_all_methods(parserView):
					result, start, stop = Mapper.call(parserView, funcname, instruction[funcname])
					parserView = Mapper.slice(start,stop)

		return (result, parserView)
	
	@staticmethod
	def ensureType(value, schema):
		if not "type" in schema:
			raise ValueError("Expect to work on a schema, but missing type. Got %s" % str(schema))

		if schema["type"] == "int":
			try:	return int(value)
			except ValueError:
				raise ValueError("Expected integer, but have %s" % str(value))
		if schema["type"] == "number":
			try:	return float(value)
			except ValueError:
				raise ValueError("Expected number, but have %s" % str(value))
		if schema["type"] == "boolean":
			try:	return bool(value)
			except ValueError:
				raise ValueError("Expected bool, but have %s" % str(value))
		if schema["type"] == "array":
			if not isinstance(value, list):
				raise ValueError("Expected list, but have %s" % str(value))
			# This is very basic recursion, covering only the basics of JSON-Schema.
			# If you do something fancy in the schema, this will probably fail.
			return map(lambda elem: ensureType(elem, schema["items"]), value)
		if schema["type"] == "object":
			raise ValueError("Not implemented")
		
	def mapSchema(self, schema, instruction_stack=[], required=True):
		"""
		Traverse the JSON-Schema file and lookup the positions of the fields
		
		This works recursively.
		"""
		
		if "anyOf" in schema or "allOf" in schema:
			pass
		
		if "old_format_prepend" in schema:
			instruction = schema["old_format_prepend"]
			if isinstance(instruction, list):
				instruction_stack = instruction_stack + instruction
			elif isinstance(instruction, dict):
				instruction_stack = instruction_stack + [ instruction ]
			else:
				raise ValueError("Don't know what to do with instruction %s" % str(instruction))
			
		if "old_format" in schema:
			# this is directly the order to map this thing, without going deeper in the hierarchy.
			instruction = schema["old_format"]
			if isinstance(instruction, list):
				instruction_stack = instruction_stack + instruction
			elif isinstance(instruction, dict):
				instruction_stack = instruction_stack + [ instruction ]
			else:
				raise ValueError("Don't know what to do with instruction %s" % str(instruction))
			
			try:
				value, _ = self.lookup(self.parserView, instruction_stack, schema)
			except TokenNotFound as e:
				if required:
					raise ValueError("Bad specfile") from e
				else:
					return None
			return self.ensureType(value, schema)

		
		if "type" in schema:
			if schema["type"] == "object":
				isrequired = lambda key: key in schema["required"] if "required" in schema else False
				if "properties" in schema: # should be there for an object, isn't it?
					log.debug("Schema has properties")
					return { key: self.mapSchema(subschema, instruction_stack, isrequired(key)) for key, subschema in schema["properties"].items() }
				else:
					log.logInfo("Don't know what to do with object %s" % str(schema))
					return None

			if schema["type"] == "array":
				if "items" in schema: # should be there for an array
					# We need some knowledge how much items there are awaited.
					# The evaluation of the current instruction stack can tell us
					# that, since it should return an array of mappers
					_, parserViews = self.lookup(self.parserView, instruction_stack, schema)
					
					def mapArrayItem(parserView):
						arrayItemFilter = lambda parserView: TokenPos(None, parserView)
						return self.mapSchema(schema["items"], instruction_stack + [arrayItemFilter])
						
					return map(mapArrayItem, parserViews)
		
			if schema["type"] in ["int", "number","boolean", "string"]:
				return None # no instruction how to map this native type
		
		raise ValueError("Cannot understand this schema: " + str(schema))

	# open tasks:
	# 1) Deal with returning None's in structures
	# 2) consistently replace - by _ in the keywords

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
  end distributed-memory
"""

p = Parser(test)
#print("Full tokens:")
#print(p.tokensWithWitespace)
print("Simple tokens:")
print(p.tokens)

m = Mapper(p)
m.mapSchema(schema)
