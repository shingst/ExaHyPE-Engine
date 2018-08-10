#!/usr/bin/env python3
#
# A token-based parser for the original Specification file format.
# It works a bit like the C++ parser used to work
#

from __future__ import print_function
import re, numbers, logging, inspect, json, pprint, argparse, sys, collections, itertools
from collections import namedtuple

from validate import schema # local namepace

logging.basicConfig(format="%(filename)s:%(lineno)s(%(funcName)s):%(levelname)s %(message)s")
log = logging.getLogger(__name__)
#log.addHandler(logging.StreamHandler()) # to stderr
log.setLevel(logging.INFO)


def grouper(iterable, n, fillvalue=None):
	"Collect data into fixed-length chunks or blocks"
	# grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
	args = [iter(iterable)] * n
	return itertools.zip_longest(*args, fillvalue=fillvalue)

class TokenOrigin(namedtuple("Token", ["name", "line", "col"])): # Parser: token + information where it originates from
	"Mimics the name it holds."
	def __eq__(self, other):
		if isinstance(other, TokenOrigin):
			return self.name == other.name
		if isinstance(other, str):
			return self.name == other
		raise NotImplemented("Comparing TokenOrigin against what?")
	def __str__(self):
		return self.name
	def __hash__(self):
		return hash(self.name)
		
TokenPos = namedtuple("TokenPos", ["token","start","stop"]) # ParserView: token + a slice of the tokenstream

def list_all_methods(object_or_class):
	"List all methods/callables from an object or class"
	return dict(inspect.getmembers(object_or_class, predicate=inspect.isroutine)).keys()

class TokenNotFound(Exception): pass
class MapperFailure(Exception): pass

def do_raise(exception_type, msg):
	raise exception_type(msg)

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
	# const: Currently kills "constants" -> tried with ^$
	seperators = r"((?:\r?\n\r?)+|[\{\}:=,\s+]|\*/|/\*|//|#|^const$)"

	def __init__(self, string=None):
		if string:
			self.tokensWithWitespace = self.tokenize(string)
			self.stringTokens = self.removeComments(self.tokensWithWitespace)

			# without information about lines and columns in original file
			self.tokens = self.makeNative(self.stringTokens)

			# with information about the origin, but as string type since it could be hard
			# to deal with native types
			self.coordinateTokens = self.removeComments(self.addCoordinates(self.tokensWithWitespace))
			
			# you can map then onto each other by the list index
			assert len(self.tokens) == len(self.coordinateTokens)
			
	
	def __str__(self):
		return "Parser(%s)" % (self.tokens)
	def __repr__(self):
		return str(self)
	
	@staticmethod
	def native(token):
		"Parse a string the exahype way"
		if token == "yes" or token == "on":
			return True
		if token == "no" or token == "off":
			return False
		try:
			return int(token)
		except ValueError:
			pass
		try:
			return float(token)
		except ValueError:
			pass
		return str(token)

	@classmethod
	def tokenize(cls, string):
		# First tokenize the whole file, keep the seperators in the stream
		return  [ token for token in re.split(cls.seperators, string) if token != "" ] # don't include empty tokens that carry no information
	
	@classmethod
	def addCoordinates(cls, tokensWithWitespace):
		# Returns tokens with information where it originates from
		line = 0
		col = 0
		origin_token_stream = []
		for token in tokensWithWitespace:
			origin_token_stream.append(TokenOrigin(token, line, col))
			plus_lines = token.count("\n") + token.count("\r")
			if plus_lines > 0:
				line += plus_lines
				col = 0
			else:
				col += len(token)
		return origin_token_stream
	
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
			if inSingleLineComment:
				if re.match(r"^(?:\r?\n\r?)+$", str(token)):
					inSingleLineComment = False
				continue
			if token == "/*":
				inMultiLineComment = True
				continue
			if inMultiLineComment:
				if token == "*/":
					inMultiLineComment = False
				continue
			if token and not re.match(cls.seperators, str(token)):
				# all non-empty strings and non-seperators are real tokens
				stringTokens.append(token)
		return stringTokens
	
	@classmethod
	def makeNative(cls, stringTokens):
		return list(map(cls.native, stringTokens))

class ParserView:
	"""
	This is view on the token stream of the Parser.
	"""

	def __init__(self, tokens, parser=None):
		self.tokens = tokens
		self.parser = parser
		if not parser: log.warn("Some methods need parser access.")
		
	def __str__(self):
		return "ParserView(%s)" % str(self.tokens)
	def __repr__(self):
		return str(self) # whatever
	
	def slice(self, start, stop=None):
		if isinstance(start, collections.Iterable):
			if not isinstance(stop, collections.Iterable):
				stop = [ stop for x in range(len(start)) ] # bring stop to same shape as start
			return [ self.slice(a,b) for a,b in zip(start,stop) ] # threaden this function, return a list of ParserViews
		else:
			return ParserView(self.tokens[slice(start,stop)], self.parser)
	
	def token_after(self, token, additional_tokens_to_skip=0, start=0):
		"""
		Lookup tokens based on a predecessor.
		
		The simplest use is to look up tokens coming after a known token, such as in
		
		> v = ParserView("foo bar baz".split())
		> after = v.token_after("bar")
		TokenPos("baz", 1, None)
		
		However, it can also be used to slice the ParserView:
		
		> v = ParserView("foo bar baz bac".split())
		> token, start, stop = v.token_after("bar")
		> v.slice(start,stop)
		ParserView(["baz", "bac"])
		
		This is especially helpful if used 
		
		If the token was not found, it throws a TokenNotFound exception.
		
		@param start At which index to start to search
		"""
		try:
			pos_of_token = self.tokens.index(token, start)
			pos_of_needle = pos_of_token + 1 + additional_tokens_to_skip
		except ValueError as e: # value not in list
			raise TokenNotFound("Token %s not in tokenlist %s" % (token, self.tokens))
		return TokenPos(self.tokens[pos_of_needle], pos_of_needle, None)
	
	def set_if_token_exists(self, token):
		"""
		Returns the token if it exists, otherwise None.
		The harmless version of token_after().
		Naming comes from optional parameters.
		
		-> TODO: Should not be required, it should be able to just call token_after.
		-> probably deprecated
		"""
		try:
			pos = self.tokens.index(token)
			return TokenPos(token, pos, None)
		except ValueError:
			return TokenPos(None, None, None)
	
	def all_numeric_tokens_after(self, token):
		"""
		Gathers all numeric tokens after a given token in a list and returns the
		slice where these tokens where found.
		"""
		try:
			start = self.tokens.index(token) + 1
		except ValueError as e:
			raise TokenNotFound("Token %s not in tokenlist %s" % (token, self.tokens))
		ret = []
		stop = start
		for token in self.tokens[start:]:
			if isinstance(token, numbers.Number):
				ret.append(token)
				stop += 1
			else:	break
		return TokenPos(ret,start,stop) # ~ TokenPos([1,2,3], 15, 19)
	
	def within_section(self, sectionname, start=0):
		"""
		Returns the slice of a given section. This can be used to subsequently
		search within that section.
		
		@param start An index where to start looking for the section beginning
		"""
		try:
			# Old variants of the 
			startPos = self.tokens.index(sectionname, start)
		except ValueError as e:
			raise TokenNotFound("Missing begin of section '%s'" % sectionname)

		log.debug("Found section %s start at %s" % tuple(map(str, (sectionname, startPos))))
		
		# search for the "end" which belongs to the start. This assumes sections to look like
		#  start sectionname
		#    .... start someting-else ... end something-else ...
		#  end sectionname
		endSectionNamePos = startPos
		while endSectionNamePos < len(self.tokens):
			endSectionName, endSectionNamePos, _ = self.token_after("end", start=endSectionNamePos)
			#endSectionNamePos += startPos # make index again relative to self.tokens
			if endSectionName != sectionname: # continue searching after the current "end"
				log.debug("Searching for end %s, found end %s at %d, looking for next end" % (sectionname, endSectionName, endSectionNamePos))
			else:
				log.debug("Begin %s at %d, found end %s at %d" % (sectionname, startPos, endSectionName, endSectionNamePos))
				return TokenPos(sectionname, startPos, endSectionNamePos)
		raise TokenNotFound("Found start of section %s but cannot find its end." % sectionname)

	def within_deprecated_section(self, sectionname):
		"""
		Try out a sectionname which was superseded. That is, if it does not find anything,
		this will return the whole parserView without modifications.
		
		-> Works, but makes not much sense as a filter if another section is looked up afterwards
		-> probably deprecated
		"""
		try:
			log.debug("Trying out to look for the deprecated sectionname %s" % sectionname)
			res = self.within_section(sectionname)
			return res
		except TokenNotFound:
			log.debug("Did not found the deprecated sectionname %s" % sectionname)
			return TokenPos(None, None, None)
	
	def within_one_of_the_sections(self, *sectionnames):
		for i,sectionname in enumerate(sectionnames):
			try:
				log.debug("Trying out different sections %s, currently %d / %d, sectionname %s" % (str(sectionnames), i, len(sectionnames), sectionname))
				res = self.within_section(sectionname)
				return res
			except TokenNotFound:
				log.debug("Did not found sectionname %s, trying next one from selection" % sectionname)
		# no section found
		raise TokenNotFound("None of the sections %s have been found in the specfile." % str(sectionnames))
	
	def within_sections(self, sectionname):
		"""
		To distinguish lists and to determine how much of them are there, within sections
		allows to count them each by looking for repeated sections
		"""
		sections = []
		currentSectionEnd = 0
		while True:
			try:
				res = self.within_section(sectionname, start=currentSectionEnd)
				currentSectionEnd = res.stop +1 
				sections.append(res)
			except TokenNotFound:
				break # there is no more section (or none at all)
		log.debug("Within sections %s gave %s" % (sectionname, str(sections)))
		#import ipdb; ipdb.set_trace()
		return TokenPos(*zip(*sections)) # return TokenPos([],[],[])

	def detect_variable_list_after(self, introduction_token):
		"""
		Tries to understand the variables instructions which can be something like
		
		  variables = foo:3,bar:2,bla:4
		  something = 3
		
		Since we have =:; as seperators, in this particular example it will be hard
		to understand where the variable list ends. Therefore, we fall back to the
		original parser to learn whether the tokens are in a single list. This is
		against the original grammar but at least a first criterion.
		"""
		token, start, _ = self.token_after(introduction_token)
		
		# Variant 1: "variables = 5"
		if isinstance(token, numbers.Number):
			return TokenPos(token, start, start+1) # no problema, only a number
		
		# We do not support variant "variables = foo,bar,baz" here because it does not
		# belong to the original standard

		variables = {} # name -> multiplicity

		assert self.parser, "Requires access to the parser for heuristics"
		token_row = self.parser.coordinateTokens[start]
		
		# Lookup the coordinates of the tokens from the parser
		# Variant: "variables = foo:3,bar:1,baz:2"
		for (i,name), (j,multiplicity) in grouper(2, enumerate(self.tokens[start:])):
			if isinstance(name, str) and isinstance(multiplicity, int):
				if token_row == self.parser.coordinateTokens[i].line and \
				   token_row == self.parser.coordinateTokens[j].line:
					variables[name] = multiplicity
				else:
					log.warn("In variable list (%s), tokens %s:%d are in different line and therefore we think they no more belong to the variable list" % (introduction_token, name, multiplicity))
					break
			else:
				# the schema no more matches
				break
		stop = j # just for bookkeeping, here
		return TokenPos(variables, start, stop)
	
	def map_token(self, **mapping):
		"""
		Maps the first token in the token string according to mapping.
		"""
		pos = 0
		token = self.tokens[pos]
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
		return TokenPos(*zip(*list_of_tokenpos)) # gives three tuples with a list each

class Required:
	"""
	An object which mimics a boolean with explanation
	"""
	def __init__(self, value, reason=""):
		if isinstance(value, list):
			self.reasons = value # "copy constructor"
		else:
			self.reasons=[]
			self.append(value, reason)
		
	def append(self, value, reason=""):
		self.reasons.append((value,reason))
		return self # chainable!
		
	def isit(self): # is it?
		return all([ i[0] for i in self.reasons ])
	
	def copy(self):
		return Required(list(self.reasons))
	
	def __str__(self):
		return "Required(%s, since [%s])" % (self.isit(),
			", ".join([ "%s -> %s" %(i[1], i[0]) for i in self.reasons])
		)

class Mapper:
	def __init__(self, parser, filterNones=True):
		self.parserView = ParserView(parser)
		self.filterNones = filterNones
	
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
		log.info("Calling %s(%s,%s) on tokenstream=%s" % (funcname, str(args), str(kwargs), str(parserView)))
		return method(*args, **kwargs)
	
	@staticmethod
	def lookup(parserViewInstance, instructions_list, schema):
		"""
		Instruction_list will be worked off,
		Schema only to pass, just in case it is needed.
		"""
		parserView = parserViewInstance
		result = None
		log.debug("Working on instructions_list = " + str(instructions_list))
		for instruction in instructions_list:
			if callable(instruction):
				# This is a hacky solution where the Mapper code himself injects callable
				# closures. Such data cannot come from JSON.
				log.debug("Calling lambda instruction (injected by the Mapper). Will overwrite parserView.")
				parserView = instruction()
				continue
			
			start = None
			stop = None

			if not isinstance(instruction, dict):
				raise ValueError("Instructions must be dictionaries.")
			if isinstance(parserView, list):
				raise TypeError("Lists of ParserViews should not be used like that")


			callable_funcnames = [funcname for funcname in sorted(instruction.keys()) if funcname in list_all_methods(parserView)]
			
			if len(callable_funcnames) > 1:
				raise ValueError("Several callable functions %s in the instructions_list %s. If you want to call several functions, put them as different items in the list, as dicts won't guarantee order." % (str(callable_funcnames), str(instructions_list)))
			
			funcname = callable_funcnames[0]
			ret = Mapper.call(parserView, funcname, instruction[funcname])
			
			if not isinstance(ret, TokenPos):
				raise TypeError("ParserView methods should always return a TokenPos which however might hold arrays.")
			
			result, start, stop = ret
			parserView = parserView.slice(start,stop) # can also be an aray of parserViews at this point
			
			log.debug("%s returned %s and new parserView %s" % (funcname, result, str(parserView)))
			
			#if isinstance(parserView, collections.Iterable):
			#	log.warn("parserView is a list, returning directly")
				# This is because somehow the list property seems to be stripped afte rthe loop and I don't know why.
			#	return (result, parserView)

		log.debug("Lookup returns result,parserView = %s,%s" % (result,parserView))
		return (result, parserView)
		
	@classmethod
	def ensureType(cls, value, schema, path, noneIsOkay=True):
		if not "type" in schema:
			# sometimes the type is externalized
			if "properties" in schema and "$ref" in schema["properties"]:
				log.info("Won't check type of %s since it uses JSON-Pointers and so on. Won't implement.")
				return value
			else:
				raise ValueError("Expect to work on a schema, but missing 'type'. Got %s" % str(schema))
		
		if noneIsOkay and not value: # early catch None types
			return value

		jsontypes = { 'integer': int, 'number': float, 'boolean': bool, 'array': list, "object": dict }
		for json_type, python_type in jsontypes.items():
			if schema["type"] == json_type:
				if not isinstance(value, python_type):
					raise ValueError("At %s, expected %s, but got '%s' from the specfile." % (path, json_type, str(value)))
				
			# This is very basic recursion, covering only the basics of JSON-Schema.
			# If you do something fancy in the schema, this will probably fail.
			if schema["type"] == "array":
				return [ cls.ensureType(elem, schema["items"], cls.join_paths(path,i)) for i,elem in enumerate(value) ]
			if schema["type"] == "object":
				return { k: cls.ensureType(v, schema["items"], cls.join_paths(path,k)) for k,v in value.items() }
		
		return value

	def notNone(self, iterable):
		if self.filterNones:
			return [ i for i in iterable if i] if isinstance(iterable,list) else { k:v for k,v in iterable.items() if v }
		else:	return iterable
		
	@staticmethod
	def join_paths(path, appendix):
		return "%s/%s" % (str(path), str(appendix))
		
	def mapSchema(self, schema, instruction_stack=[], path="/", required=None):
		"""
		Traverse the JSON-Schema file and lookup the positions of the fields
		
		This works recursively.
		"""
		
		if required == None:
			required = Required(True, "Root Schema requires value")
		
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
			
			log.debug("Looking up scalar schema %s" % schema)
			try:
				value, _ = self.lookup(self.parserView, instruction_stack, schema)
				log.debug("Got scalar return value %s" % str(value))
			except TokenNotFound as e:
				if required.isit():
					log.exception(e)
					log.error("Bad specification file, TokenNotFound by the instructions %s" % str(instruction_stack))
					raise MapperFailure("Bad specification file, TokenNotFound by the instructions %s" % str(instruction_stack))
				else:
					# and we even now the reason why this key is not required.
					# Could store it for debugging.
					return None
			return self.ensureType(value, schema, path)

		
		if "type" in schema:
			if schema["type"] == "object":
				if "properties" in schema: # should be there for an object, isn't it?
					log.debug("Mapping dict-schema %s" % schema)
					return self.notNone({ key: self.mapSchema(
						subschema, instruction_stack,
						required=required.copy().append(
							"required" in schema and key in schema["required"],
							"Value in Schema %s" % path),
						path=self.join_paths(path,key)
					) for key, subschema in schema["properties"].items() })
				else:
					log.logInfo("Don't know what to do with object %s" % str(schema))
					return None

			if schema["type"] == "array":
				if "items" in schema: # should be there for an array
					log.debug("Mapping array-schema %s" % schema)
					# We need some knowledge how much items there are awaited.
					# The evaluation of the current instruction stack can tell us
					# that, since it should return an array of mappers
					result, parserViews = self.lookup(self.parserView, instruction_stack, schema)
					if not result:
						return None
					
					if not isinstance(parserViews, collections.Iterable):
						raise MapperFailure("For an array-type schema, expecting the mapping instructions to produce a list of parserViews. However, got (%s,%s)" % (str(result), str(parserViews)))
						
					return self.notNone([ self.mapSchema(
						schema["items"],
						instruction_stack + [lambda: parserView],
						path=self.join_paths(path,i)
					) for i,parserView in enumerate(parserViews)])
		
			if schema["type"] in ["integer", "number","boolean", "string"]:
				return None # no instruction how to map this native type
		# sometimes the type is externalized
		elif "properties" in schema and "$ref" in schema["properties"]:
			log.info("Schema %s is currently not supported as it uses JSON-Pointers and so on.")
			return None

		
		# The requirement detection will be wrong if the requirement is part of an "anyOf" or similar.
		# This is because we don't want to implement a code which really understands JSON-Schema.
		raise ValueError("Cannot understand this schema: %s. It is assumably required (%s)." % (str(schema),str(required)))

	# open tasks:
	# 1) Deal with returning None's in structures
	# 2) consistently replace - by _ in the keywords

class Frontend:
	def start_interactive(self, parser, mapper):
		#import ipdb; ipdb.set_trace()
		
		print("You are given the objects 'parser' and 'mapper'.")
	
		# arguments are exposed to the shell
		try:
			# show fancy IPython console
			from IPython import embed
			embed()
		except ImportError:
			# show standard python console
			import readline, code
			variables = globals().copy()
			variables.update(locals())
			shell = code.InteractiveConsole(variables)
			shell.interact()
	
	def add_arguments(self, parser):
		parser.add_argument("-d", "--debug", action="store_true", help="Show debug messages and, in case of errors, stack traces")
		parser.add_argument("-s", "--shell", action="store_true", help="Start a shell to play with the parser interactively")
		parser.add_argument("-c", "--compact", action="store_false", help="Print output compact in a single line (default is pretty print)")

		parser.add_argument('specfile', nargs="?", type=argparse.FileType('r'), default=sys.stdin, help="The specification file to read in (*.exahype), default stdin")
		parser.add_argument('jsonfile', nargs="?", type=argparse.FileType('w'), default=sys.stdout, help="The output file to write to (*.exahype), default stdout")
		
	def run(self, parser=None):
		if not parser:
			parser = argparse.ArgumentParser(
				description="A translator for the old ExaHyPE specification file format to the new JSON format",
				#epilog="See http://www.exahype.eu and the Guidebook for more help."
			)
			self.add_arguments(parser)

		# also, for interactive colorful output on stderr, check this one:
		# if not sys.stdout.isatty(): # Not an interactive device.

		args = parser.parse_args()
		
		if args.debug:
			log.setLevel(logging.DEBUG)
		
		specfile_as_string = args.specfile.read()
		try:
			parser = Parser(specfile_as_string)
		except Exception as e:
			if args.debug:
				log.exception(e)
				raise e # rerise or so
			else:
				log.fatal(str(e))
		
		mapper = Mapper(parser.tokens, parser)
		
		if args.shell:
			self.start_interactive(parser,mapper)
			sys.exit(0)
			
		json_specfile = mapper.mapSchema(schema)
		
		json_specfile_as_string = json.dumps(
			json_specfile,
			sort_keys = True,
			indent = 4 if args.compact else None
		)
		
		print(json_specfile_as_string, file=args.jsonfile)



if __name__ == '__main__':
	Frontend().run()
