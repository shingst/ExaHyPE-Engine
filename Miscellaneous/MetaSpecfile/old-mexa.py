#!/usr/bin/python
# Regular python 2

"""
This is Mexa, the ExaHyPE meta specification file format. Actually this is
mexa.py, the reference implementation for reading and converting this file
format.

The Mexa file format is a simple and generic hierarchic parameter file format
inspired by the Cactus parameter format.

This python implementation allows to read in a Mexa file and dump it's content
in various other configuration formats.

Written by SvenK in Nov 2017.
"""

# batteries:
import os, re, sys, ast, inspect, argparse, base64, pprint, itertools, collections
from argparse import Namespace as namespace
from collections import namedtuple
from itertools import izip
#from future.utils import raise_from # no batteries

###
### Helper functions and shorthands
###

baseflags = re.IGNORECASE
match = lambda pattern,string,flags=0: re.match(pattern,string,baseflags+flags)
unpack = lambda f: lambda p: f(*p)
# The identity function. Don't mix up with pythons internal "id"
idfunc = lambda x: x
# A nicer functional list access thing
first = lambda x: x[0]
last = lambda x: x[-1]
# Removes all None from a list
removeNone = lambda l: [ x for x in l if x is not None ]
removeFalse = lambda l: [ x for x in l if x ]
# removes the list wrapper if it contains only one item
unlistOne = lambda l: l[0] if len(l)==1 else l
# Gives the list of duplicates from a list
dupes = lambda a: [x for n, x in enumerate(a) if x in a[:n]]
# Insert list b at position i in list a, replacing element a[i]
replace_with_list_at = lambda i, a, b: a[:i]+b+a[i+1:]
# Remove common prefix of two strings, i.e. remove_prefix("abcdef","abc")=="def"
remove_comstr_prefix = lambda text, prefix: text[text.startswith(prefix) and len(prefix):]
# flatten a 2d list
flatten2d = lambda l: [item for sublist in l for item in sublist]
# unique items
unique = lambda l: list(set(l))
# unique items while preserve the order
def unique_preserve(seq):
	seen = set()
	seen_add = seen.add
	return [x for x in seq if not (x in seen or seen_add(x))]
# invert dictionary
invdict = lambda d: {v: k for k, v in d.iteritems()}
# an own isninstance function which also allows mapping if types is a list.
# Reads "is a" as the Perl isa function. Mimics isinstance.
def isa(obj, types):
	if isinstance(types, list):
		return any([isa(obj, t) for t in types])
	return isinstance(obj, types)
# Gives an iterator over a dict or list
anykey = lambda obj: obj.keys() if isa(obj,dict) else xrange(len(obj))
# walks a complex dict or list structure
def mapComplex(func, node):
	for key in anykey(node):
		item = node[key]
		if isa(item,dict) or isa(item,list):
			mapComplex(func, item)
		else:
			node[key] = func(item)
	return node
# remove elements from list by indexlist
withoutIndices = lambda lst, indices: [i for j, i in enumerate(lst) if j not in indices]
# raise an exception. Useful from lambdas.
def raise_exception(e):
	raise e
# returns string until first occurance of character, not including the character
untilCharacter = lambda txt, char: first(txt.partition(char))
#
def quoted_printable(s, escape='%'):
	"""
	Returns a quoted printable version of the string s with escape character %, no maximum line length
	(i.e. everything in a single line) and everything non-alphanumeric replaced.
	"""
	# sourcecode inspired by quopri, https://github.com/python/cpython/blob/2.7/Lib/quopri.py
	HEX = '0123456789ABCDEF'
	quote = lambda c: escape + HEX[ord(c)//16] + HEX[ord(c)%16] # quote a single character
	needsquote = lambda c: not ('0' <= c <= '9' or 'a' <= c <= 'z' or 'A' <= c <= 'Z') # Whether character needs to be quoted
	return "".join([ quote(c) if needsquote(c) else c for c in s])

###
### Definition of the symbol classes
###

class symbol:
	"""
	A symbol is an hierarchical identifier, similar to symbols in LISP and Mathematica.
	Symbols are the atoms of the mexa language.
	Symbols are treated as immutable: There is no method to change them after construction.
	In Mexa, symbols always appear as the LHS of an relation (operation).
	"""
	
	# primitive types which we allow to store on the RHS.
	primitives = [int,float,str,unicode,bool]

	def __init__(self,name=''):
		""""
		Creates a symbol from a name such as Foo/Bar/Baz or a path such
		as ['Foo','Bar','Baz']. Without argument, this gives the root symbol
		symbol().
		Since a symbol is immutable, we can store a string and list version at the
		same time for lookup efficiency.
		"""
		if isa(name,(list,tuple)):
			self._path = map(str.lower, map(str.strip, name))
		elif isa(name, symbol):
			self._path = name._path # kind of copy constructor, but only references to other list.
		else:
			name = str(name).lower() # in order to get case-insensitive
			self._path = map(str.strip, re.split(lang.symb_split_or, name))
		# in case of roots and merged paths and so on
		self._path = removeFalse(self._path)
		# make the path immutable
		self._path = tuple(self._path)
		# Compute a string representation, useful for path computations
		self._canonical = "/".join(self._path)
		
	def canonical(self):
		"Canonical string representation with a specific seperator"
		return self._canonical
	def __repr__(self):
		return "symbol(%s)" % self.canonical()
	def isRoot(self):
		return len(self._path) == 0
	
	def node(self):
		if self.isRoot():
			return self
		return symbol(self._path[-1])
	def parent(self):
		if self.isRoot():
			raise ValueError("Symbol is already root")
		return symbol(self._path[:-1])
	def ancestors(self):
		"""
		Lists all parenting symbols, for instance for /a/b/c it is the
		list [/, /a, /a/b].
		"""
		return [symbol(self._path[:i]) for i,_ in enumerate(self._path)]
	
	def add_prefix(self, other, inplace=False): # should be named add_prefix
		"Returns a new symbol wich is prefixed by self."
		if not isa(other,symbol): other = symbol(other)
		newpath = other._path + self._path;
		if inplace:
			self._path = newpath
			return self
		else:	return symbol(newpath)
	
	def remove_prefix(self, other, inplace=False):
		"Returns a new symbol which has removed the common prefix"
		if not isa(other,symbol): other = symbol(other)
		newpath = remove_comstr_prefix(self.canonical(), other.canonical())
		if inplace:
			self._path = newpath
			return self
		else:	return symbol(newpath)

	# allow symbols to be dict keys:
	def __hash__(self):
		return hash(tuple(self._path))
	def __eq__(self, other):
		return self.canonical() == other.canonical()
	def __lt__(self, other): # sortable
		return self.canonical() < other.canonical()

# deprecated:
def hasnosymbol(rhs):
	# todo: assert that lists are no more supported.
	if type(rhs) in symbol.primitives:
		return True
	if type(rhs) == list:
		return all(map(iscomplete,rhs))
	if type(rhs) == dict:
		raise ValueError("Dictionaries are not supported as right hand sides in Operation %s" % str(self))

class VariableContext:
	"Context for evaluating variables."
	def __init__(self, op, operations):
		self.op = op
		self.operations = operations
	
	def resolve_symbol(self, symbol):
		return self.operations.resolve_symbol(symbol, relative_to=self.op.lhs, src=self.op.rhs)

def rhs2python(text, context=None):
	"""
	Parses a right hand side into int, float, boolean, symbol, whatever.
	This does return functions in some cases. They are supposed to be executed with
	an operations instance then for variable substitution.
	If this returns a list or if the returned and subsequently evaluated function
	returns a list, it is subject to the caller to deal with this list accordingly.
	
	context shall be a VariableContext instance or None.
	"""
	text = text.strip()
	
	# as a kind of preprocessor, support "@foo"-like variables. In case they
	# are detected (we uglily try to detect them not in comments), we return
	# instead a promise function for later evaluation with a correct variable
	# context.
	text_decommented = first(text.partition("#"))
	if re.search(lang.stringvarsimple, text_decommented) or re.search(lang.stringvarcomplex, text_decommented):
		#import ipdb; ipdb.set_trace()
		def promise_rhs2python(context, text=text): # text=text is a workaround for py2 missing closures
			replmatch = lambda matchobj: context.resolve_symbol(matchobj.group(1))
			text = re.sub(lang.stringvarsimple, replmatch, text, count=0, flags=baseflags)
			text = re.sub(lang.stringvarcomplex, replmatch, text, count=0, flags=baseflags)
			return rhs2python(text, context=context)
		promise_rhs2python.__repr__ = "evalstr(%s)" % text
		# if we already have an operations context, directly evaluate the variables.
		# Otherweise, return the closure.
		return promise_rhs2python(context) if context else promise_rhs2python
	
	# We use python's parser for both understanding the data structure
	# and removing pythonic line comments such as "#".
	try:
		typed = ast.literal_eval(text)
		if isa(typed,symbol.primitives):
			return typed
		elif isa(typed,[tuple,list]):
			# We said we support lists.
			return list(typed)
		else:
			# we probably have a created a weird type, for instance we
			# parsed a dictionary, which we do not want to support.
			# Instead, thus, return a text. This may not correctly strip comments,
			# thought.
			return text
	except (ValueError, SyntaxError):
		# This is most likely a string or something weird
		
		if re.match(r"^(Yes|True|On)"+lang.linecomment, text, re.IGNORECASE):
			return True
		if re.match(r"^(No|False|Off)"+lang.linecomment, text, re.IGNORECASE):
			return False
		
		# A symbol or a list of symbols.
		symbs = re.match(r"^("+lang.symb+")(\s+"+lang.symb+")*"+lang.linecomment, text, re.IGNORECASE)
		# If only one symbol is detected, return the symbol. Otherwise
		# return the list of symbols which has to be understood correctly.
		if symbs:
		      return unlistOne(map(symbol, removeNone(symbs.groups())))
		
		return text
	return text

class operation:
	"""
	An operation is a relationship which relates lhs and rhs together. We make use of
	Python types here extensively, i.e. operationships are differed by their type.
	Operations shall always have symbols as their lhs and should have sourced instances
	as their right hand sides, i.e. an operation is of type tuple<symbol,sourced>.
	"""
	
	def __init__(self, lhs, rhs):
		self.lhs, self.rhs = lhs, rhs
	def __repr__(self):
		return "%s(%s, %s)" % (self.__class__.__name__, self.lhs, self.rhs)
	def iscomplete(self):
		return hasnosymbol(self.rhs)
	
	def as_tuple(op, symbolmapper): # was op2tuple
		""""
		Extract an op to a plain python object, stripping the objects and relation in it
		@args symbolmapper a Function which maps a symbol to something else
		"""
		def rhsValue2py(irhs):
			if   isa(irhs, sourced): # unwrapping
				return rhsValue2py(irhs.value)
			elif isa(irhs, list): # threading
				raise ValueError("Lists no more supported.") # this is a hint for a bad parser
				return map(rhsValue2py, irhs)
			elif isa(irhs, symbol):
				return symbolmapper(irhs)
			else:
				return irhs
			
		fst = op.lhs.canonical()
		if isa(op,assign) and isa(op.rhs,list):
			raise ValueError("Assignments no more supported in get_tuple") # also a hint for a bad parser
			snd = map(rhsValue2py, op.rhs)
		else:
			snd = rhsValue2py(op.rhs.value)
		return (fst,snd)

	def as_str(op, symbolmapper=None): # was op2str
		"""
		Converts an operation to a string, looking similar to the parsed input of the operator
		In contrast to op2tuple (which strips the operator), this one does not handle lists on the
		rhs. It expects them to be removed at some calling step, if neccessary. Otherwise they just
		get strings.
		"""
		if not symbolmapper: symbolmapper = symbol.canonical
		symbop = invdict(operation.symbols)
		s_lhs = op.lhs.canonical()
		s_op = symbop[type(op)]
		if isinstance(op.rhs.value, symbol):
			s_rhs = symbolmapper(op.rhs.value)
		else:
			delim = '"' if type(op.rhs.value) in [str,unicode] else ''
			s_rhs = delim+str(op.rhs.value)+delim
		return "%s %s %s" % (s_lhs, s_op, s_rhs)
	
	def add_prefix(self, path):
		"Return a new operation where lhs is prefixed with path"
		return self.__class__( self.lhs.add_prefix(path), self.rhs)
	
	def remove_prefix(self, path):
		"Return a new operation list where each lhs has a common prefix with path removed"
		return self.__class__( self.lhs.remove_prefix(path), self.rhs)
	
	def strip_source(self):
		"Return a new operation where the source is dropped, if present"
		value = self.rhs.value if isa(self.rhs, sourced) else self.rhs
		return self.__class__( self.lhs, value)
	
	@classmethod
	def from_textline(cls, srcline): # was line2operation
		"""
		Parses a single line. Input is a sourceline object. You can test lines quickly with
		> line2operation(sourceline.from_unknown("foo=10"))
		"""
		if re.match(r"^"+lang.comment+"|^\s*$", srcline.text): # comments
			return None

		parts = re.match(r"^(?P<lhs>(?:"+lang.symb+"|))\s*(?P<op>(?:"+lang.opor+"))\s*(?P<rhs>.+)\s*$", srcline.text, re.IGNORECASE)
		if not parts:
			raise ValueError("Don't understand line %d in %s: '%s'"%(srcline.linenum,srcline.fname,srcline.text))

		# split lhs into name path
		name = symbol(parts.group('lhs'))
		# parse rhs as symbol or whatever
		value = rhs2python(parts.group('rhs'))

		return operation.symbols[parts.group('op')](name, sourced(value, srcline))
	
	def evaluate(self, operations):
		"""
		This method may take the operations as read only input parameter and otherwise
		return a list of operations which replace this very command.
		"""
		raise ValueError("Evaluate not implemented for %s" % str(self))
	
	#@classmethod
	#def register(cls, operations):
	#	"""
	#	Called once in an operations instance: Register this type of relationships
	#	within the operations.
	#	"""
	#	raise ValueError("Registration not implemented for %s" % str(cls))
	
	def evaluate_rhs(self, context):
		# since rhs2python is already called in from_textline
		if callable(self.rhs.value):
			self.rhs.value = self.rhs.value(context)
		return self.rhs.value
	#	new_value = rhs2python(self.rhs.value, context=context)
		#self.rhs.value = new_value
	#	return new_value

	
class let(operation):
	"overwritable: can be overwritten. We want to get rid of this type"
	pass

class equals(operation):
	"""
	Equalities. They have interesting properties:
	* They are immutable: Set only once.
	* They are real identities which go in both directions,
	  i.e. a=b will affect not only a but also b if both
	  are nodes (trees), not leafs.
	"""
	def evaluate(op, operations):
		# init:
		if not hasattr(operations, "equals_persistent"):
			# prepare persistent data structures on operations
			operations.equals_persistent = namespace()
			operations.equals_persistent.node_classes = equivalence_classes()
			#unnneeded so far# equals_node_isListed = dict() # maps symbol -> bool
		node_classes = operations.equals_persistent.node_classes # abbreviation
		
		# work:
		if isa(op.rhs.value, list):
			# lists get created by rhs2python i.e. with something like a = (1,2,3)
			def listToAppends(i,vi): #for i,vi in enumerate(op.rhs.value):
				if isa(vi,symbol.primitives):
					return equals(
						symbol("l%d"%i).add_prefix(op.lhs),
						sourced(vi, "list assignment rule"))
				else:
					raise ValueError("Found illegal non-primitive value in RHS of assignment operation "+str(op))
			return map(listToAppends, enumerate(op.rhs.value))
		elif isa(op.rhs.value, symbol):
			# sort equality into operation
			node_classes.add_operation(op)
			# flag operation for being analyzed in a later sweep.
			# The undirected equivalence is encoded in the resolution of it.
			return [
				# the ordinary interpretation "a=b means b sets a"
				#extends(op.lhs, node_classes.get_sourced_resolver(op.rhs)),
				extends_symmetric(op.lhs, op.rhs)
				# the reverse and mathematical interpretation, "a=b means also a sets b"
				###extends(op.rhs.value, op.rhs.derive(value=node_classes.get_resolver(op.lhs), src="LHS")),
			]
		elif isa(op.rhs.value,symbol.primitives):
			# no more replacement.
			return [define(op.lhs,op.rhs)]
		else:
			raise ValueError("Illegal non-primitive assignment at RHS of "+str(op))

class subsets(operation):
	""""
	An equality in only one direction:
	a < b means a gets all properties from b but not the other way around such as a = b.
	Similarly to equals(), this property is transitive: a < b < c means also a < c. Therefore,
	we also come up with a directed equivalence class (hierarchy graph).
	"""
	
	def evaluate(op, operations):
		# init:
		if not hasattr(operations, "subsets_persistent"):
			# prepare persistent data structures on operations
			operations.subsets_persistent = namespace()
			operations.subsets_persistent.node_graph = hierarchy_graph()
		node_graph = operations.subsets_persistent.node_graph # abbreviation
		
		
		# ensure that we extend only from symbols
		if not isa(op.rhs.value, symbol):
			raise ValueError("We only can extend from other symbols. In %s" % str(op))

		node_graph.add_operation(op)
		# flag operation for being analyzed in a later sweep
		return [
			extends(op.lhs, node_graph.get_sourced_resolver(op.rhs))
		]

		### DO THIS IN A SECOND STEP: 
		### include another tree structure
		##queried_symbol = op.rhs.value
		##ext_oplist = operations.query(queried_symbol)
		##return op.lhs.prefix_oplist(queried_symbol.remove_prefix_oplist(ext_oplist))

class extends(operation):
	"""
	extends(a,b) means that a is extended by b.
	This relationship is generated by equals(a,b) => [ extends(a,b), extends(b,a) ]
	and                              subsets(a,b) => [ extends(a,b) ]
	
	When evaluating the extend operation, we make use of the delayed equivalence resolver, i.e.
	this is a seconds-step operation.
	"""
	def evaluate(op, operations):
		ret = []
		
		if not isa(op.rhs.value, equivalence_resolver.delayed_resolver):
			raise ValueError("extend excepts a delayed equivalence resolver. Got instead: %s" % str(op.rhs))
		
		for symb in op.rhs.value.get():
			# basically inherit from all these symbols
			#import ipdb; ipdb.set_trace()
			ret += operations.query(symb, exclude_root=True).remove_prefix(symb).add_prefix(op.lhs).ops
		return ret

class extends_symmetric(operation):
	def evaluate(op, operations):
		ret = []
		
		if not isa(op.rhs.value, symbol):
			raise ValueError("extend_symmetric expects a symbol: %s" % str(op.rhs))
		
		# list<sourced<symbol>>
		equals = operations.equals_persistent.node_classes.get(op.lhs)
		
		has_leaf = any([ operations.is_leaf(symb) for symb in equals ])
		all_empty = all([ operations.is_empty_node(symb) for symb in equals ])

		if has_leaf and all_empty:
			for symb in equals:
				if operations.is_leaf(symb):
					ret += [define(op.lhs, operations.resolve_leaf(symb))]
		elif (has_leaf and not all_empty):
			raise ValueError("In the equivalence class of %s, there is a leaf but also nodes: %s" % (op.lhs, str(equals)))
		elif (has_leaf and all_empty):
			raise ValueError("There is no symbol for %s. All I have is: %s" % (op.lhs, str(equals)))
		else:
			for symb in equals:
				# inherit like "a=b means b sets a"
				ret += operations.query(symb, exclude_root=True).remove_prefix(symb).add_prefix(op.lhs).ops
				# inherit the other way around -- quickly for the time being here.
				# the reverse and mathematical interpretation, "a=b means also a sets b"
				##### ret += operations.query(symb, exclude_root=True).remove_prefix(symb).add_prefix(op.rhs.value).ops
		
		if not len(ret):
			raise ValueError("Replacement %s brought no result. Equivalence class contains %s. Maybe you forgot to define %s?" % (op.lhs, equals, op.rhs))
		
		return ret


class define(operation):
	"""
	An actual definition. One-of-a-kind. The end poduct of equality evaluation.
	"""
	pass


class append(operation):
	"""
	Appending is just an abbreviation for assigning.
	Therefore this operation has a simple 1:1 replacement rule.
	"""
	def evaluate(op, operations):
		# This algorithm reads as rule: "a/b += c/d" => "a/b/d = c/d"
		
		if isa(op.rhs.value,list):
			# we do support multiple RHS values, i.e. a syntax like
			#  a += b c d  <=> equals(a, sourced([b,c,d], ...))
			return flatten2d([o.evaluate(operations) for o in op.rhs.value])
		
		# name of the node to create below the lhs. This is by
		# definition the node name.
		if not isa(op.rhs.value,symbol):
			raise ValueError("Only support symbol for appending, got %s in %s" % (type(op.rhs.value), op))
		name = op.rhs.value.node()
		return [ equals(op.lhs.add_prefix(name), op.rhs) ]
	
class include(operation):
	"Represents an inclusion"
	
	def evaluate(op, operations):
		# include a file
		fname = op.evaluate_rhs(VariableContext(op,operations))

		if not isa(fname, [str,unicode]):
			raise ValueError("For file inclusion, only strings are supported. In %s" % str(op))
		
		#print "DEBUGGING:"
		#import ipdb; ipdb.set_trace()
		# WITH evaluation.
		# inc_oplist = operations.from_filename(fname).evaluate(reduce_let=False).add_prefix(op.lhs).ops
		# PROBLEM of evaluation in the other context is that we miss all equivalence classes etc.
		inc_oplist = operations.from_filename(fname).add_prefix(op.lhs).ops
		#print "Oplist to include:"
		#pprint.pprint(inc_oplist)
		#print "End of included oplist."
		return inc_oplist

# The registered operation symbols
operation.symbols = { '=': equals, '<=': subsets, '<<': include, '+=': append, ':=': let }

# Regexps for defining the language
lang = namespace()
lang.opor = "|".join(map(re.escape, operation.symbols.keys())) # regex detecting operators
lang.symb_split = ("::", "/")
lang.symb_split_or = "|".join(map(re.escape, lang.symb_split))
lang.symb = r"[a-z_][a-z_0-9:/]*" # regex defining a LHS symbol
lang.comment = r"#" # comment character
lang.linecomment = r"(?:" + lang.comment + r".*)?$" # allows a comment until end of line
lang.stringvarchar = "@" # variable escape character
lang.stringvarsimple = lang.stringvarchar + "([a-z_0-9]+)" # only alphanumeric
lang.stringvarcomplex = lang.stringvarchar + "\{("+lang.symb+")\}"

class sourced:
	"""
	Represents an object (value) which is attached source information to.
	The source information shall tell where it comes from.
	More source information can be appended during the lifetime of this object in
	order to allow chaining of source information (for instance inclusion of files).
	"""
	def __init__(self, value, src):
		self.value = value
		self.sources = [src]
	def __repr__(self):
		return "%s(%s, %s)" % (self.__class__.__name__, self.value, self.sources_as_str())
	def add_source(self, src):
		"Add a source information to this source. src should be a string or sourceline instance."
		self.sources += [src]
	def derive(self, value=None, src=None):
		"Derive a new sourced object which contains a ref to value but a new source list"
		if not value: value = self.value
		src = self.sources if not src else flatten2d([self.sources, [src]])
		return sourced(value, src)
	def sources_as_str(self):
		if len(self.sources)==1:
			return self.sources[0]
		else:
			# TODO: Make this nicer
			return " -> ".join(map(str, self.sources))
		
	# unused?:
	def eval_string_within(self, operations):
		return operations.evalstring(self.value, src)
	
	def __hash__(self): # usable as dict keys
		return hash(self.value)
	def __eq__(self, other): # comparable
		return self.value == other.value
	def __lt__(self, other): # sortable
		return self.value < other.value

class sourceline(namedtuple('sourceline', 'fname linenum text')):
	"""
	Sourceline represents the source of something. It is a key class to allow
	transparent traceback of operations back to their source.
	"""
	
	def __repr__(self):
		return '%s:%s' % (self.fname, self.linenum)
	def verbose(self):
		return '%s line %s: `%s`' % (self.fname, self.linenum, self.text.strip())

	@classmethod
	def from_python(cls):
		"Create a sourceline instance from a python caller position."
		previous_frame = inspect.currentframe().f_back
		(filename, line_number, function_name, lines, index) = inspect.getframeinfo(previous_frame)
		return cls(filename, line_number, lines[0])
	
	@classmethod
	def from_unknown(cls, text="unknown"):
		"Quickly create an instance without any known location"
		return cls('unknown', 0, text)

class sourcemap:
	"""
	Allows to map an iterable with a given source name to inject sourceline objects.
	"""
	def __init__(self,iterable,source_name=None):
		self.iterable = iterable
		# todo: test with StringIO
		if not source_name:
			try:
				self.source_name = iterable.name
			except AttributeError:
				self.source_name = str(iterable) # hopefully short description
		else:
			self.source_name = source_name
	def iter(self):
		"Returns an iterator over the file"
		offset = 1 # python starts with 0 for line counting but humans start with 1
		return (sourceline(self.source_name,linenum+offset,text) for linenum, text in enumerate(self.iterable))
	def map(self,func):
		return map(func, self.iter())

class equivalence_resolver:
	"""
	Base class for the directed and undirected relationship classes
	"""
	def __init__(self):
		self.edges = list()
	def __repr__(self):
		return "%s(%s)" % (self.__class__.__name__, pprint.pformat(self.edges,width=1))
	def add(self, a):
		raise ValueError("please implement")
	def get(self, a):
		raise ValueError("please implement")
	def add_operation(self, op):
		#print "%s.add(%s,%s)" % (self.__class__.__name__, str(op.lhs), str(op.rhs.value))
		assert isa(op.lhs, symbol)
		assert isa(op.rhs.value, symbol)
		return self.add(op.lhs, op.rhs.value)
	def get_resolver(self, a):
		"""
		Get a RHS object for delayed resolving.
		"""
		return equivalence_resolver.delayed_resolver(a, self)
	def get_sourced_resolver(self, sourced_obj):
		if isa(sourced_obj,sourced):
			return sourced(self.get_resolver(sourced_obj.value), sourced_obj.sources)
		else:
			raise ValueError("Got %s, expected sourced() instance. Please use get_resolver instead." % str(sourced_obj))
	class delayed_resolver:
		"This is a functor or future or whatever with a readable repr."
		def __init__(self, a, resolver):
			self.a = a
			self.resolver = resolver
		def __repr__(self):
			return "delayed:%s:resolver(%s)" % (self.resolver.__class__.__name__, pprint.pformat(self.a,width=1))
		def get(self):
			return self.resolver.get(self.a)

class equivalence_classes(equivalence_resolver):
	"""
	Represents equivalence classes between elements. Equal elements are identified by their
	equivalence class. This concrete algorithm here is not in particular fast.
	Naming: We call "classes" = "edges".
	"""
	def add(self, a, b):
		# first, add a new equivalence class:
		self.edges.append({a,b})
		# then cleanup over all cells and merge
		for cx, cy in itertools.combinations(self.edges, 2):
			if cx & cy: # nonempty intersection: merge two equivalence classes
				# print "equal: %s, %s, %s" % (cx,cy,cx&cy)
				self.edges.remove(cx)
				cy.update(cx)
				#print "Updating cy="+str(cy)
		return self # chainable
	def equivalent_nodes(self, a):
		for cls in self.edges:
			if a in cls:
				return cls
		return {a} # equivalence class with its own
	def get(self, a):
		"Return all equivalent nodes except the node itself"
		assert isa(self.equivalent_nodes(a),set), "Broken equivalence_classes: "+str(self)
		return self.equivalent_nodes(a) - {a}


class hierarchy_graph(equivalence_resolver):
	"""
	The directed version of the equivalence_classes.
	Again, this tries to preserve the order of data how they were inserted.
	"""
	def add(self, a, b):
		"insert an edge (a -> b) which we denote as (a <= b)"
		edge = (a,b)
		if not edge in self.edges:
			self.edges.append(edge)
		return self
	def parent(self, a):
		"Returns the list of vertices which point to a"
		return [ b for ai,b in self.edges if ai == a ]
	def ancestors(self, a):
		"""
		Returns the list of all vertices which have a way to a.
		The list includes him own (cf. symbol.ancestors) and doublers in case of a diamant graph form.
		"""
		if not self.parent(a):
			return [a]
		else :
			return [a] + flatten2d([ self.ancestors(ap) for ap in self.parent(a) ])
	def get(self, a):
		"""
		Resolve the graph: Get all ancestors of a, i.e. the list of v where a < v.
		"""
		vertlist = self.ancestors(a)
		vertlist.remove(a) # don't include a itself
		return unique_preserve(vertlist) # remove doublers

class operations:
	"""
	Wraps a list of operations. Basically represents a mexafile.
	This is of type list<operation>.
	"""
	
	def __init__(self, ops=[]):
		self.ops = ops # this is by intention not a copy, but just a reference
		
	def __repr__(self):
		return '%s(%s)' % (self.__class__.__name__, pprint.pformat(self.ops,width=1))
	
	def __iter__(self):
		return iter(self.ops)
	def __len__(self):
		return len(self.ops)
	
	def get_symbols(self):
		"Returns the list of symbols which this operation list holds"
		return [op.lhs for op in self.ops]
	
	def is_empty(self):
		return len(self.ops) == 0
		
	@classmethod
	def from_filehandle(cls, fh):
		# The ordered list of operations as they appear in the file
		return cls(removeNone(sourcemap(fh).map(operation.from_textline)))
		# Prepend the list with the rootbase
		# self.ops = self.rootbase + self.ops # no more.

	@classmethod
	def from_filename(cls, fname):
		"Create an instance from a filename instead of filehandle"
		with open(fname, 'r') as fh:
			return cls.from_filehandle(fh)
	
	@classmethod
	def from_textlines(cls, oplines, source_name="unoriginated textlines"):
		"""
		* List of single lines: Are parsed, source_name used as source name.
		oplines must be iterable (i.e. list of strings or an open file handle).
		"""
		if not isa(oplines, list):
			oplines = [oplines]
		oplist = removeNone(sourcemap(oplines, source_name).map(operation.from_textline))
		return cls(oplist)
	
	def count_optypes(self):
		"""
		Returns a Counter for the kind of optypes which are available in this operations list.
		Useful for debugging and statistics/short output.
		Class names are strings for readability.
		"""
		return collections.Counter([ op.__class__.__name__ for op in self.ops ])
		
	def query(self, root='', exclude_root=False):
		"""
		Get the assignment tree based on some root (which may be string, list, symbol).
		Returns a new operations instance.
		Maintains the order of operations as they appear.
		"""
		root = symbol(root)
		# Search for all symbols which are *below* the root
		# on symbols:
		# tree = { lhs : rhs for lhs,rhs in self.symbols.iteritems() if root in lhs.ancestors() }
		# on the oplist:
		return operations([ op for op in self.ops 
			if root in op.lhs.ancestors() # root="foo", include "foo/bar" and "foo/bar/baz"
			or (root == op.lhs and not exclude_root)  # root="foo", include "foo" itself.
		])
	
	def is_node(self, path, typefilter=[equals,define]):
		"""
		Can decide on a path or operation whether it is a node in the tree.
		Each path is either node, leaf or ill.
		"""
		path_symbs = self.query(path).where_symbol(typefilter).get_symbols()
		return len(path_symbs) > 1 and not symbol(path) in path_symbs
		
	def is_leaf(self, path, typefilter=[equals,define]): # was isLeaf
		"Can decide on a path or operation whether it is a leaf in the tree, i.e. has no more children"
		path_symbs = self.query(path).where_symbol(typefilter).get_symbols()
		return len(path_symbs) == 1 # and symbol(path) in path_symbs
	
	def is_ill(self, path, typefilter=[equals,define]):
		"Something which is neither leaf nor node: c in c=2,c/a=2"
		path_symbs = self.query(path).where_symbol(typefilter).get_symbols()
		return len(path_symbs) > 1 and symbol(path) in path_symbs

	def get_leafs(self):
		"Returns the list of symbols which are leafs"	
		return filter(self.is_leaf, self.get_symbols())
	
	def get_nodes(self):
		"""
		Returns the list of symbols which are nodes.
		The list is ordered in the order the associated operations appear.
		This list also includes the root symbol, symbol().
		"""
		return unique_preserve(flatten2d([ op.lhs.ancestors() for op in self.ops ]))

	def set_symbols(self):
		"""
		Get a dictionary with all assignments. Ensures no assignment appears two times.
		You can choose with values what you want to get as values:
		  * rhs: Just the rhs objects from the operations
		  * index: The index position where the operation (lhs) was (first) defined
		"""
		self.symbols = self.get_dict()

	def get_dict(self, values='rhs', typefilter=[let,equals,define]):
		"""
		was assigndict: Gave only the dictionary of assignments (as a filter).
		"""
		symbols = {} # could use an OrderedDict here if needed.
		for i,op in enumerate(self.ops):
			if isa(op,typefilter):
				if values == 'rhs': val = op.rhs
				elif values == 'index': val = i
				else:
					raise ValueError("Allowed values for 'values' are 'rhs' and 'index', given: '%s'"%values)
				
				# only assign may be used once, let can be overwritten
				if isa(op,equals) and op.lhs in symbols:
					from pprint import pprint #debugging
					pprint(self.ops)
					raise ValueError(
						"Double definition of %s. Was first defined as %s and secondly as %s. List is printed above" % (
							op.lhs,
							symbols[op.lhs],
							op.rhs)
						)
				else:
					symbols[op.lhs] = val
		return symbols
	
	def check_tree_structure(self):
		"""
		Ensure nodes and leafs are distinct groups.
		"""
		for symb in filter(self.is_ill, self.get_symbols()):
			raise ValueError(
				"Symbol %s appears as node but also has children. This is not allowed:\n %s"
				% (symb, self.query(symb))
			)
		
	def check_doublings(self):
		"""
		Check whether some symbols are given multiple times.
		"""
		get_dict()
		
	def add_source(self, fname, line):
		"""
		In-place add a source line to each entry in this operations list.
		The variable names are historic in their meaning.
		"""
		for i,op in enumerate(self.ops):
			op.rhs.add_source(sourceline(fname, i, line))
	
	def evaluate_symbol(self, symbol, inplace=False, eliminate=True, inplace_max_iters=10):
		"""
		Evaluate a questioned symbol, where symbol is a class, for instance `append`.
		Returns new operations object or does the evaluation inplace.
		"""
		ret_oplist = []
		for op in self.ops:
			if isa(op,symbol):
				new_oplist = op.evaluate(self)
				# add/chain backtrace information where the oplist entries comes from
				for i,op in enumerate(new_oplist):
					op.rhs.add_source("evaluation:%s" % symbol.__name__)
				ret_oplist += new_oplist
			else:
				ret_oplist.append(op) # pass throught
		
		if eliminate and any([isa(op,symbol) for op in ret_oplist]):
			# There are still instances of symbol in the oplist and we were
			# asked to eliminate all of them. Recursively call ourselves,
			# expecting that they vanish.
			if inplace_max_iters == 0:
				remaining = operations([op for op in ret_oplist if isa(op,symbol)])
				raise ValueError("While trying to evaluate %s, reached maximum number of iterations. The user probably included cyclic links. These symbols remain: %s" % (str(symbol),str(remaining)))
			return self.evaluate_symbol(symbol, inplace=inplace, eliminate=eliminate, inplace_max_iters=inplace_max_iters-1)
		
		if inplace:
			self.ops = ret_oplist
		else:
			return operations(ret_oplist)

	def evaluate_all_rhs(self):
		"""
		Evaluates all RHS, this works inplace.
		Lists might survive as RHS values. In this case, it is up to the symbol
		evaluation to determine what happens with this list.
		"""
		for op in self.ops:
			op.evaluate_rhs(VariableContext(op,self))

	def evaluate(self, reduce_let=True):
		"Evaluate is chainable and shall always be called after the constructor"
		
		#self.check_tree_structure()

		# eliminate the main symbols
		for symbol in [include, append, equals, subsets]:#, equivalent_node]:
			self.evaluate_symbol(symbol, inplace=True, eliminate=True)
			# this is no more possible after elimination:
			assert not any([isa(op,symbol) for op in self.ops]), "An %s has survived" % symbol.__class__.__name__
			# For early failure:
			# self.check_tree_structure() # checks for tree consistency
			#  -> produces problems in case  a = b
			#                                a/c = 2
			self.set_symbols() # checks for doublings and enables search
		
		# the operations list consists now only of
		#  defines = final primitives and
		#  extends = variables to be resolved.
		
		# next and remaining steps is to replace all the extends.
		### self.evaluate_symbol(extends, inplace=True)

		self.check_tree_structure() # check now instead
		self.symbols = self.set_symbols() # needed for variable access
		self.evaluate_all_rhs()
	
		# evaluate() is chainable:
		return self
	
	def add_prefix(self, path):
		"Return a new operations list where each lhs is prefixed with path"
		if not isa(path,symbol): path = symbol(path)
		return operations([ op.add_prefix(path) for op in self.ops ])
	
	def remove_prefix(self, path):
		"Return a new operation list where each lhs has a common prefix with path removed"
		if not isa(path,symbol): path = symbol(path)
		return operations([ op.remove_prefix(path) for op in self.ops ])
	
	def strip_source(self):
		"""
		Returns a copy where all sourced instances from the RHS are stripped. Useful
		for quickly looking into the data
		"""
		return operations([ op.strip_source() for op in self.ops ])
	
	def where_symbol(self, symbol): # replaces filtertype
		"""
		Returns a new instance where only operations with type symbol are given.
		symbol may be a list of symbols.
		"""
		return operations([op for op in self.ops if isa(op,symbol)])
	
	def is_empty_node(self, varname):
		"""
		Returns whether this is a node without subnodes. May be a leaf, may also be an operation
		which resolves to something once resolved.
		"""
		return self.query(varname, exclude_root=True).is_empty()
	
	def resolve_leaf(self, varname):
		"""
		Resolves a symbol to a leaf. It if is not defined or not a leaf, raises exception.
		Use is_leaf() to check whether it is a leaf before. Returns the RHS.
		"""
		sv = symbol(varname)
		value = [ op for op in self.ops if op.lhs == sv ]
		if len(value) == 1:
			return value[0].rhs
		raise ValueError("Symbol %s is not a leaf but there match these objects: %s" % (str(sv), str(value)))
	
	def resolve_symbol(self, varname, relative_to=symbol(), src=sourceline.from_unknown()):
		"""
		Looks up the value of a variable such as "foo" or "foo/bar" or "foo::bar::baz"
		in the symbols dictionary 'replacements'. In case of errors, src is spilled out.
		
		Note that the return value will *always* be a string. We cast all primitives to
		strings. If the result is *not* a primitive, we print errors.
		Use rhs2python() to convert the string then to a native datatype.
		"""
		sv = symbol(varname)
		if not hasattr(self,'symbols'):
			self.set_symbols()
		if not sv in self.symbols:
			raise ValueError("Variable '%s' not defined but used in %s" % (varname,src.verbose()))
		value = self.symbols[sv].value
		if not isa(value,symbol.primitives):
			raise ValueError("Variable %s value is '%s' and type %s which cannot be inserted at this place. We can only insert a primitive value like %s at this place (source: %s)" % (sv, value, type(value), symbol.primitives),src.verbose())
		return str(value)

	# todo: Move this function where it belongs to
	def get_value_by_key(self, key_name):
		"Returns a RHS value by the symbol key name, stripping source information"
		return self.resolve(self.symbols[symbol(key_name)])
	
	# todo: Move this function where it belongs to
	def resolve(self, irhs):
		""""
		Resolve a RHS object to some plain python object (also lists). For instance 
		   rhs(value=8, src=test.par:86) => 8
		   rhs(value='foobar', src=test.par:84) => 'foobar'
		   [rhs(a),rhs(b),rhs(c)] => [a,b,c]
		"""
		if type(irhs) == list:
			return map(self.resolve, irhs)
		elif isinstance(irhs, sourced):
			return self.resolve(irhs.value)
		else:
			if type(irhs) in symbol.primitives:
				return irhs
			elif isinstance(irhs, symbol):
				# here is the work: Resolve a symbol.
				tree = self.query(irhs)
				if irhs in self.symbols:
					return self.resolve(self.symbols[irhs])
				elif tree:
					# we actually asked for inheritance
					pass
				else:
					raise ValueError("Symbol '%s' not defined" % str(irhs.canonical()))
			else:
				raise ValueError("Don't understand rhs type of %s" % str(irhs))
	
	def outgoing_edges(self, node): # or: neighbours
		# gives all edges starting from node, i.e. if there is find all rel(lhs,rhs) with lhs=node.
		if not isa(node,symbol): node = symbol(node)
		if isa(node,operation): node = node.lhs # as a service
		return operations([ op for op in self if op.lhs == node ])
	
	def paths_from(self, op):
		edges = self.outgoing_edges(op.rhs.value)
		if edges:
			return { op: [ opi for opi in self.paths_from(edges) ] }
		else:
			return [ op ]
	
	
	def paths_from(self, node):
		"""
		Returns all paths starting from a node, i.e. a lhs in real(lhs,rhs).
		Each path is a tuple of operations. A list of paths is returned
		"""
		edges = self.outgoing_edges(node)
		prepend_each = lambda fst, lst: [fst] + lst
		if len(edges):
			return [ tuple(prepend_each(op, self.outgoing_edges(op.rhs.value).ops)) for op in edges ]
		else:
			return []



###
### Output
###

#
# Todo: Clean all the output and representation stuff.
#       We will no more have symbols in the oplist, so complexity reduces massively.
#

class mexafile(operations):
	# TODO: Collect interesting output plotters here, too.
	pass
	
# implemented styles for structured data output (xml, json, yaml)
native_styles = ["tree","tree-backref","linear"]
def native(self, style='linear', root='', symbol_resolver=None):
	"Returns python native data (nested dicts, lists, ...)"
	if   style=='tree':   return self.tree_native  (root,symbol_resolver)
	elif style=='tree-backref': return self.tree_backref_native(root,symbol_resolver)
	elif style=='linear': return self.linear_native(root,symbol_resolver)
	else: raise ValueError("Style: Only linear and tree supported, %s given" % style)

def tree_native(self, root='', symbol_resolver=None, backref=False):
	"""
	Like query, but will result in a nested dictionary structure
	"""
	oplist_absolute = self.query(root) # with absolute paths
	oplist_relative = symbol(root).remove_prefix_oplist(oplist_absolute) # relative paths to root
	# a list with the general structure of the tree, especially parents come before childs
	# so we can make a dict tree out of it
	outline = sorted(unique(flatten2d([op.lhs.ancestors() for op in oplist_relative])))
	
	# setup the tree outline (nodes)
	tree = {}
	for sym in outline:
		subtree = tree
		for symp in sym.path:
			if not symp in subtree:
				subtree[symp] = {}
			subtree = subtree[symp]
	
	# as a placeholder, should be related to a context
	if not symbol_resolver:
		symbol_resolver = lambda sym: "REF="+sym.canonical()
		
	# backref: Include the full path for each node (not leaf)
	backref_key = "$path"
		
	# fill the tree with leafs
	for abs_op,op in izip(oplist_absolute,oplist_relative):
		leaf, rhs = op2tuple(op, symbol_resolver)
		
		if op.lhs.isRoot():
			tree = rhs
		else:
			parent = reduce(dict.get, op.lhs.parent().path, tree) 
			leaf_name = op.lhs.node().canonical()
			parent[leaf_name] = rhs
			
			if backref and not backref_key in parent:
				# we do *not* put the parent in a symbol_resolver but instead
				# use the canonical description.
				# Note that the absolute path is the original one in the defining file,
				# not taking into account the actual mapping (inclusion) of the path.
				parent[backref_key] = abs_op.lhs.parent().canonical()

	return tree

def tree_backref_native(self, root='', symbol_resolver=None):
	"A variant of tree_native with full path references on each node"
	return self.tree_native(root, symbol_resolver, backref=True)

def linear_native(self, root='', symbol_resolver=None):
	if not symbol_resolver:
		symbol_resolver = lambda sym: "REF="+sym.canonical()
	return dict(map(lambda op: op2tuple(op, symbol_resolver), self.query(root)))

def dump_mexa(self, root='', simple_mexa=False, return_list=False, opfilter=idfunc):
	"""
	Print the oplist in a form which is again in the mexa format, is the most compact
	form (no redundancies, no copies, instead still references). In contrast to the input document,
	1) all overwritable definitions have been resolved
	2) all strings have been evaluated
	3) the file is checked for consistency
	4) order is unchanged
	5) data are normalized (boolean, also numbers)
	6) no more comments or newlines
	"""
	ret = []
	
	#rec_simplemexa = lambda root: self.dump_mexa(root=root, simple_mexa=True, return_list=True)#, opfilter=symbol(root).prefix_oplist)
	# kind of functools.partial
	#op2strMaker = lambda symbolmapper: lambda op: op2str(op, symbolmapper)

	for op in opfilter(self.query(root)):
	#	symbolmapper = rec_simplemexa if simple_mexa else symbol.canonical
	#	myop2str = lambda op: op2str(op, symbolmapper) 
		
		# assert isa(op,assign), "There is a "+str(op)
		if isa(op.rhs, list):
			# reconstruct the append operation from the list assignment
			if simple_mexa:
				# resolve the list to it's meaning.
				for irhs in op.rhs:
					assert isa(irhs.value, symbol), "We only support the += operator on symbols: "+str(irhs)
					ret += self.dump_mexa(root=irhs.value, simple_mexa=True, return_list=True, opfilter=op.lhs.prefix_oplist)
			else:
				local_oplist = [ append(op.lhs,irhs) for irhs in op.rhs ]
				ret += map(op2str, local_oplist)
		elif isa(op.rhs.value,symbol) and simple_mexa:
			opfilter = lambda oplist: op.lhs.prefix_oplist(op.rhs.value.remove_prefix_oplist(oplist))
			ret += self.dump_mexa(root=op.rhs.value, simple_mexa=True, return_list=True, opfilter=opfilter)
		else:
			ret.append(op2str(op))
		
		

	trailing_newline = [""] # to get \n the last character in the output
	return ret if return_list else "\n".join(ret+trailing_newline)

def simple_mexa(self, root=''):
	"""
	Simple mexa: Only assignments, nothing else. Redundancy will happen, thought.
	"""
	return self.dump_mexa(root,simple_mexa=True)

encodings = {
	# base64 is the well known base64
	'base64': lambda f: base64.b64encode(f),
	# base16 is only the hex
	'base16': lambda f: base64.b16encode(f),
	# urlencode/quoted printable:
	'quotedprintable': lambda f: quoted_printable(f)
}
encoding_header = "##MEXA-simple configuration file"
def dump_encoded_mexa(self, encoding, prepend_header=encoding_header, root='', opfilter=idfunc):
	"""
	Return an encoded mexa file. Several encodements are available, see the
	encodings table.
	@arg encoding: A value in encodings
	@arg prepend_header: A single header line (incl. comment sign) to be prepended before the file
	@arg simple_mexa: Whether to create simple mexa output
	"""
	fcontent = self.dump_mexa(root=root, simple_mexa=True, opfilter=opfilter, return_list=False)
	fcontent = prepend_header + "\n" + fcontent
	return self.encodings[encoding](fcontent)

def json(self, style='linear', root=''):
	"""
	Translate the oplist to a more human readable one, omitting all the classy
	meta data. In Json without hierarchy
	"""
	
	# JSON pointer reference
	json_pointer = lambda irhs: { '$ref': '#/'+irhs.canonical() }
	data = self.native(style, root, json_pointer)
	
	import json
	return json.dumps(data)

def yaml(self, style='linear', root=''):
	# placeholder
	yaml_pointer = lambda irhs: { 'YAML-REF': 'towards->'+irhs.canonical() }
	data = self.native(style, root, yaml_pointer)
	
	# this has to be installed
	import yaml
	return yaml.dump(data,default_flow_style=False)

def xml(self, style='linear', query_root=''):
	# xml dumps straight the classy structure
	# This is extremely verbose and more suitable as a tech-demo
	# Could also implement a hierarchy writer here
	
	# lxml and xml.sax are included in python (batteries)
	from lxml import etree
	from xml.sax.saxutils import escape as xml_escape
	
	root = etree.Element(self.__class__.__name__)
	root.set('style', style) # actually, the style is ignored here.
	root.set('query_root', query_root)
	
	def visit(obj, parent, tagname=None):
		if not tagname: tagname = obj.__class__.__name__
		me = etree.SubElement(parent, tagname)
		
		if type(obj) == list:
			for oi in obj:
				visit(oi, me)
		elif hasattr(obj, '_asdict'): # namedtuples instable API gives OrderedDict
			for namedtuple_fieldname, value in obj._asdict().iteritems():
				visit(value, etree.SubElement(me, namedtuple_fieldname))
		elif isinstance(obj, symbol):
			for p in obj.path:
				etree.SubElement(me, 'pathseg').text = p
		else:
			me.text = xml_escape(str(obj)).strip()

	visit(self.query(query_root), root)
	return etree.tostring(root, pretty_print=True)

#def exaspecfile(self, root="exahype", tplfile='./jinja-test.txt'):#'./exa-specfile-tpl.exahype'):
def exaspecfile(self, root="exahype", embed_config=True, tplfile='./exa-specfile-tpl.exahype'):
	"Return an ExaHyPE specfile which corresponds to these data"
	
	copy_resolver = lambda symb: self.tree_backref_native(root=symb,symbol_resolver=copy_resolver)
	ctx = self.tree_backref_native(root,copy_resolver)
	
	# assuming we only have one exahype project in ctx
	ctx = {'exahype_projects': [ctx] }
	
	# replace True and False by "on" and "off"
	exaBool = { True: "on", False: "off" }
	replBool = lambda item: exaBool[item] if isa(item,bool) else item
	ctx = mapComplex(replBool, ctx)

	# for debugging:
	if False:
		pprint.pprint(ctx)
	
	mexa_path = mf.get_value_by_key("mexa")
	
	# this has to be installed:
	import jinja2
	
	jinja_env = jinja2.Environment(
		loader=jinja2.FileSystemLoader(mexa_path),
		undefined=jinja2.StrictUndefined
	)
	
	# provide further jinja functions:
	# w decorators: https://stackoverflow.com/a/47291097
	
	def jinja_filter(func):
		jinja_env.filters[func.__name__] = func
		return func
	
	@jinja_filter
	def dimlist(comp_domain, field):
		"Compute the ExaHypE computational domain string (width_x, width_y, width_z?) "
		fieldnames = [field+"_"+i for i in "xyz"]
		fieldnames = fieldnames[:comp_domain['dimension']]
		fields = [str(comp_domain[fn]) for fn in fieldnames]
		return ", ".join(fields)
	
	@jinja_filter
	def count_variables(variable_list):
		"""Count the variables in a list 'x,y,z' or 'x:1,y:5,z:17' to 3 or 23, respectively"""
		matches = re.findall(r"([a-zA-Z-_]+)(?:[:](\d+))?", variable_list)
		return sum([1 if count == "" else int(count) for label,count in matches ])
	
	@jinja_filter
	def resolve_path(node, prepend='/'):
		"""
		Lookup the tree_backref_native inserted backreference and embed it.
		In order to have this be understood as a file path by the ExaHyPE specfile parser,
		prepend a slash.
		"""
		# print out the query xpath (backref)
		backref_key = "$path"
		if backref_key in node:
			return prepend + node[backref_key]
		else:
			raise ValueError("Missing backref key '%s' in node '%s'" % (backref_key, str(node)))
		
	@jinja_filter
	def embed(node, encoding):
		# embed the configuration as a string, suitable for the specfile constants
		# parameters
		root = resolve_path(node)
		# embedded: Remove the prefix since we do not need the fully featured address list
		opfilter=symbol(root).remove_prefix_oplist # => does not work, opfilter is broken
		#opfilter = idfunc
		return self.dump_encoded_mexa(encoding=encoding, root=root, opfilter=opfilter)
	
	@jinja_filter
	def link_in_list(node):
		# This can be either resolve_path or embed. Here we do both for simplicity
		encoding = 'quotedprintable'
		return 'mexa:embedded,mexaref:%s,mexaformat:%s,mexacontent:%s' % (resolve_path(node), encoding, embed(node,encoding) )
	
	@jinja_filter
	def as_float(txt):
		"""
		Ensure a number does look like a float (2.) instead of an int (2).
		This is what the ExaHyPE parser is sensitive on...
		"""
		return str(float(txt))
	
	# function (not a filter):
	jinja_env.globals['error'] = lambda msg: raise_exception(ValueError("Template stopped: "+msg))
	
	#try:
	return jinja_env.get_template(tplfile).render(ctx)
	# without exception chaining, we loose the stack trace:
	#except Exception as e:
	#	print "Debuggin context:"
	#	pprint.pprint(ctx)
	#	print "Exception:"
	#	print str(e)
	#	raise e	

encode_languages=["json", "yaml", "xml"]
def encode(self, lang='json', style='linear', root=''):
	"Encode in some language"
	if   lang=="json": return self.json(style,root)
	elif lang=="yaml": return self.yaml(style,root)
	elif lang=="xml":  return self.xml (style,root)
	else: raise ValueError("Language %s not supported" % lang)
	

# TODO: Add a graphviz output option.

# a basic frontend:

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('expressions', type=str, nargs='*',# action='append',
	metavar='foo=bar', help="Mexa expressions to validate. Each argument represents one line. Arguments are appended to files read in.")
parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
	metavar='FILENAME', help="Input file to read. If no file is given, read from stdin.")
parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
	metavar='FILENAME', help="Output file to write to. If no file is given, write to stdout.")
parser.add_argument('--outformat', choices=['mexa','simple-mexa','encode-mexa','yaml','json','xml','exahype'], 
	help="Language format for output", default='mexa')
parser.add_argument('--style', choices=native_styles, default='linear', help="Style applied if outformat in "+str(encode_languages))
parser.add_argument('--root', default='', help='Queried root container')
args = parser.parse_args()

mf = mexafile.from_filehandle(args.infile)

append_oplist = removeNone(args.expressions)
if len(append_oplist):
	#print "Parsing: " + str(append_oplist)
	mf.append(append_oplist, "Command line argument expression")

# could do here also --env which then adds all or requested environment variables below env/

# appending stuff to any file, as a root
mf.ops += [
	# mexa = path to the current script directory
	let(symbol("mexa"), sourced(os.path.dirname(os.path.realpath(__file__)), sourceline.from_python())),
]

sys.exit() # stop here for the time being

# todo: Could offer here also an option to inject various or all environment
# variables, just by passing "--env". Could setup manual assign() operations
# similar to mexafile.rootbase.
# Could also get rid of "let" statements by moving the rootbase outside the class.

mf.evaluate()

write = lambda t: args.outfile.write(t)

# 1. Create list of operations
# mf = mexafile(open("test.par","r"))
#for o in mf.oplist: print o

if args.outformat in mexafile.encode_languages:
	write(mf.encode(lang=args.outformat, style=args.style, root=args.root))
elif args.outformat == 'mexa':
	write(mf.dump_mexa(root=args.root))
elif args.outformat == 'simple-mexa':
	write(mf.simple_mexa(root=args.root))
elif args.outformat == 'encode-mexa':
	write(mf.dump_encoded_mexa('quoted_printable',root=args.root))
elif args.outformat == 'exahype':
	# could pass root, but then must ensure that exahype information exist, i.e.
	# at least one exahype project or so.
	write(mf.exaspecfile())
else:
	write("Argformat %s not yet implemented" % args.outformat)
