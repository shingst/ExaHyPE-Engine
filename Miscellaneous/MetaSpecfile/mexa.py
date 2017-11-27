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
import os, re, sys, ast, inspect, argparse, base64, pprint
from argparse import Namespace
from collections import namedtuple
from itertools import izip
#from future.utils import raise_from # no batteries

# helpers and shorthands:
baseflags = re.IGNORECASE
match = lambda pattern,string,flags=0: re.match(pattern,string,baseflags+flags)
unpack = lambda f: lambda p: f(*p)
# The identity function. Don't mix up with pythons internal "id"
idfunc = lambda x: x
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

# operations.
let = namedtuple('let', 'lhs rhs') # overwritable = can be overwritten
assign = namedtuple('assign', 'lhs rhs') # immutable = set only once
append = namedtuple('append', 'lhs rhs')
extend = namedtuple('extend', 'lhs rhs') # extend datastructure
include = namedtuple('include', 'lhs rhs') # include file
opsymbol = { '=': assign, '<=': extend, '<<': include, '+=': append, ':=': let }

# reduce operation (no more lhs):
rhs = namedtuple('rhs', 'value src')
# -> should collect: hasnosymbol, parseRHS

lang = Namespace()
lang.opor = "|".join(map(re.escape, opsymbol.keys())) # regex detecting operators
lang.symb = r"[a-z_][a-z_0-9:/-]*" # regex defining a LHS symbol
lang.comment = r"#" # comment character
lang.linecomment = r"(?:" + lang.comment + r".*)?$" # allows a comment until end of line
lang.stringvarchar = "@" # variable escape character
lang.stringvarsimple = lang.stringvarchar + "([a-z_0-9]+)" # only alphanumeric
lang.stringvarcomplex = lang.stringvarchar + "\{("+lang.symb+")\}"

# atoms:
class symbol:
	# primitive types which we allow to store. We also allow lists of
	# these types
	primitives = [int,float,str,unicode,bool]

	def __init__(self,name):
		""""
		Creates a symbol from a name such as Foo/Bar/Baz or a path such
		as ['Foo','Bar','Baz'].
		"""
		if type(name) == list:
			self.path = map(str.lower, map(str.strip, name))
		elif isinstance(name, symbol):
			self.path = name.path # kind of copy constructor
		else:
			name = str(name).lower() # in order to get case-insensitive
			self.path = map(str.strip, re.split('::|/', name))
		# in case of roots and merged paths and so on
		self.path = removeFalse(self.path)
		
	def canonical(self):
		"Canonical string representation with a specific seperator"
		return "/".join(self.path)
	def __repr__(self):
		return "symbol(%s)" % self.canonical()
	def isRoot(self):
		return len(self.path) == 0
	
	def node(self):
		if self.isRoot():
			return self
		return symbol(self.path[-1])
	def parent(self):
		if self.isRoot():
			raise ValueError("Symbol is already root")
		return symbol(self.path[:-1])
	def ancestors(self):
		"""
		Lists all parenting symbols, for instance for /a/b/c it is the
		list [/, /a, /a/b].
		"""
		return [symbol(self.path[:i]) for i,_ in enumerate(self.path)]
	
	def prefix_symbol(self, other):
		"Returns a new symbol wich is prefixed by self."
		return symbol(self.path + other.path)

	def prefix_oplist(self, oplist):
		"Prefix each lhs on the oplist with "
		return [ type(op)(self.prefix_symbol(op.lhs), op.rhs) for op in oplist ]
	
	def remove_prefix_oplist(self, oplist):
		"Remove the prefix of self on every item in oplist"
		return [ type(op)(symbol(remove_comstr_prefix(op.lhs.canonical(), self.canonical())), op.rhs) for op in oplist]
	
	# allow symbols to be dict keys:
	def __hash__(self):
		return hash(tuple(self.path))
	def __eq__(self, other):
		return self.canonical() == other.canonical()
	def __lt__(self, other): # sortable
		return self.canonical() < other.canonical()
	
def hasnosymbol(rhs):
	if type(rhs) in symbol.primitives:
		return True
	if type(rhs) == list:
		return all(map(iscomplete,rhs))
	if type(rhs) == dict:
		raise ValueError("Dictionaries are not supported as right hand sides in Operation %s" % str(self))

class operation:
	def __init__(self, lhs, rhs, srcmap):
		self.lhs, self.rhs = lhs, rhs
		self.srcmap = srcmap
	def __repr__(self):
		return "%s(%s, %s)" % (type(self), lhs,rhs)
	def iscomplete(self):
		return hasnosymbol(self.rhs)
	def eval(self):
		raise NotImplementedError("Subclasses should implement this")
		
# allow transparent traceback of operations back to their source
class sourceline(namedtuple('sourceline', 'fname linenum text')):
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


# First step: Parse parfile into a list of operations

def parseRHS(text):
	text = text.strip()
	# We use python's parser for both understanding the data structure
	# and removing pythonic line comments such as "#".
	try:
		typed = ast.literal_eval(text)
		if type(typed) in symbol.primitives:
			return typed
		else:
			# we probably have a created a weird type, for instance when
			# parsed an object such as [1,2,3] which we don't want to support.
			# Instead, thus, return a text. This may not correctly strip comments,
			# thought.
			return text
	except (ValueError, SyntaxError):
		# This is most likely a string or something weird
		
		if re.match(r"^(Yes|True|On)"+lang.linecomment, text, re.IGNORECASE):
			return True
		if re.match(r"^(No|False|Off)"+lang.linecomment, text, re.IGNORECASE):
			return False
		
		# probably a link to something or several things. This *always* yields
		# a list.
		symbs = re.match(r"^("+lang.symb+")(\s+"+lang.symb+")*"+lang.linecomment, text, re.IGNORECASE)
		if symbs:
		      return unlistOne(map(symbol, removeNone(symbs.groups())))
		
		return text
	return text

def line2operation(srcline):
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
	value = parseRHS(parts.group('rhs'))

	return opsymbol[parts.group('op')](name, rhs(value, srcline))

def oplist2assigndict(oplist, values='rhs'):
	symbols = {} # could use an OrderedDict here if needed.
	for i,op in enumerate(oplist):
		if isa(op,[let, assign]):
			if values == 'rhs': val = op.rhs
			elif values == 'index': val = i
			else:
				raise ValueError("Allowed values for 'values' are 'rhs' and 'index', given: '%s'"%values)
			
			# only assign may be used once, let can be overwritten
			if isa(op,assign) and op.lhs in symbols:
				from pprint import pprint #debugging
				pprint(oplist)
				raise ValueError(
					"Double definition of %s. Was first defined as %s and secondly as %s. List is printed above" % (
						op.lhs,
						symbols[op.lhs],
						op.rhs)
					)
			else:
				symbols[op.lhs] = val
	return symbols

def op2tuple(op, symbolmapper):
	""""
	Extract an op to a plain python object, stripping the objects and relation in it
	@args symbolmapper a Function which maps a symbol to something else
	"""
	def rhsValue2py(irhs):
		if   isa(irhs, rhs): # unwrapping
			return rhsValue2py(irhs.value)
		elif isa(irhs, list): # threading
			return map(rhsValue2py, irhs)
		elif isa(irhs, symbol):
			return symbolmapper(irhs)
		else:
			return irhs
		
	fst = op.lhs.canonical()
	if isa(op,assign) and isa(op.rhs,list):
		snd = map(rhsValue2py, op.rhs)
	else:
		snd = rhsValue2py(op.rhs.value)
	return (fst,snd)

def op2str(op, symbolmapper=None):
	"""
	Converts an operation to a string, looking similar to the parsed input of the operator
	In contrast to op2tuple (which strips the operator), this one does not handle lists on the
	rhs. It expects them to be removed at some calling step, if neccessary. Otherwise they just
	get strings.
	"""
	if not symbolmapper: symbolmapper = symbol.canonical
	symbop = invdict(opsymbol)
	s_lhs = op.lhs.canonical()
	s_op = symbop[type(op)]
	if isinstance(op.rhs.value, symbol):
		s_rhs = symbolmapper(op.rhs.value)
	else:
		delim = '"' if type(op.rhs.value) in [str,unicode] else ''
		s_rhs = delim+str(op.rhs.value)+delim
	return "%s %s %s" % (s_lhs, s_op, s_rhs)

def quoted_printable(s, escape='%'):
	"""
	Returns a quoted printable version of the string s with escape character %, no maximum line length
	(i.e. everything in a single line) and everything non-alphanumeric replaced.
	"""
	# sourcecode inspired by quopri, https://github.com/python/cpython/blob/2.7/Lib/quopri.py
	HEX = '0123456789ABCDEF'
	def quote(c):
		"""Quote a single character."""
		i = ord(c)
		return escape + HEX[i//16] + HEX[i%16]
	def needsquote(c):
		"""Whether character needs to be quoted"""
		return not ('0' <= c <= '9' or 'a' <= c <= 'z' or 'A' <= c <= 'Z')
	return "".join([ quote(c) if needsquote(c) else c for c in s])


class mexafile:
	# These operations are prepended to any read in file.
	rootbase = [
		# mexa = path to the current script directory
		let(symbol("mexa"), rhs(os.path.dirname(os.path.realpath(__file__)), sourceline.from_python())),
	]
	
	def __init__(self, fh):
		# The ordered list of operations as they appear in the file
		self.oplist = removeNone(sourcemap(fh).map(line2operation))
		# Prepend the list with the rootbase
		self.oplist = self.rootbase + self.oplist
		# run the actual operator evaluator
		# self.evaluate() ## shall be invoked manually
		
	def append(self, oplines, source_name="unkown"):
		"""
		Quickly add operations (list of single lines) into this oplist.
		"""
		if not isa(oplines, list):
			oplines = [oplines]
		self.oplist += sourcemap(oplines, source_name).map(line2operation)

	@classmethod
	def from_filename(cls, fname):
		"Create an instance from a filename instead of filehandle"
		with open(fname, 'r') as fh:
			return cls(fh)
		
	def getsymbols(self, values='rhs'):
		"""
		Get a dictionary with all assignments. Ensures no assignment appears two times.
		You can choose with values what you want to get as values:
		  * rhs: Just the rhs objects from the operations
		  * index: The index position where the operation (lhs) was (first) defined
		"""
		return oplist2assigndict(self.oplist, values)
		
	def evaluate(self, reduce_let=True):
		"Evaluate is chainable and shall always be called after the constructor"
		
		# 2. Preliminary symbol dictionary, for string expansion in includes.
		self.symbols = self.getsymbols()

		# 3. Evaluate the includes recursively in-place:
		while len(self.filtertype(include)):
			for i,op in enumerate(self.oplist):
				if isa(op,include):
					# just ensure here that inheritance works correctly
					if not isa(op.rhs.value, [str,unicode]):
						raise ValueError("For file inclusion, only strings are supported. In %s" % str(op))
					
					# include a file
					fname = self.evalrhsstring(op.rhs)
					inc_oplist = mexafile.from_filename(fname).evaluate(reduce_let=False).oplist
					#print "Oplist to include:"
					#pprint.pprint(inc_oplist)
					#print "End of included oplist."
					# allow including at any point
					inc_oplist = op.lhs.prefix_oplist(inc_oplist)
					self.oplist = replace_with_list_at(i, self.oplist, inc_oplist)
					#print "Included in place:"
					#pprint.pprint(self.oplist)
					#print "Done"
					break # go to next iteration of "which has_includes"
		
		# 4. Evaluate data structure extensions in-place:
		while len(self.filtertype(extend)):
			self.symbols = self.getsymbols() # update for the query search
			for i,op in enumerate(self.oplist):
				if isa(op,extend):
					# ensure that we extend only from symbols
					if not isa(op.rhs.value, symbol):
						raise ValueError("We only can extend from other symbols. In %s" % str(op))
					
					# include another tree structure
					queried_symbol = op.rhs.value
					ext_oplist = self.query(queried_symbol)
					ext_oplist = op.lhs.prefix_oplist(queried_symbol.remove_prefix_oplist(ext_oplist))
					self.oplist = replace_with_list_at(i, self.oplist, ext_oplist)
					break
		
		# 4. Update symbol dictionary for string evaluation
		self.symbols = self.getsymbols()

		# Shorthands to update the RHS
		def updateRhs(i, newrhs):
			"Since we have immutable tuples, we need to create new ones"
			optype = type(self.oplist[i])
			self.oplist[i] = optype(self.oplist[i].lhs, newrhs)
		def updateOplistRHS(updater):
			for i,op in enumerate(self.oplist):
				if(isa(op.rhs,list)):
					newrhs = map(updater, op.rhs)
					if all(newrhs): # check for None
						updateRhs(i,newrhs)
				elif(isa(op.rhs,rhs)):
					newrhs = updater(op.rhs)
					if newrhs:
						updateRhs(i, newrhs)
				else:
					raise ValueError("Unsupported type for %s in operation %s" % (op.rhs,op))
		
		# 5. Evaluate all strings.
		def evaluateRhsString(irhs):
			if isa(irhs.value,[str,unicode]):
				newval = self.evalstring(irhs.value, irhs.src)
				newrhs = rhs(newval, irhs.src)
				#if newval != irhs.value:
				#	print newrhs
				return newrhs
			else:	return None # for clarity
		updateOplistRHS(evaluateRhsString)

		# 6. Collect the appends and join them to the list of assignments,
		#    preserving the order of the definitions.
		# Create a new list in order to avoid modifying while looping.
		new_oplist = []
		for i,op in enumerate(self.oplist):
			if isa(op,append):
				# they change due to the subsequent oplist modifications
				#opindices = self.getsymbols('index') # returns only assign indices
				opindices = oplist2assigndict(new_oplist, 'index') # returns only assign indices

				# find the first assignment occurance in the file
				if op.lhs in opindices.keys():
					opindex = opindices[op.lhs]
					# update it
					if isa(new_oplist[opindex].rhs,list):
						new_oplist[opindex].rhs.append(op.rhs)
					else:
						old = new_oplist[opindex]
						new_oplist[opindex] = assign(old.lhs, [old.rhs, op.rhs])
					# remove the append operation
					#assert i != opindex
					#del self.oplist[i]
				else:
					# there is not yet/even an assignment operation.
					# Replace the append with an assign.
					#opindex = i
					#self.oplist[opindex] = assign(op.lhs, [ op.rhs ])
					newop = assign(op.lhs, [op.rhs])
					#print "Making new "+str(newop)
					new_oplist.append(newop)
			else:
				new_oplist.append(op)
		# TODO: Clean the comments and leftovers from the old inplace algorithm
		self.oplist = new_oplist
		#pprint.pprint(self.oplist); print "after step 6"
		
		# now all appends should be removed from the list
		for op in self.oplist:
			assert not isa(op,append), "An append has survived"
	
		# once again, update the symbols dict
		self.symbols = self.getsymbols()

		# Old version of evaluating append operations, but on the symbols dict instead
		# on the ordered list.
		if False:
			for op in self.filtertype(append):
				if op.lhs in self.symbols:
					if type(self.symbols[op.lhs]) == list:
						self.symbols[op.lhs].append(op.rhs)
					else:
						self.symbols[op.lhs] = [ self.symbols[op.lhs], op.rhs ]
				else:
					self.symbols[op.lhs] = [ op.rhs ]
		
		# 7. Resolve inheritance and resolve any symbols on the RHS
		#def evaluateRhsSymbols(irhs):
		#	return irhs
		#updateOplistRHS(evaluateRhsSymbols)
		# --> this is instead done getplain().
		
		# 8. Could also ensure that the configuration clearly seperates into
		#    Nodes  (holding subnodes)
		#    Leafs  (holding RHS data)
		
		# 9. Reduce let statements.
		if reduce_let:
			#[ assign(lhs, rhs) for lhs, rhs in self.symbols ]
			#for op in self.filtertype(let):
			lets, lets_positions = self.filtertype(let, return_indices=True)
			lets_names = unique([ op.lhs for op in lets ])
			# delete all the lets
			# self.oplist = without_indices(self.oplist, lets_positions)
			# -> done anyway already below.
			# and prepend them as assignments before everything else
			lets_values = [ self.symbols[lhs] for lhs in lets_names ]
			resolved_lets = [ assign(lhs, irhs) for lhs, irhs in izip(lets_names, lets_values) ] 
			self.oplist = resolved_lets + [op for op in self.oplist if not isa(op,let) ]
			# Could improve this by positioning the assign() in place of the last occurance of
			# a let() with same lhs.
			
		# evaluate() is chainable:
		return self

	def filtertype(self, optype, return_indices=False):
		"Typical optypes are assign or include. optype can be a list of optypes"
		items = [op for op in self.oplist if isa(op,optype)]
		indices = [i for i,op in enumerate(self.oplist) if isa(op,optype)]
		return (items,indices) if return_indices else items
	
	def getvar(self, varname, src=sourceline.from_unknown()):
		"""
		Looks up the value of a variable such as "foo" or "foo/bar" or "foo::bar::baz"
		in the symbols dictionary 'replacements'. In case of errors, src is spilled out.
		Note that the return value will *always* be a string. We cast all primitives to
		strings. If the result is *not* a primitive, we print errors.
		"""
		sv = symbol(varname)
		if not sv in self.symbols:
			raise ValueError("Variable '%s' not defined but used in %s" % (varname,src.verbose()))
		value = self.symbols[sv].value
		if not type(value) in symbol.primitives:
			raise ValueError("Variable %s value is '%s' and type %s which cannot be inserted at this place. We can only insert a primitive value like %s at this place." % (sv, value, type(value), symbol.primitives))
		return str(value)

	def evalstring(self, text, src=sourceline.from_unknown()):
		"""
		Evaluate strings such as "bla @foo @{bar}baz".
		@src shall be a sourceline object if problems occur
		"""
		# find only all variables:
		#vars = re.findall(lang.stringvarsimple, text, baseflags)
		#vars = re.findall(lang.stringvarcomplex, text, baseflags)
		
		# replace the variables:
		replmatch = lambda matchobj: self.getvar(matchobj.group(1), src)
		text = re.sub(lang.stringvarsimple, replmatch, text, count=0, flags=baseflags)
		text = re.sub(lang.stringvarcomplex, replmatch, text, count=0, flags=baseflags)
		return text
	
	def evalrhsstring(self, rhs):
		return self.evalstring(rhs.value, rhs.src)
	
	def get_value_by_key(self, key_name):
		return self.resolve(self.symbols[symbol(key_name)])
	
	def resolve(self, irhs):
		""""
		Resolve a RHS object to some plain python object (also lists). For instance 
		   rhs(value=8, src=test.par:86) => 8
		   rhs(value='foobar', src=test.par:84) => 'foobar'
		   [rhs(a),rhs(b),rhs(c)] => [a,b,c]
		"""
		if type(irhs) == list:
			return map(self.resolve, irhs)
		elif isinstance(irhs, rhs):
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
	
	def query(self, root=''):
		"""
		Get the assignment tree based on some root (which may be string, list, symbol)
		"""
		root = symbol(root)
		# Search for all symbols which are *below* the root
		# on symbols:
		# tree = { lhs : rhs for lhs,rhs in self.symbols.iteritems() if root in lhs.ancestors() }
		# on the oplist:
		return [ op for op in self.oplist 
			if root in op.lhs.ancestors() # root="foo", include "foo/bar" and "foo/bar/baz"
			or root == op.lhs             # root="foo", include "foo" itself.
		]
	
	def isLeaf(self, path):
		"Leaf: Has no more children"
		return bool(self.query(path))
	
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
parser.add_argument('--style', choices=mexafile.native_styles, default='linear', help="Style applied if outformat in "+str(mexafile.encode_languages))
parser.add_argument('--root', default='', help='Queried root container')
args = parser.parse_args()

mf = mexafile(args.infile)

append_oplist = removeNone(args.expressions)
if len(append_oplist):
	#print "Parsing: " + str(append_oplist)
	mf.append(append_oplist, "Command line argument expression")

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
