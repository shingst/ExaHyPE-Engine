# The Mexa Language

We start inverted with Level 3 down to Level 1.

## Level 3: Hierarchical dictionaries

We hereby give a cartoon how a hierarchical mapping of parameters
could look like:

```
.
├── Name: "AdvectionTest"
├── domain
│   ├── size
│   │   ├── width:  10
│   │   └── height: 10
│   └── grid
│       ├── dx
│       │   ├── min:  1e-3
│       │   └── max:  0.5
│       └── dy
│           ├── min:  1e-3
│           └── max:  0.5
├── boundary_conditions
│   ├── left: "exact"
│   ├── right: "exact"
│   ├── top: "exact"
│   └── bottom: "exact"
└── initial_data
    ├── name: "BlackHole"
    └── bhcode
        ├── pos
        │   ├── x:  0
        │   └── y:  0
        ├── M: 1.0
        └── a: 0.9
```

There are a number of obvious facts in such a tree of assignments:

* There are either directories/folders/nodes or data/leafs. One cannot be both.
* This is a graph which is a tree with a root, here denoted with a point.
* No key may appear two times in some directory.
* There is no such thing as a "plain list" (think of the black hole position vector)
* The order of entries might be relevant or not, we cannot say.
* Dictionaries appear as universal concept here. They either map to (nested) dictioniaries or primitive data.
* Thanks to the unique naming, we can address every item with an absolute path.

Different data types can be represented by this data structure. For instance,
ordered or unordered lists can just display as alphabetically ordered dicts. The
list `[1,3,7]` is thus represented by `{i0:1,i1:3,i2:7}` or `{x:1,y:3,z:7}`.
In case an ordering is not imporatnt, the keys are completely arbitrary.

Similarly, one can think of multiline strings as lists and subsequently as maps.
The mutliline message `foo\nbar\nbaz` then gets `{s1:'foo',s2:'bar,s3:'baz'}`.

Obviously, being a data structure, there is nothing beyond the data, such as coments.
Instead, everything is data and comments are *meta data*. In the exa approach, there
is no schema (c.f. the documentation file `CONCEPTS.md`). We think of this as
*there cannot be too much data*. Any tree could be enriched with keys such as
`{_desc:'an unordered list', _date:'2018-01-27', u:1, z:7, d:3}`. In this paradigm,
more data is always better then less data.

## Level 2: Simple-Mexa

The `simple-mexa` language is a simple and straightforward way to represent the given
information from the above example. It could be written down as

```
# This is a single line comment
Name = "AdvectionTest"
domain/size/width = 10  # comments can be anywhere
domain/size/height = 10
grid/dx/min = 1e-3
grid/dx/max = 0.5
grid/dy/min = 1e-3
grid/dy/max = 0.5
boundary_conditions/left = "exact"
boundary_conditions/right = "exact"
boundary_conditions/top = "exact"
boundary_conditions/bottom = "exact"
initial_data/name = "BlackHole"
initial_data/bhcode/pos/x = 0
initial_data/bhcode/pos/y = 0
initial_data/bhcode/M = 1.0
initial_data/bhcode/a = 0.9
```

The following rules apply:

 * Every command is always a single line.
 * All assignments are immutable (pi=3.14 is always 3.14)
 * The order of the lines is interchangable.
 * Variables can be sorted hierarchically.
 * Native support for Int, Float, String.
 * Full path verbosity implements principle of least surprise.
 
Furthermore, note that due to the linearization of the tree, an order is imposed
and simple-mexa files thus can represent sorted dictioniaries.

Keys in Mexa can be any strings, but the operator symbol `=` is reserved as well
as the dictionary seperator. It might not only be `/` but also `::` is allowed.

### Primitive types in Simple-Mexa

On the right hand side, the data types Integer, floating point number, boolean
and String are accepted. These are refered to as *primitives*. All numbers can be
given with leading sign, `+0` is equally correct as `+2.3e+7`. Booleans are
indicated as keywords, `True / Yes / On` means logical one while
`False / Off / No` means logical zero.

Strings have to be enclosed in single quotes or double quotes. This has several
advantages: They can be distinguished from booleans, comment characters can occur
in strings and leading and trailing whitespace is clearly defined. One could think
of a variant of simple-mexa with non-enclosed strings. Such variants exist, but
they are ill-defined for these very reasons.

Note that multiline strings are not supported by the language as first class citizens. 
Instead, they should be represented by dictionaries, as proposed in the level 3:
  
```
another/file/l0 = "this could be the"
another/file/l1 = "content of another"
another/file/l2 = "file"
```

## Level 1: Functional mexa

The `full-mexa` is an extension of the simple mexa file format. Every simple mexa
file is a valid full mexa file but not the other way around. Full mexa adds two
different kinds of features: Abbreviations aka syntactic sugar as well as the concept
of *references*. This is one way how the above example could be written, exploiting
some of the `full-mexa` features:

```
name = "AdvectionTest"  # this is the name of the simulation

domain/size/width = 10
domain/size/height = 10

di/min = 1e-3
di/max = 0.5

grid/dx = di
grid/dy = di

e = "exact"
bc/left = e
bc/right = e
bc/top = e
bc/bottom = e

id/name = "BlackHole"
bhcode/pos = (0, 0)
bhcode/M = 1.0
bhcode/a = 0.9

initial_data = id
boundary_conditions = bc
```

Notice the appearance of non-quoted text on the right hand side of the equality.
This is one of the many features of the mexa language. For a full language
reference, see below.

# Functional Mexa Language Reference

This language reference seperates into the following sections:

  1. Basic language structure
  2. Four basis operations equality, extension, appending and inclusion
  3. Syntactic sugar for the right hand side
  

## Basic structure and evaluation

Every mexa file is also a simple-mexa file. Each line therefore reads as
`LHS OP RHS`, up to comments. We encourage to write this also as `OP(LHS,RHS)`
which resembles the internal representation notation in languages such as C++
and Python.

The left hand side (LHS) always defines a (hierarchic) *symbol*. The right
hand side (RHS) is subject to the interpreter. Mexa files are interpreted down to
simple-mexa files. This step typically contains three points:

  1. Parse the mexa file to some list of operations
  2. Evaluate the syntactic sugar down to operations, sort them in
  3. Evaluate the equalitites and extensions
  
Without doubt, due to the structure of the language, step 3 is most complex. In
fact, it turns out that a graph representation of the file is very helpful at
this step and graph traversing algorithms for resolving variables are suitable.

At this last step, the `mexa.py` sets up a graph with two kind of relations: The
directed extension `E(a,b)` and the path parentalship relationship `P(a,b)`. In
this language

  * a symbol `a/b/c` translates to the edge list `P(a,b),P(b,c)`,
  * the expression `a=b` translates to `E(a,b), E(b,a)`,
  * the expression `a<b` translates to `E(a,b)`
  
It is then subject to very carefully crafted algorithms to resolve this graph and
detect all valid equalitites by not missing any faulty ones.

## Operations (Relations)

Mexa understands to establish four different kind of relationships between graph
edges:

   1. The immutable equality `a = b`
   2. Datastructure extension (inheritance): `a <= b`  (abbreviated as `a < b`)
   3. The list creation (appending): `a += b`
   4. Inclusion of other files at any level: `a << b`

### Equality (assignment) for leafs

The assignment on two nodes is easy to understand:

```
a = 5
b = a
```

just means that `b=5`. Note however that in mexa, the order plays no role. This
means we have *lazy evaulation* of such expressions. Indeed, you might also write

```
b = a
a = 5
```

and get the same result, `a=5` and `b=5`. This feature gives the author of a mexa file
a lot of freedom to group relevant stuff.

Furthermore, assignments in mexa are immutable. Therefore, the assignment hasthe
mathematical role of an *equality*. `a` and `b` are *always* `5` and never change
their value.

### Equality for nodes

The example

```
foo/a = 5
foo/b = 6
bar/c = 7
bar = foo
```

reads *foo and bar are equal*. Subsequently, in the end both foo and bar have three
childs with keys `a, b, c` and values `5, 6, 7`. If both `foo` and `bar` would have
defined similarly named children, the file would not compile. This is the big
difference between equality and assignment as in an imperative language. Note that
again the order of rows plays no role at all.

Node equality also traverses down to subnodes and applies on all levels, i.e.

```
dog = cat
cat/location/lat = 57
cat/location/lng = 50
cat/name = "Peter"
dog/color = "black"
```

now both dog and cat have a location and a name.

### Extension

The extension operation is a directed version of the equality. In fact, the equality
`a = b` is evaluated as `a <= b` and `b <= a`. Extension vs equality makes no
difference for leafs. For nodes, however, it does:

```
triangle/x = 1
triangle/y = 2
triangle/z = 3
rectangle/w = 0
rectangle <= triangle
```

This does assign the rectangle the nodes `w,x,y,z` while the triangle is not affected.

With the extension operation, complex data structures can be built while no
**information is lost** in the middle. This is a non-destructive version of the typical
imperative approach of *default values* which are *overwritten* subsequently. The lack
of default values makes it obvious when an overwrite did not take place already at
compile time.

### List creation

The Mexa lang is very much about building up dictionary trees, but is also much about
building up lists. Therefore, it allows to write

```
lst += a
lst += b
lst += c
```

to mean

```
lst/a = a
lst/b = b
lst/c = c
```

In case the RHS symbol contains a path, only the node name is taken. In such cases
it can quickly become unclear what kind of naming happens, in case it is relevant
one should not use the list creation operator but instead regular assignments.
An example:

```
lst += foo/bar  # evaluates to lst/bar = foo/bar
lst += bla/baz/foo  # evaluates to lst/foo = bla/baz/foo
```

In case the RHS is not a symbol but a primitive, it is up to the interpreter to
come up with a reasonable key name. The only rule is that it must alphabetically
preserve the ordering. That is,

```
lst/x = 5
lst += 6    # may evaluate to lst/x1 = 6
lst += 7    # may evaluate to lst/x2 = 7
lst/y = 8
```

shall result in a list of upcounting numbers.

The list creation operator also simplifies the creation of list of strings. They
just read

```
foo += "bla"  # could evaluate to foo/s1 = "bla"
foo += "baz"  # could evaluate to foo/s2 = "baz"
foo += "bim"  # could evaluate to foo/s3 = "bim"
```

It is up to the interpretation of the mexa file to understand a list of strings as
a multiline string, thought.

### Inclusion operation

The `a << "something"` operation is somewhat special and has a similar meaning as in
Mathematica: It allows to place code from another file at exactly this node. On the
right hand side of this operator only strings are accepted (one could think of
accepting symbol names, but this destroys the independency of line evaluation order).

If `ext.txt` contains the lines `a=1` and `b=2`, then `foo << "ext.txt"` results in
`foo/a=1` and `foo/b=2`.

## Mexa high level syntactic sugar

For the right hand sides (RHS), mexa allows a number of convenient expressions which
allow to write more compact expressions.

### Vectors

Vectors can be easily noted with round or square brackets as
`a = (0,1,2,3,4)` or `a = [0,1,2,3,4]`. Curly brackets are reserved at the moment.
Such an expresison shall evaluate to an ordered dictionary with alphabetically
ordered keys. The actual naming of the key is up to the interpreter. It could be
for instance

```
a/a = 0
a/b = 1
a/c = 2
a/d = 3
a/e = 4
```

or also `i0,i1,i2,...` are common schemes. For vectors of length `3`, `x,y,z` is
a naming which could be useful.

Vectors can also be specified as RHS of append operations. The implication is
just the logical evaluation in the order of occurance, i.e. `a+=0`, then
`a+=(1,2)`, then `a+=(3,4)` is equal to our example above.

### Shorter Appending

As a special rule of vectors, without delimiters and brackets, on the list creation
operator multiple items may appear on the right hand side:

```
a += b c d
```

evaluates to `a += b`, `a += c`, `a += d` before the individual creation operations
are evaluated. This allows to save even more lines.

### String variable expansion

RHS variable expansion is well-defined for strings: In the string `"foo @bar baz"`,
the `@bar` is replaced by the variable `bar`. The lookup is done following the
same rules as a reference to the `bar` symbol at this very place. In case `bar`
does not resolve to a string, the input form as given in the file shall be used.

Symbol names can be taken into curly braces, such as `"foo@{bar}baz"`. This syntax
also allows names such as `@{bar/baz}` or `@{bar::baz}`, respectively.

Variable names in strings in single quotes shall not be evaluated.
