# The Mexa configuration attemp concepts

The following concepts can be read top-bottom, i.e. they successively
develop the idea of *mexa*. These concepts shall be sorted into the
threefold of what mexa can be:

1. A schemaless hierarchical key-value storage
2. A very simple file format for parameters based on (1)
3. A purely functional language based on (2) to create parameter trees

### Associative key: Key-Value mapping

Configuration parameters means a certain configuration *key* is
assigned a *value*. Implementations of hash maps exist in all programming
languages. In Mexa, the key-value mapping is always scalar, i.e.
the types of the values are *primitives* such as Integers,
floating point numbers, booleans or character strings.

### Hierarchical directory tree

Many configuration systems are structured in a hierarchy, think of 
the ini format, LDAP or the Microsoft Windows Registry. In such a
hierarchy, there is a clear distinction between *folders* and *files*,
in the language of graph theory *nodes* and *leafs*.

Some programming languages, notably untyped script languages, allow
to represent hierarchical (nested) directories (hash maps) directly
in code, most notably Python and PHP.

### Syntax-free

So far, the concepts of *key-value-pairs* and **hierarchy* are pretty
common. They are not only quickly representable in most computing languages
but also in most configuration file formats, markup languages or data
serialization/exchange languages such as XML, Yaml and JSON. Therefore,
we want to actually *remove* the syntax from the idea of *mexa*. It is
thus *syntax-free*, the contents shall be representable in any way.

*Note:* This idea goes along the way how [Peano](http://www.peano-framework.org)
treats file formats: There is the *syntax*, for instance ASCII, HDF5
or a home-brewn IEE-754-embedding-file format, and there is the
*semantics* which could be thought of one layer *above* the physical
file format.

Actually, the *syntax-free* keyword is not generally aplicable.

### Verbosity and Extensivity

The `plain-mexa` language is verbose: The full path of the key is always
given in every assignment line. This is by intention and contradicts
principles such as *don't repeat yourself*.

On the other hand, `full-mexa` is an extensive language which allows you
to avoid *magic numbers* by any means. Instead, you are encouraged to
assign every possible primitive data (numbers and trees) a name and use
these names extensively.

### Convention/Contract vs. Schema

A schema (think of DTD/XSD or Cactus parameter.ccl files) allows to
*statically* check the semantical correctness of a parameter file. For
instance, it can enforce a given quantity to

- be in a number interval
- match a regexp
- be a valid file name
- have certain subentries following rules, each
- or have some default values

It can be furthermore used to document the parameter structure. Without
doubt, a schema is helpful.

Mexa was constructed without such a schema in mind. Instead, parameters
are given by convention and should be *complete*, which follows the
*principle of least surprise*: Every possible parameter is expressed
in the parameter file.

*Note*: This principle can hardly be followed in the practice as it allows to
change all existing parameter files when a new parameter is introduced
in the code.

### General-purpose

Needless to say, while mexa fulfills the idea of a domain specific language
(DSL), i.e. collecting steering parameters -- "a bunch of numbers" -- for
scientific simulations, it has nevertheless the claim to be a general-purpose
configuration language.

This claim fails pretty quickly when it comes to text representation: While
there has been a focus on scientific number and vector syntax, text
representation is treated poorly. It is possible to embed arbitrary multiline
strings (i.e. text files) in mexa files. The paradigm to do so is to
embed it naturally into a dumb dictionary which represents a list of lines.

### Embeddable

Mexa is primarily designed for being embedded in other configuration languages.
Therefore, standards have been defined in which way a mexa file can be encoded
and readers should support it.