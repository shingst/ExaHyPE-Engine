# The ExaHyPE meta specfile language (smexa)

This directory holds the exahype meta specfile language, also abbreveated
as *mexa*. The mexa ecosystem covers quite a broad range of ideas and
implementations.

## Rationale

Mexa was born due to the limitations of the parameter language shipped with
ExaHyPE. It is supposed as an on-top solution, to embed sophisticated
configurations into the mentioned parameter language. To do so, a syntax
was designed which ends up in a purely functional approach to parameter 
handling. However, mexa tries not to replace one restrictive configuration
language by another. Instead, the tooling tries to enhance interoperability,
openness and compatibility to other data representation formats. Users are
encouraged to choose the kind of file format and scripting language which
fits their needs best.

## The mexa language hierarchy

The most sophisticated part of the *mexa* ecosystem is the full-blown
immutable parameter language. It looks a bit like the Cactus parameter
language but allows variables and arbitrary ordering. Here is an
kitchen sink example:

```
# This is a comment.
# Arbitrarily structured parameters:
ExaHyPE::project::name = "GRMHD"
ExaHyPE::paths::peano_kernel_path = "./Peano"
ExaHyPE::architecture = "noarch"
ExaHyPE::logfile = "whatever.log"

Solver::type = "Finite-Volumes"
Solver::patch-size = 8

BasePlotter::time = 0.0
BasePlotter::repeat = 0.001
Paths::output_folder = "./vtk-output"

# inheritance:
IntegralPlotter <= BasePlotter
IntegralPlotter::output = "@{Paths/output_folder}/foo.txt"

# include files:
Solver::plotter_list << "solver-list.txt"

# lists
Domain::exclude  = "left"
Domain::exclude += "right"
Domain::exclude += "top"

# references
ExaHyPE::Solver = Solver
ExaHyPE::plotter_list = Solver::plotter_list

# hierarchies: "Parameters" gets three children
Parameters += ExaHyPE
Parameters += Domain
Parameters += IntegralPlotter
```

The language is specified at `doc/MEXALANG.md`. We refer to this as
**Level 1**. The simple version of this language is called `mexa-simple`.
It is just another (hierarchical) configuration language but very easy
to parse: 

```
initialdata::type = "blackhole-puncture"
initialdata::blackhole::mass = 1
initialdata::blackhole::isSpinning = True
initialdata::blackhole::spin = 0.9
```

The simple mexa language is again specified at `doc/MEXALANG.md`. We
refer to this as **Level 2**.

Contents of a simple mexa file can be represented in any other parameter
language such as Yaml, JSON or XML without any information loss. Especially
they can be embedded straight forwardly. Thus we approach the third and
most abstract level of this hierarchy: The hierarchical key-valued
data structure. It can be easily represented in any programming language.
We added this level to stress that parameter handling is a very basic form
of data managament which can *meet* on the lowest common denominator. We
refer to this as **Level 3**.

## Overview of this directory


**cpp/**: An implementation of a Level 3 code which has bindings to a
C++ Standard template Library (STL) representation, a simple-mexa file
parser as well as decodes BASE64/Quoted-printable simple-mexa files.

**mexa.py**: The swiss army knife reference implementation of the Level 1
mexa language. As thus, it can convert between the mexa level 1, level 2
and other file formats such as Yaml, XML and JSON. It also generates
ExaHyPE specification files. More documentation on this tool is given
below.

**examples**: Example mexa Level 1 files to test the interpreter. They
primarily serve as test cases.

**doc**: Documentation of the mexa paradigm, file formats and much more.
We encourage you to browse throught these files. They explain many
concepts in-depth.

## mexa.py Usage

The `mexa.py` reference implementation serves as a swiss army knife
to convert between the mexa file format and something else. We support
the creation of both *JSON* and *YAML* and both *structured* and *flat*
output files. Here, *structured* means that a nested dictioniary structure
is created while *flat* means the parameters are in one long list.

Thanks to JSON, one can for instance use the `jq` query tool to write
Queries to ask for certain values in the parameter file (*todo: examples*).

Of course, also ExaHyPE specification files can be generated. We either
propose to *serialize* (i.e. embed) the extensive parameters which can
not be stored in ExaHyPE in certain specfile fields or to let a trimmed
down Mexa file go along a specification file.

In any case, we probably want to write a simple C parser for Mexa parameters.

Mexa can also dump an *evaluated* version of the Mexa files where all
operations *except* the assignment where carried out, including string
variable substitution. Thus, the remaining Mexa file format can be parsed
trivially line by line without taking any logic into account.


## License

The mexa approach was invented by SvenK in Dec 2017 for ExaHyPE. All code
and texts is released under GPL/BSD. (c) 2017, 2018 http://exahype.eu
