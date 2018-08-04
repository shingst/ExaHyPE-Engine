#!/usr/bin/env python3

"""
This is the exahype.toolkit package and holds the code of the "toolkit2", also called
"toolkit NG", the successor of the java-based ExaHyPE toolkit which itself is the
main glue code generator for ExaHyPE supposed to be runned by users.

You can run the frontend on the command line with something like

$ python3 -m exahype.toolkit.frontend
"""


# Python does not like me when it comes to packages.
# Maybe we can get rid of this...
import sys
from os.path import join, dirname, basename
sys.path.append(join(dirname(__file__),"..","toolkit"))

from helper import BadSpecificationFile
from directories import DirectoryAndPathChecker
from solver import SolverGenerator
from plotter import PlotterGenerator
from kernelcalls import KernelCallsGenerator

#__all__ = [ "helper", "directories", "solver", "plotter", "kernelcalls"  ]
