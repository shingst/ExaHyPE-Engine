"""
exahype.specfiles is a Python3 package for dealing with ExaHyPE specification
files. It shall provide functions to Read/Convert and Write/Export files
from and to the following formats:

 * ExaHyPE Specification v1 files
 * ExaHyPE Specification v2 files (aka JSON files)
 * INI files? If you want...
 * YAML, if you want...

Parts of the codes may be just copied from https://bitbucket.org/svek/mexa
"""

# Python does not like me when it comes to packages.
# Maybe we can get rid of this...
import sys
from os.path import join, dirname, basename
sys.path.append(join(dirname(__file__),"..","specfiles"))

from validate import validate
