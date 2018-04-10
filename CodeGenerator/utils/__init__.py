import os
import sys

# import current dir into sys.path for correct import of dependencies with symlinks
currDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currDir)
