#!/bin/sh
#
# This script is a mini guide how to compile and run ExaHyPE.
# Also look into README.md as well as the latest guidebook
# which you can find at http://dev.exahype.eu/guidebook.pdf
#

# stop on error
set -e

# download the grid generation framework Peano
# and a number of python dependencies
# (use -h for information on different ways to clone the repositories.)
bash Submodules/updateSubmodules.sh

# Now you are ready to follow compile and run an ExaHyPE application
# according to the guidebook:
Toolkit/toolkit.sh Demonstrators/EulerFV/EulerFV.exahype

# set build parameters
export COMPILER=Gnu
export SHAREDMEM=None
export DISTRIBUTEDMEM=None
#export TBB_INC=/usr/include/tbb
#export TBB_SHLIB="-L/usr/lib -ltbb"

# build sample application (use -j for more make threads)
( cd Demonstrators/EulerFV && make clean && make -j24 )

# run sample application
( cd Demonstrators/EulerFV && ./ExaHyPE-EulerFV EulerFV.exahype )



