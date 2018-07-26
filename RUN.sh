#!/bin/sh
#
# This script is a mini guide how to compile and run ExaHyPE.
# Also look into README.md as well as the latest guidebook
# which you can find at http://dev.exahype.eu/guidebook.pdf
#

# stop on error
set -e

# download the grid generation framework Peano.
./Peano/updatePeano || ./Peano/fetchSnapshot.sh

# compile the Java toolkit
bash ./Toolkit/build.sh

# Optional: Install libxsmm

# git clone https://github.com/hfp/libxsmm.git Libxsmm
# cd Libxsmm/
# make generator

# Now you are ready to follow compile and run an ExaHyPE application
# according to the guidebook:
java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ApplicationExamples/EulerFlow.exahype

# set build parameters
export CC=GNU
#export SHAREDMEM=TBB
#export TBB_INC=/usr/include/tbb
#export TBB_SHLIB="-L/usr/lib -ltbb"

# build sample application
cd ApplicationExamples/EulerFlow && make clean && make -j $(nproc)

# run sample application
./ExaHyPE-EulerFlow ../EulerFlow.exahype



