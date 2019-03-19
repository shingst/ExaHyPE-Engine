# helper script for Dominic; will be removed later
rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels  2> /dev/null
../../../Toolkit/toolkit.sh -d ../GRMHD_Limiting.exahype
ln -sf ../../../ExternalLibraries/TOVSolver/src/ tovsolver_lib
