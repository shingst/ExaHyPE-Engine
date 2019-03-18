module switch intel/18.0
module switch tbb/2018
module switch gcc/5
module load gsl/2.4

export EXAHYPE_CC="mpicc"
export EXAHYPE_FC="mpif90"
export COMPILER=Intel
export MODE=Release
export DISTRIBUTEDMEM=MPI
export SHAREDMEM=TBB
