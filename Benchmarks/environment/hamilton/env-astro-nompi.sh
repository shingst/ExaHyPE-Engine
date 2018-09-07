module purge
module load python/3.6.3
module load java
module load slurm
module load intel/xe_2018.2
module load gcc/4.9.1
module load gsl/intel

# Not required
# export TBB_SHLIB="-L/ddn/apps/Cluster-Apps/intel/xe_2017.2/tbb/lib/intel64/gcc4.7 -ltbb"

export EXAHYPE_CC=icpc

export COMPILER_LFLAGS=" -lgsl -lgslcblas -lm "

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=None
export GPROF=off
