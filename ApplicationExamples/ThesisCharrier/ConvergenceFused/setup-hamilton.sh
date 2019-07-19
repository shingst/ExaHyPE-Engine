module load slurm
module load intelmpi/intel/2018.2
module load intel/xe_2018.2
module load gcc/7.3.0
module load python/3.6.8
module unload gcc/8.2.0

export COMPILER=Gnu
export MODE=Release
export DISTRIBUTEDMEM=None
export SHAREDMEM=TBB

../../../Toolkit/toolkit.sh -sd Euler.exahype2
