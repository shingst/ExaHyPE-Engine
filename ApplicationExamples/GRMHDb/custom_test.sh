#module load lrztools

export COMPILER=Intel
export MODE=Release
export SHAREDMEM=TBB
export DISTRIBUTEDMEM=None
export OMP_NUM_THREADS=72              # numero di threads per processo
#export USE_IPO=off
# this is for IBM MPI
#export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS   # numero di cores assegnati ai threads
#if [ "$OMP_NUM_THREADS" -gt 1 ] ; then
#  module load mpi_pinning/hybrid_blocked
#else
#  module load mpi_pinning/mpp
#fi
export I_MPI_PIN=on
export I_MPI_PIN_DOMAIN=omp:compact     #   i.e. =OMP_NUM_THREADS
