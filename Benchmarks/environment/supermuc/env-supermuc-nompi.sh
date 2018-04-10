module switch python/3.5_intel
module switch java/1.8
module switch intel/18.0
module switch tbb/2018
module switch gcc/5

export EXAHYPE_CC=icpc

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=None
export ARCHITECTURE=hsw
export GPROF=off

# optimised kernels
export USE_IPO=on
