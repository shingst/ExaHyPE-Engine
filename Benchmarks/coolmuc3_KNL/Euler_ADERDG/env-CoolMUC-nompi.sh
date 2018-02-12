source /etc/profile.d/modules.sh

module unload gcc
module switch python/3.5_intel
module switch java/1.8
module load gcc/4.9
module switch intel/17.0 
module switch tbb/2017

export EXAHYPE_CC=icpc

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=None
export GPROF=off
export USE_IPO=on
