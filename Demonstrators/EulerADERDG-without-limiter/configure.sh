echo "Configure project for multicore scaling test (no output)."
mkdir multicore/results
rm -rf *.o cfiles.mk ffiles.mk cipofiles.mk kernels
../../Toolkit2/toolkit.sh EulerADERDG-without-limiter.exahype
