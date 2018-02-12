echo "Configure project for multicore scaling test (no output)."
mkdir multicore/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk Fortran/*.o
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/coolmuc3_KNL/CCZ4/multicore/CCZ4-no-output.exahype )
