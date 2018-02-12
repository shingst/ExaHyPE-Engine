echo "Configure project for multicore scaling test (no output)."
directory=a2.0

mkdir $directory/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk Fortran/*.o
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/CCZ4/$directory/CCZ4-no-output.exahype )
