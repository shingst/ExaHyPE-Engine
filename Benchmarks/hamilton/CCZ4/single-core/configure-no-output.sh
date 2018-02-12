echo "Configure project for single-core optimisation studies. (no output)."
mkdir single-core/results
rm -rf *.o cfiles.mk ffiles.mk cipofiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/CCZ4/single-core/CCZ4-no-output.exahype )
