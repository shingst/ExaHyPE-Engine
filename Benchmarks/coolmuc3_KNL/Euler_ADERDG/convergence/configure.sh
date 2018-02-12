echo "Configure project for convergence test."
mkdir convergence/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/coolmuc3_KNL/Euler_ADERDG/convergence/Euler_ADERDG.exahype )
