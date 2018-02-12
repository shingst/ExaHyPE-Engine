echo "Configure project for plenty-nodes scaling test (output)."
mkdir plenty-nodes/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/coolmuc3_KNL/Euler_ADERDG/plenty-nodes/Euler_ADERDG-output.exahype )
