echo "Configure project for plenty-nodes scaling test (no output)."
mkdir plenty-nodes/results
rm -r *.o cfiles.mk ffiles.mk cipofiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/coolmuc2/Elastic3D/plenty-nodes/Elastic3D-no-output.exahype )
