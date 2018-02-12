echo "Configure project for single-core optimisation studies. (no output)."
mkdir single-core/results
rm -rf *.o cfiles.mk ffiles.mk cipofiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Elastic3D/single-core/Elastic3D-no-output.exahype )
