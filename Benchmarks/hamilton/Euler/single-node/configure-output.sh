echo "Configure project for single-node scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler/single-node/Euler-output.exahype )
mkdir single-node/results
rm *.o
