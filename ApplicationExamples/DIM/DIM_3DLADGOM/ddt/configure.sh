rm -rf *.o cipofiles.mk ffiles.mk cfiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ApplicationExamples/DIM/DIM_3DLADGOM/ddt/DIM-debug.exahype )
