INCLUDES="-I../libs/easi/include"
INCLUDES=$INCLUDES" -I.."
INCLUDES=$INCLUDES" -I../libs/yaml-cpp/include"
INCLUDES=$INCLUDES" -I../libs/ImpalaJIT/include"
INCLUDES=$INCLUDES" -I../libs/ASAGI/include"
INCLUDES=$INCLUDES" ${MPI_INC}"


LIBRARIES=$LIBRARIES" -L../libs/yaml_build -lyaml-cpp"
LIBRARIES=$LIBRARIES" -L../libs/asagi_build -lasagi_nompi"
LIBRARIES=$LIBRARIES" -L ../libs/ImpalaJIT/lib -limpalajit"


set -x
module load gcc/5
module load mpi.intel
icpc -g -o test.prog --std=c++11 -DNOMPI -DUSE_ASAGI $INCLUDES easi_test.cpp $LIBRARIES $LIB 
