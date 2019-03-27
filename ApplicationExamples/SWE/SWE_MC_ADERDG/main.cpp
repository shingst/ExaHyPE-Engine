#include "exahype/main.h"
#include <iostream>


#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"


int main(int argc, char** argv){
    MPI_Init( &argc, &argv );
    std::cout << "Im in the main" << std::endl;
    exahype::main(argc,argv);
    exahype::main(argc,argv);
    MPI_Finalize();
}
