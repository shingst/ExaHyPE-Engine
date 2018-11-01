#include "PrimitiveWriter.h"
#include "Fortran/PDE.h"


GRMHD::PrimitiveWriter::PrimitiveWriter(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}



GRMHD::PrimitiveWriter::PrimitiveWriter(GRMHDSolver_ADERDG&  solver) {
  // @todo Please insert your code here
}


GRMHD::PrimitiveWriter::PrimitiveWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

GRMHD::PrimitiveWriter::~PrimitiveWriter() {
  // @todo Please insert your code here
}


void GRMHD::PrimitiveWriter::startPlotting(double time) {
  numberOfC2PFailures = 0;
  allConversions = 0;
}


void GRMHD::PrimitiveWriter::finishPlotting() {
  if(numberOfC2PFailures > 0) {
    static tarch::logging::Log _log("GRMHD::PrimitiveWriter");
    logInfo("finishPlotting", "Counted " << numberOfC2PFailures << " Cons2Prim failures during plotting (" << (allConversions/numberOfC2PFailures*100) << "%)");
  }
  startPlotting(0);
}


void GRMHD::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
    int err;
    pdecons2prim_(outputQuantities, Q, &err);
    if(err != 0) {
	    //printf("Cons2Prim Failure in PrimitiveWriter!!!");
	    numberOfC2PFailures++;

	    constexpr int c2pFailureIndicatingMagicNumberForDensity = -1;
	    Q[0] = c2pFailureIndicatingMagicNumberForDensity;
    }
    allConversions++;
}


