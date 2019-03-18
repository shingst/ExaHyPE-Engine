#include "ConservedWriter.h"


GRMHD::ConservedWriter::ConservedWriter(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}

GRMHD::ConservedWriter::ConservedWriter(GRMHDSolver_ADERDG&  solver) {
  // @todo Please insert your code here
}

GRMHD::ConservedWriter::ConservedWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
}


GRMHD::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void GRMHD::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void GRMHD::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void GRMHD::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<1; i++){ 
    outputQuantities[i] = Q[i];
  }

  // as always, for debugging: 
  /*

  const int start=0;
  outputQuantities[start+1] = x(0);
  outputQuantities[start+2] = x(1);
  outputQuantities[start+3] = DIMENSIONS==3 ? x(2) : -1;
  outputQuantities[start+4] = pos(0);
  outputQuantities[start+5] = pos(1);
  outputQuantities[start+6] = DIMENSIONS==3 ? pos(2) : -1;
  outputQuantities[start+7] = offsetOfPatch(0);
  outputQuantities[start+8] = offsetOfPatch(1);
  outputQuantities[start+9] = DIMENSIONS==3 ? offsetOfPatch(2) : -1;
  */
}


