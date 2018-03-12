#include "H5Writer.h"


GRMHD::H5Writer::H5Writer(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}


GRMHD::H5Writer::~H5Writer() {
  // @todo Please insert your code here
}


void GRMHD::H5Writer::startPlotting(double time) {
  // @todo Please insert your code here
}


void GRMHD::H5Writer::finishPlotting() {
  // @todo Please insert your code here
}


void GRMHD::H5Writer::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<25; i++){ 
    outputQuantities[i] = (i<19 ? Q[i] : 0.0);
  }
  /*
  // abuse some output quantities:
  int start = 14;
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


