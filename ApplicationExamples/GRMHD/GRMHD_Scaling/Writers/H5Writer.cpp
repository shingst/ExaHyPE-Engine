#include "H5Writer.h"

/*
 * In case you wonder, H5Writer has nothing to do with GRMHD
 * but is a tool for testing plotting devices.
 */


GRMHD::H5Writer::H5Writer(GRMHDSolver_FV&  solver) {}
GRMHD::H5Writer::H5Writer(GRMHDSolver_ADERDG&  solver) {}
GRMHD::H5Writer::~H5Writer() {}
void GRMHD::H5Writer::startPlotting(double time) {}
void GRMHD::H5Writer::finishPlotting() {}

void GRMHD::H5Writer::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  int start = -1;
  outputQuantities[start+1] = x(0);
  outputQuantities[start+2] = x(1);
  outputQuantities[start+3] = DIMENSIONS==3 ? x(2) : -1;
  outputQuantities[start+4] = pos(0);
  outputQuantities[start+5] = pos(1);
  outputQuantities[start+6] = DIMENSIONS==3 ? pos(2) : -1;
  outputQuantities[start+7] = offsetOfPatch(0);
  outputQuantities[start+8] = offsetOfPatch(1);
  outputQuantities[start+9] = DIMENSIONS==3 ? offsetOfPatch(2) : -1;
}


