#include "SRHDSolver_Plotter1.h"


SRHD::SRHDSolver_Plotter1::SRHDSolver_Plotter1(SRHDSolver&  solver) {
  // @todo Please insert your code here
}


SRHD::SRHDSolver_Plotter1::~SRHDSolver_Plotter1() {
  // @todo Please insert your code here
}


void SRHD::SRHDSolver_Plotter1::startPlotting(double time) {
  // @todo Please insert your code here
}


void SRHD::SRHDSolver_Plotter1::finishPlotting() {
  // @todo Please insert your code here
}


extern "C" {
// Fortran function, implemented in C2P-SRHD.f90
void pdecons2prim_(double* const V, double* const Q, int* iErr);
}

void SRHD::SRHDSolver_Plotter1::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  // Convert Q to Primitive quantities
  int iErr;
  pdecons2prim_(outputQuantities, Q, &iErr);
}


