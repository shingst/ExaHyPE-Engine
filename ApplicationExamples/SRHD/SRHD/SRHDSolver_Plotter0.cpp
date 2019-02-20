#include "SRHDSolver_Plotter0.h"


SRHD::SRHDSolver_Plotter0::SRHDSolver_Plotter0(SRHDSolver&  solver) {
  // @todo Please insert your code here
}


SRHD::SRHDSolver_Plotter0::~SRHDSolver_Plotter0() {
  // @todo Please insert your code here
}


void SRHD::SRHDSolver_Plotter0::startPlotting(double time) {
  // @todo Please insert your code here
}


void SRHD::SRHDSolver_Plotter0::finishPlotting() {
  // @todo Please insert your code here
}


void SRHD::SRHDSolver_Plotter0::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  for (int i=0; i<5; i++){ 
    // Conserved quantities, just copy
    outputQuantities[i] = Q[i];
  }
}


