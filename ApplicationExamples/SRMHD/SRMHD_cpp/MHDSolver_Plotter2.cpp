#include "MHDSolver_Plotter2.h"

#include "MHDSolver.h"

MHDSolver::MHDSolver_Plotter2::MHDSolver_Plotter2(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::MHDSolver_Plotter2::~MHDSolver_Plotter2() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter2::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter2::finishPlotting() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter2::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  // MHD primitive variable plotter
  int i;
  MHDSolver::Cons2Prim(Q, outputQuantities);
}


