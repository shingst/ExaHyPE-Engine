#include "MHDSolver_Plotter3.h"

#include "MHDSolver.h"

MHDSolver::MHDSolver_Plotter3::MHDSolver_Plotter3(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::MHDSolver_Plotter3::~MHDSolver_Plotter3() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter3::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter3::finishPlotting() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter3::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  // convert tarch::la::Vector to double*  
  double xpos[DIMENSIONS];
  double V[10];
  for(int i=0; i<DIMENSIONS; i++) xpos[i] = x[i];

  // Caveat: I don't properly use initialdatabyexahypespecfile here as it doesn't pass
  // the time. However, the Alfen Wave is the only exact solution we currently have, anyway
  MHDSolver::AlfvenWave(xpos, V, timeStamp);
  MHDSolver::Cons2Prim(outputQuantities, V);
}


