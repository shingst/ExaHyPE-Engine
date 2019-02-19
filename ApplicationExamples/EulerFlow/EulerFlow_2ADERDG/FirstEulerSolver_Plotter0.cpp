#include "FirstEulerSolver_Plotter0.h"


Euler::FirstEulerSolver_Plotter0::FirstEulerSolver_Plotter0(FirstEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::FirstEulerSolver_Plotter0::~FirstEulerSolver_Plotter0() {
  // @todo Please insert your code here
}


void Euler::FirstEulerSolver_Plotter0::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::FirstEulerSolver_Plotter0::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::FirstEulerSolver_Plotter0::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
}


