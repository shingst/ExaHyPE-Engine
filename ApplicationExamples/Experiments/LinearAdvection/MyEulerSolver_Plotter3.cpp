#include "MyEulerSolver_Plotter3.h"
#include "InitialData.h"


Euler::MyEulerSolver_Plotter3::MyEulerSolver_Plotter3(MyEulerSolver& solver) {
  // @todo Please insert your code here
}


Euler::MyEulerSolver_Plotter3::~MyEulerSolver_Plotter3() {
  // @todo Please insert your code here
}


void Euler::MyEulerSolver_Plotter3::startPlotting(double time) {
  this->time = time;
}


void Euler::MyEulerSolver_Plotter3::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::MyEulerSolver_Plotter3::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  /**
   * This is the plotter for the exact solutions, given
   * as primitive Variables
   **/
  
  double xpos[DIMENSIONS];
  for(int i=0; i<DIMENSIONS; i++) xpos[i] = x[i];

  MovingGauss2D(xpos, outputQuantities, this->time);
  
  //cons2prim(outputQuantities, Q);

  /*
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
  */
}


