#include "MHDSolver_Plotter4.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"


MHDSolver::MHDSolver_Plotter4::MHDSolver_Plotter4(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::MHDSolver_Plotter4::~MHDSolver_Plotter4() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter4::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter4::finishPlotting() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter4::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const int nDim = MHDSolver::MHDSolver::nDim;
  const int nVar = MHDSolver::MHDSolver::nVar;
  
  // convert tarch::la::Vector to double*  
  double xpos[nDim];
  double exact[nVar];
  for(int i=0; i<nDim; i++) xpos[i] = x[i];

  alfenwave_(xpos, exact, &timeStamp);
  
  // compute the difference
  for (int i=0; i<nVar; i++){ 
    outputQuantities[i] = ( Q[i] - exact[i] );
  }
}


