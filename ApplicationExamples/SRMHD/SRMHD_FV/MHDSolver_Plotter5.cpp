#include "MHDSolver_Plotter5.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"

#include <limits> // double max
#include <cmath> // C++11: std::signbit

MHDSolver::MHDSolver_Plotter5::MHDSolver_Plotter5(MHDSolver&  solver) {
  // @todo Please insert your code here
}


MHDSolver::MHDSolver_Plotter5::~MHDSolver_Plotter5() {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter5::startPlotting(double time) {
  // @todo Please insert your code here
}


void MHDSolver::MHDSolver_Plotter5::finishPlotting() {
  // @todo Please insert your code here
}

// There is also http://stackoverflow.com/a/25909704 not relying on C++11
template<typename T>
T sign(T arg) {
  // Returns the signum of arg. Works especially with NaNs and zeros.
  // signbit(arg) is true if arg is negative
  return std::signbit(arg) ? T(-1) : T(+1);
}

void MHDSolver::MHDSolver_Plotter5::mapQuantities(
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
  for (int i=0; i<nVar; i++) {
    outputQuantities[i] = ( Q[i] - exact[i] ) / exact[i];
    // Make NaN to finite value. Basically works like numpy.nan_to_num.
    if(outputQuantities[i] != outputQuantities[i] || outputQuantities[i] == std::numeric_limits<double>::infinity()) {
        outputQuantities[i] = sign(outputQuantities[i]) * std::numeric_limits<double>::max();
    }
  }
}


