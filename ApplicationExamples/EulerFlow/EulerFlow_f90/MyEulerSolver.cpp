#include "MyEulerSolver.h"

#include <memory>

void Euler::MyEulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // This function is called by the constructor.
  // You can access spec file parameters as well as command line arguments (argv as std::vector).
  // @todo Please implement/augment if required.
}

//************************************************* 
//for FORTRAN kernels the fluxes and eigenvalues 
//have to be implemented in the file ./PDE.f90 and ./typesDef 
//************************************************* 






void Euler::MyEulerSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 3
  // Number of variables    = 5 (#unknowns + #parameters)


  // @todo Please implement
  // fluxOut
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  // stateOut
  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
}






bool Euler::MyEulerSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  // @todo Please implement
  if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    return true;
  }
  return false;
}



void Euler::MyEulerSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement
  if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    const double GAMMA = 1.4;
    // @todo Please implement
    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    Q[4] =
        1. / (GAMMA - 1.) +
        std::exp(-sqrt((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) +
                   (x[2] - 0.5) * (x[2] - 0.5)) /((0.25)*(0.25)));
    //  Q[4] = 2.5;
  }
}



exahype::solvers::Solver::RefinementControl Euler::MyEulerSolver::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}









