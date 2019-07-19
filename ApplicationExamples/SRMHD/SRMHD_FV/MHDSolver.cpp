#include "MHDSolver.h"
//#include "fortran.h" _ltob

#include "InitialDataAdapter.h"
#include "PDE.h"

#include <memory>
#include <cstring>
#include <stdio.h>
#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h" // matrix indexing

/* This is the MHDSolver.cpp binding to Fortran functions, as done in SRHD. */


//void MHDSolver::MHDSolver::init(std::vector<std::string>& cmdargs, exahype::parser::ParserView& _constants) {
//  // just pass the pointer to the crazy Fortran glue code. Should be improved.
//  constants = &_constants;
//}

void MHDSolver::MHDSolver::flux(const double* const Q, double** const F) {
  // Caveats: Fortran accepts a uniform array of size (nVar*nDim), however C passes an array of pointers.
  // This Fortran interface works only if F is a continous array and F[1]==F[nDim+1] etc!
  pdeflux_(F[0], Q);
}



void MHDSolver::MHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* const lambda) {
  double nv[3] = {0.};
  nv[normalNonZeroIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}



bool MHDSolver::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return (t < 1e-10);
}



void MHDSolver::MHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  // Fortran call:
  adjustedsolutionvalues_(x, &w, &t, &dt, Q);
}

void MHDSolver::MHDSolver::algebraicSource(const double* const Q, double* const S) {
  pdesource_(S, Q);
}



exahype::solvers::Solver::RefinementControl MHDSolver::MHDSolver::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void MHDSolver::MHDSolver::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double* const stateIn, double* const stateOut) {

//  std::cout << "boundary=" << faceIndex <<  ", t=" << t <<  ", x2D={"<< x[0] << "," << x[1] << "}" << std::endl;
    
  // faceIndex: LEFT=0, RIGHT=1, BOT=2, TOP=3
  if(faceIndex == 0 || faceIndex == 1) {

    // FV Goudonov scheme exact BC: No integration neccessary, no fluxes neccessary!
    alfenwave_(x, stateOut, &t);
    return;
  } else {
    // outflow bc
    for (int i=0; i<nVar; ++i) {
       stateOut[i] = stateIn[i];
    }
    return;
  }

//  const double* normalVelocity = stateIn+1+normalNonZero;
//  const double sign = (faceIndex-2*normalNonZero)==0 ? -1.0 : 1.0;
//  const double outwardDirectedVelocity = sign * (*normalVelocity);
//
//  if (outwardDirectedVelocity>0) { // outflow; take inside value
//    // do nothing
//  } else { // inflow; take outside values
//    stateOut[0] = 1.0; // density for InitialBlast
//    stateOut[4] = 0.1; // pressure for InitialBlast
//  }

  // debugging stuff:
  //printf("stateOut[%d]=%.5e == stateIn[%d]=%.5e at t=%f, dt=%f x=%.5e, y=%.5e, Equals=%e\n", statem, stateOut[statem], statem, stateIn[statem], t, dt, x[0], x[1], stateOut[statem] - stateIn[statem]);
  //printf("FOut[%d]=%e == Fin[%d]=%e\n", statem, fluxOut[statem], statem, fluxIn[statem]);
}


//void MHDSolver::MHDSolver::nonConservativeProduct(const double* const Q, const double* const gradQ, double* const BgradQ) {
//	std::memset(BgradQ, 0, nVar * sizeof(double));
//}
//
//void MHDSolver::MHDSolver::coefficientMatrix(const double* const Q, const int normalNonZero, double* const Bn) {
// std::memset(Bn, 0, nVar * nVar * sizeof(double));
//}


