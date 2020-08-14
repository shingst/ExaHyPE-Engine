/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "MyEulerSolver.h"
#include "InitialData.h"

#include "MyEulerSolver_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include "peano/utils/Loop.h"

#include "kernels/KernelUtils.h"

#include <memory>
#include <cstring> // memset
#include <string>
#include <stdlib.h> // getenv

// global Parameters for Euler::MyEulerSolver.


void Euler::MyEulerSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // This guy usually should go into the header but it does not make sense to
  // commit a whole header only because of one file.
  static tarch::logging::Log _log("MyEulerSolver::init");

  // This function is called inside the generated constructor.
  // @todo Please implement/augment if required

  // Demonstration how to access parameters:
  logInfo("init(...)", "EulerFlow was called with the following parameters:");
  for(size_t i=0; i<cmdlineargs.size(); i++) {
    logInfo("init(...)", "- " << cmdlineargs[i].c_str() );
  }

  // Until Parameter access works, we use ENVIRONMENT variables
  const char* _id = std::getenv("EXAHYPE_INITIALDATA");
  const char* _bc = std::getenv("EXAHYPE_BOUNDC");

  std::string id("MovingGauss"), bc("outflow");

  if(_id) id=_id; else logInfo("ID", "Loading default Initial Data");
  if(_bc) id=_bc; else logInfo("BC", "Loading default Boundary Conditions");

  logInfo("ID", std::string("Loading Initial data: '")+id+std::string("'"));
  if(id == "ShuVortex") idfunc = ShuVortex2D;
  if(id == "MovingGauss") idfunc = MovingGauss;
  if(id == "DiffusingGauss") idfunc = DiffusingGauss;
  if(!idfunc) {
      logError("ID", "Cannot understand requested ID.");
      exit(-1);
  }

  logInfo("BC", std::string("Applying Boundary Conditions: '")+id+std::string("'"));
  if(id == "outflow") {
  }

//   idfunc = MovingGauss; //  TODO(Dominic): Use this function for the dyn amr test.
}

void Euler::MyEulerSolver::flux(const double* const Q, double** const F) {
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  f.rho ( vars.j()                                 );
  f.j   ( irho * outerDot(vars.j(),vars.j()) + p*I );
  f.E   ( irho * (vars.E() + p) * vars.j()         );
}

void Euler::MyEulerSolver::eigenvalues(const double* const Q,
                                       const int normalNonZeroIndex,
                                       double* const lambda) {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double GAMMA = 1.4;
  const double irho = 1./vars.rho();
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho()=u_n - c;
  eigs.E()  =u_n + c;
  eigs.j(u_n,u_n,u_n);
}


void Euler::MyEulerSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    // pass the time for exact initial data as t is not exactly 0.
//    idfunc(x, Q, 0.0);
    DiffusingGauss(x,Q,t);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::MyEulerSolver::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
//  if (t>1e-6)
//    return exahype::solvers::Solver::RefinementControl::Keep;


//  double largestRho  = -std::numeric_limits<double>::max();
//  double smallestRho =  std::numeric_limits<double>::max();
//
//  kernels::idx3 idx_luh(Order+1,Order+1,NumberOfVariables);
//  dfor(i,Order+1) {
//    ReadOnlyVariables vars(luh + idx_luh(i(1),i(0),0));
//
//    largestRho  = std::max (largestRho,  vars.rho());
//    smallestRho = std::min (smallestRho, vars.rho());
//  }
//
////  if ((largestRho-smallestRho)/largestRho > 0.9e-1) { // This doesn't work with the erasing
//  if (largestRho > 0.65) {
//    return exahype::solvers::Solver::RefinementControl::Refine;
//  }
//
//  if (dx[0] < getMaximumMeshSize()/3.)
//    return exahype::solvers::Solver::RefinementControl::Erase;

  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Euler::MyEulerSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

    // Compute boundary state.
//  DiffusingGauss(x, stateOut, t);
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];

    // Compute flux and
    // extract normal flux in a lazy fashion.
    double fi[DIMENSIONS][NumberOfVariables], *F[DIMENSIONS];
    for(int d=0; d<DIMENSIONS; d++) F[d] = fi[d];
    // it could also be done "more effective" with something like
    // fi[normalNonZero] = reinterpret_cast<double[5]>(fluxOut);
    F[normalNonZero] = fluxOut; // This replaces the double pointer at pos normalNonZero by fluxOut.
    flux(stateOut, F);

    // This copy is not neccessary as we do have one component of
    // F already pointing to fluxOut.
    /*
    for (int i=0; i<5; i++) {
      fluxOut[i] = F[normalNonZero][i];
    }
    */

  // The problem with these definitions is that in a simulation
  // with a global nonzero velocity (as in MovingGauss2D), errnous
  // values move into the simulation domain very quickly. So these
  // boundary conditions are not good at all. Instead, we should
  // have per default "vacuum" boundary conditions, so that vacuum
  // values enter the grid as soon as matter moves away.

  /*
  //  fluxOut
  //  //@todo Please implement
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  //  // stateOut
  //  // @todo Please implement
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
  */
}
