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

#include "EulerSolver_ADERDG.h"
#include "EulerSolver_ADERDG_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include <algorithm>

#include <string>

#include <math.h>

#include "peano/utils/Loop.h"

#include "kernels/KernelUtils.h"

#include "kernels/GaussLegendreQuadrature.h"

tarch::logging::Log Euler::EulerSolver_ADERDG::_log("Euler::EulerSolver_ADERDG");

void Euler::EulerSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}

void Euler::EulerSolver_ADERDG::flux(const double* const Q, double** const F) {
  constexpr double gamma = 1.4;
  const double irho      = 1./Q[0];
  const double j2        = Q[1]*Q[1]+Q[2]*Q[2];
  const double p         = (gamma-1) * (Q[3] - 0.5*irho*j2);
  
  // col 1
  F[0][0] = Q[1];
  F[0][1] = irho*Q[1]*Q[1] + p;
  F[0][2] = irho*Q[2]*Q[1];
  F[0][3] = irho*(Q[3]+p)*Q[1];
  
  // col 2
  F[1][0] = Q[2];
  F[1][1] = irho*Q[1]*Q[2];
  F[1][2] = irho*Q[2]*Q[2] + p;
  F[1][3] = irho*(Q[3]+p)*Q[2];
}

void Euler::EulerSolver_ADERDG::eigenvalues(const double* const Q,
    const int direction,
    double* const lambda) {
  constexpr double gamma = 1.4;
  const double irho      = 1./Q[0];
  
  const double u_n = Q[direction + 1] * irho;
  
  const double j2  = Q[1]*Q[1]+Q[2]*Q[2];
  const double p   = (gamma-1) * (Q[3] - 0.5*irho*j2);
  const double c   = std::sqrt(gamma * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n + c;
}

/**
 * (Smooth solution)
 *
 * Entropy wave is a moving Gaussian matter distribution where it is simple
 * to give an analytic result.
 *
 * See also chapter 7.13.2 in "I do like CFD, VOL.1" by Katate Masatsuka.
 */
void referenceSolution(const double* const x,const double t,double* const Q) {
  constexpr double gamma = 1.4;
  constexpr double p     = 1.0;
  constexpr double v0    = 0.5;
  constexpr double width = 0.3;

  const double distX    = x[0] - 0.5 - v0 * t;
  const double distY    = x[1] - 0.5;
  const double distance = std::sqrt(distX*distX + distY*distY); 
  
  Q[0] = 0.5 + 1.0 * std::exp(-distance / std::pow(width, DIMENSIONS));
  Q[1] = Q[0] * v0;
  Q[2] = 0;
  // total energy = internal energy + kinetic energy
  Q[3] = p / (gamma-1)  +  0.5*Q[0] * (v0*v0); 
}

void Euler::EulerSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt, double* const Q) {
  if (tarch::la::equals(t, 0.0)) {
    referenceSolution(x,t,Q);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::EulerSolver_ADERDG::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  double largestRho  = -std::numeric_limits<double>::infinity();
  for (int i=0; i<(Order+1)*(Order+1); i++) { // loop over all nodes
    const double* Q = luh + i*NumberOfVariables;
    largestRho  = std::max (largestRho,   Q[0]);
  }
  
  if ( largestRho/1.5 > 0.75 ) { // maximum is 1.5
    return RefinementControl::Refine;
  } else if ( largestRho/1.5 < 0.6 ) {
    return RefinementControl::Erase;
  } else {
    return RefinementControl::Keep;
  }
}

void Euler::EulerSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  referenceSolution(x,t+0.5*dt,stateOut);
  
  double _F[3][NumberOfVariables]={0.0};
  double* F[3] = {_F[0], _F[1], _F[2]};
  flux(stateOut,F);
  for (int i=0; i<NumberOfVariables; i++) {
    fluxOut[i] = F[direction][i];
  }
}
