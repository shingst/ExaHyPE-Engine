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
#include <cstring> // memset

#include <string>

#include <math.h>

#include "peano/utils/Loop.h"

#include "kernels/KernelUtils.h"

#include "kernels/GaussLegendreBasis.h"

tarch::logging::Log Euler::EulerSolver_ADERDG::_log("Euler::EulerSolver_ADERDG");

void Euler::EulerSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}

void Euler::EulerSolver_ADERDG::flux(const double* const Q, double** const F) {
  
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];
  #if DIMENSIONS==3
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3];
  #else
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);

  // col 1
  F[0][0] = Q[1];
  F[0][1] = irho*Q[1]*Q[1] + p;
  F[0][2] = irho*Q[2]*Q[1];
  F[0][3] = irho*Q[3]*Q[1];
  F[0][4] = irho*(Q[4]+p)*Q[1];

  // col 2
  /*
  F[1][0] = Q[2];
  F[1][1] = irho*Q[1]*Q[2];
  F[1][2] = irho*Q[2]*Q[2] + p;
  F[1][3] = irho*Q[3]*Q[2];
  F[1][4] = irho*(Q[4]+p)*Q[2];
  */
  
  // done in ncp
  F[1][0] = 0.;//Q[2];
  F[1][1] = 0.;//irho*Q[1]*Q[2];
  F[1][2] = 0.;//irho*Q[2]*Q[2] + p;
  F[1][3] = 0.;//irho*Q[3]*Q[2];
  F[1][4] = 0.;//irho*Q[4]*Q[2]+irho*p*Q[2];

  #if DIMENSIONS==3
  // col 3
  F[2][0] = Q[3];
  F[2][1] = irho*Q[1]*Q[3];
  F[2][2] = irho*Q[2]*Q[3];
  F[2][3] = irho*Q[3]*Q[3] + p;
  F[2][4] = irho*(Q[4]+p)*Q[3];
  #endif
}

void Euler::EulerSolver_ADERDG::nonConservativeProduct(const double* const Q, const double* const gradQ, double* const BgradQ) {
  
  //ncp is div(F(Q))
  
  std::memset(BgradQ,0,NumberOfVariables*sizeof(double));
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];

  // x
  
  // done in flux
  
  // y
  const double* const Q_dy = gradQ+NumberOfVariables;
  const double irho_dy = -Q_dy[0]*irho*irho;
  #if DIMENSIONS==3
  const double j2      = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
  const double j2_dy   = 2*Q[1]*Q_dy[1] + 2*Q[2]*Q_dy[2] + 2*Q[3]*Q_dy[3];
  #else
  const double j2      = Q[1]*Q[1] + Q[2]*Q[2];
  const double j2_dy   = 2*Q[1]*Q_dy[1] + 2*Q[2]*Q_dy[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);
  const double p_dy = (gamma-1) * (Q_dy[4] - 0.5*(j2_dy*irho+j2*irho_dy));
  
  BgradQ[0] += Q_dy[2];
  BgradQ[1] += irho_dy*Q[1]*Q[2] + irho*Q_dy[1]*Q[2] + irho*Q[1]*Q_dy[2];
  BgradQ[2] += irho_dy*Q[2]*Q[2] + irho*Q_dy[2]*Q[2]*2                  + p_dy;
  BgradQ[3] += irho_dy*Q[3]*Q[2] + irho*Q_dy[3]*Q[2] + irho*Q[3]*Q_dy[2];
  BgradQ[4] += irho_dy*Q[4]*Q[2] + irho*Q_dy[4]*Q[2] + irho*Q[4]*Q_dy[2] + irho_dy*p*Q[2] + irho*p_dy*Q[2] + irho*p*Q_dy[2];

  #if DIMENSIONS==3
  // z
  
  // done in flux
  #endif

}

void Euler::EulerSolver_ADERDG::eigenvalues(const double* const Q,
    const int direction,
    double* const lambda) {
  
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];
  #if DIMENSIONS==3
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3];
  #else
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);

  const double u_n = Q[direction + 1] * irho;
  const double c   = std::sqrt(gamma * p * irho);

  std::fill_n(lambda,5,u_n);
  lambda[0] -= c;
  lambda[4] += c;
}

void Euler::EulerSolver_ADERDG::entropyWave(const double* const x,double t, double* const Q) {
  const double gamma         = 1.4;
  constexpr double width     = 0.3;
  constexpr double amplitude = 0.3;

  #if DIMENSIONS==2
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1]);
  tarch::la::Vector<DIMENSIONS,double> v0(0.5,0.5);
  tarch::la::Vector<DIMENSIONS,double> x0(0.5,0.5);
  #else
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1],x[2]);
  tarch::la::Vector<DIMENSIONS,double> v0(0.5,0.5,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(0.5,0.5,0.5);
  #endif
  const double distance  = tarch::la::norm2( xVec - x0 - v0 * t );

  Q[0] = 0.5 + amplitude * std::exp(-distance / std::pow(width, DIMENSIONS));
  Q[1] = Q[0] * v0[0];
  Q[2] = Q[0] * v0[1];
  Q[3] = 0.0;
  // total energy = internal energy + kinetic energy
  const double p = 1.;
  Q[4] = p / (gamma-1)  +  0.5*Q[0] * (v0[0]*v0[0]+v0[1]*v0[1]); // v*v; assumes: v0[2]=0
}


void Euler::EulerSolver_ADERDG::referenceSolution(const double* const x,double t, double* const Q) {
  entropyWave(x,t,Q);
}

void Euler::EulerSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt, double* const Q) {
  if (tarch::la::equals(t, 0.0)) {
    referenceSolution(x,0.0,Q);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::EulerSolver_ADERDG::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  double maxE = -std::numeric_limits<double>::max();
  double minE = +std::numeric_limits<double>::max();
  
  const int nodes = std::pow((Order+1), DIMENSIONS);
  for(int i=0;i<nodes;i++) {
    maxE = std::max(maxE, luh[i*NumberOfVariables+NumberOfVariables-1]);
    minE = std::min(minE, luh[i*NumberOfVariables+NumberOfVariables-1]);
  }

  if ( maxE/minE > 1.05 ) {
    return RefinementControl::Refine;
  } else if ( level > getCoarsestMeshLevel() ) {
    return RefinementControl::Erase;
  } else {
    return RefinementControl::Keep;
  }
}

void Euler::EulerSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  double Q[NumberOfVariables]     = {0.0};
  double _F[3][NumberOfVariables] = {0.0};
  double* F[3] = {_F[0],_F[1],_F[2]};

  // initialise
  std::fill_n(stateOut, NumberOfVariables, 0.0);
  std::fill_n(fluxOut,  NumberOfVariables, 0.0);
  for (int i=0; i<Order+1; i++) {
    const double ti = t + dt * kernels::legendre::nodes[Order][i];
    referenceSolution(x,ti,Q);
    flux(Q,F);
    for (int v=0; v<NumberOfVariables; v++) {
      stateOut[v] += Q[v]            * kernels::legendre::weights[Order][i];
      fluxOut[v]  += F[direction][v] * kernels::legendre::weights[Order][i];
    }
  }
}

void Euler::EulerSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
    std::copy_n(Q, NumberOfVariables, observables);
}


bool Euler::EulerSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  if (observablesMin[0] <= 0.0) return false;
  if (observablesMin[4] < 0.0) return false;
  return true;
}
