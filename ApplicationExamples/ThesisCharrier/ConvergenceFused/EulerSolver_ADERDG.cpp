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


#include "kernels/GaussLegendreBasis.h"

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
 * We run the Ma=3.5 test described in [2]. Domain is [-1,1]. 
 *
 * [1] F. J. Hindenlang and G. J. Gassner, On the order reduction of entropy stable DGSEM for the compressible Euler equations, arXiv:1901.05812 [math], Jan. 2019.
 * [2] J. Chan, On discretely entropy conservative and entropy stable discontinuous Galerkin methods, Journal of Computational Physics, vol. 362, pp. 346-374, Jun. 2018.
 */
void Euler::EulerSolver_ADERDG::referenceSolution(const double* const x,const double t,double* const Q) {
  constexpr double gamma = 1.4;
  constexpr double p     = 1.0;
  constexpr double v[2]  = {2.5,2.4}; 
  constexpr double width = 0.3;

  Q[0] = 1 + 0.1 *std::sin(M_PI*((x[0]-v[0]*t) + (x[1]-v[1]*t))); 
  Q[1] = Q[0] * v[0];
  Q[2] = Q[0] * v[1];
  //// total energy = internal energy + kinetic energy
  Q[3] = p / (gamma-1)  +  0.5*Q[0] * (v[0]*v[0]+v[1]*v[1]); 
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

void Euler::EulerSolver_ADERDG::resetGlobalObservables(GlobalObservables& globalObservables) const {
  globalObservables.eL1()   = 0;
  globalObservables.eL2()   = 0;
  globalObservables.eLInf() = 0;
}

void Euler::EulerSolver_ADERDG::mapGlobalObservables(
    GlobalObservables&                          globalObservables,
    const double* const                         luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const double t,
    const double dt) const {
  constexpr int QuadOrder = 9;
  const auto cellOffset = cellCentre - 0.5 * cellSize;
  for (int iy=0; iy<QuadOrder+1; iy++) { // loop over all nodes
    for (int ix=0; ix<QuadOrder+1; ix++) { // loop over all nodes
      // quadrature values
      int    i[DIMENSIONS] = {ix,iy};
      double x[DIMENSIONS] = {0.0};
      double J_w = 1.0;
      for (int d=0; d<DIMENSIONS; d++) {
	x[d] = cellOffset[d] + cellSize[d] * kernels::legendre::nodes[QuadOrder][i[d]]; 
        J_w *= kernels::legendre::weights[QuadOrder][i[d]] * cellSize[d];
      } 
      // reference values
      double Qana[NumberOfVariables] = {0.0};
      referenceSolution(x,t,Qana);
      
      // solution values
      double Q[NumberOfVariables] = {0.0};
      for (int v=0; v < NumberOfVariables; v++) {
        Q[v] = kernels::legendre::interpolate(
          cellOffset.data(),
          cellSize.data(),
          x,
          NumberOfVariables,
          v,
          Order,
          luh
        );
      }
      
      // errors 
      double eL1   = 0;
      double eL2   = 0;
      double eLInf = 0;
      for (int v=0; v<NumberOfVariables; v++) {
        const double dQ = std::abs(Q[v]-Qana[v]);
        eL1  += dQ*J_w;
        eL2  += dQ*dQ*J_w;
        eLInf = std::max(eLInf,dQ);
      }

      globalObservables.eL1()   += eL1;
      globalObservables.eL2()   += eL2;
      globalObservables.eLInf() = std::max(globalObservables.eLInf(),eLInf);
    }
  }
}

void Euler::EulerSolver_ADERDG::mergeGlobalObservables(
    GlobalObservables&         globalObservables,
    ReadOnlyGlobalObservables& otherObservables) const {
  globalObservables.eL1()   += otherObservables.eL1();
  globalObservables.eL2()   += otherObservables.eL2();
  globalObservables.eLInf()  = std::max(globalObservables.eLInf(),otherObservables.eLInf());
}

void Euler::EulerSolver_ADERDG::wrapUpGlobalObservables(GlobalObservables& globalObservables) const {
  globalObservables.eL2() = std::sqrt(globalObservables.eL2());
}
  
void Euler::EulerSolver_ADERDG::beginTimeStep(const double minTimeStamp,const bool isFirstTimeStepOfBatchOrNoBatch) {
  ReadOnlyGlobalObservables observables = getGlobalObservables();
  logInfo("beginTimeStep(...)","sweep/reduce/time=" <<minTimeStamp);
  logInfo("beginTimeStep(...)","sweep/reduce/eL1="  <<observables.eL1());
  logInfo("beginTimeStep(...)","sweep/reduce/eL2="  <<observables.eL2());
  logInfo("beginTimeStep(...)","sweep/reduce/eLInf="<<observables.eLInf());
}
