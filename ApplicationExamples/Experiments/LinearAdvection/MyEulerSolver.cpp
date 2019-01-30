
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

#include <memory>

void Euler::MyEulerSolver::init() {
  // This function is called inside the generated constructor.
  // @todo Please implement/augment if required
}

void Euler::MyEulerSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)

  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
#ifdef Dim2
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
#elif Dim3
  const double p =
      (GAMMA - 1) *
      (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
#else
#error Dim2 or Dim3 must be defined
#endif

  double* f = F[0];
  double* g = F[1];

  // @todo Please implement
  // f
  f[0] = Q[0]*Q[1]; // TODO Attention: I multiply here with Q[0] in contrast  to Euler.
  f[1] = 0; // irho * Q[1] * Q[1] + p;
  f[2] = 0; // irho * Q[1] * Q[2];
  f[3] = 0; // irho * Q[1] * Q[3];
  f[4] = 0; //irho * Q[1] * (Q[4] + p);
  // g
  // @todo Please implement
  g[0] = Q[0]*Q[2]; // TODO Attention: I multiply here with Q[0] in contrast  to Euler.
  g[1] = 0; // irho * Q[2] * Q[1];
  g[2] = 0; //irho * Q[2] * Q[2] + p;
  g[3] = 0; //irho * Q[2] * Q[3];
  g[4] = 0;  //irho * Q[2] * (Q[4] + p);

#ifdef Dim3
  double* h = F[2];
  // h
  // @todo Please implement
  h[0] = Q[0]*Q[3]; // TODO Attention: I multiply here with Q[0] in contrast  to Euler.
  h[1] = 0; //irho * Q[3] * Q[1];
  h[2] = 0; //irho * Q[3] * Q[2];
  h[3] = 0; //irho * Q[3] * Q[3] + p;
  h[4] = 0; //irho * Q[3] * (Q[4] + p);
#endif
}

void Euler::MyEulerSolver::algebraicSource(const double* const Q, double* S) {
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}

void Euler::MyEulerSolver::eigenvalues(const double* const Q,
                                       const int normalNonZeroIndex,
                                       double* lambda) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];

#ifdef Dim2
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
#elif Dim3
  double p = (GAMMA - 1) *
             (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
#else
#error Dim2 or Dim3 must be defined
#endif

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = 0; // std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = 0; //u_n;
  lambda[2] = 0; //u_n;
  lambda[3] = 0; //u_n;
  lambda[4] = 0; //u_n + c;
}

bool Euler::MyEulerSolver::hasToAdjustSolution(
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}

void Euler::MyEulerSolver::adjustedSolutionValues(const double* const x,
                                                  const double w,
                                                  const double t,
                                                  const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    // pass the time for exact initial data as t is not exactly 0.
    InitialData(x, Q, t);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::MyEulerSolver::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  // @todo Please implement
  if (dx[0] > getMaximumMeshSize()/9.) {
    if (center[1] < 1.0 && center[1] > 0.2) {
      if (center[0] < 0.66 && center[0] > 0.0) {
        return exahype::solvers::Solver::RefinementControl::Keep;
      }
    }
  }

  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Euler::MyEulerSolver::boundaryValues(const double* const x, const double t,
                                          const int faceIndex,
                                          const int normalNonZero,
                                          const double* const fluxIn,
                                          const double* const stateIn,
                                          double* const fluxOut, double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  //  InitialData(x,stateOut,t);

  // skip fluxes for the time being as it crashes

  // Compute boundary state.
  InitialData(x, stateOut, t);

  // Compute flux and
  // extract normal flux in a lazy fashion.
  double f[5];
  double g[5];
  #if DIMENSIONS == 2
  double* F[2];
  F[0] = f;
  F[1] = g;
  #else
  double h[5];
  double* F[3];
  F[0] = f;
  F[1] = g;
  F[2] = h;
  #endif
  F[normalNonZero] = fluxOut;
  flux(stateOut, F);

  for (int i=0; i<5; i++) {
    fluxOut[i] = F[normalNonZero][i];
  }

//  //  fluxOut
//  //  //@todo Please implement
//  fluxOut[0] = fluxIn[0];
//  fluxOut[1] = fluxIn[1];
//  fluxOut[2] = fluxIn[2];
//  fluxOut[3] = fluxIn[3];
//  fluxOut[4] = fluxIn[4];
//  //  // stateOut
//  //  // @todo Please implement
//  stateOut[0] = stateIn[0];
//  stateOut[1] = stateIn[1];
//  stateOut[2] = stateIn[2];
//  stateOut[3] = stateIn[3];
//  stateOut[4] = stateIn[4];
}
