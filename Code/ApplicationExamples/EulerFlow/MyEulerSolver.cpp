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

#include <memory>

Euler::MyEulerSolver::MyEulerSolver(int kernelNumber, int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler)
  : exahype::solvers::ADERDGSolver(
      "MyEulerSolver", 5, 0, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}

void Euler::MyEulerSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)

  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  double* f = F[0];
  double* g = F[1];

  // @todo Please implement
  // f
  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);
  // g
  // @todo Please implement
  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);
  
#if DIMENSIONS == 3
  double* h = F[2];
  // h
  // @todo Please implement
  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
#endif
}



void Euler::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2/3
  // Number of variables    = 5 (#unknowns + #parameters)
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];

#if DIMENSIONS == 2
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);
#else
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);
#endif

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}



bool Euler::MyEulerSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}



void Euler::MyEulerSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    const double GAMMA = 1.4;

    Q[0] = 1.;
    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
#if DIMENSIONS == 2
    Q[4] =
        1. / (GAMMA - 1) +
        std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
                 (0.05 * 0.05)) *
            1.0e-3;
#else
    Q[4] =
        1. / (GAMMA - 1) +
        std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) 
                  + (x[2] - 0.5) * (x[2] - 0.5)) / (0.05 * 0.05 * 0.05)) *
            1.0e-3;
#endif  
  }
}



exahype::solvers::Solver::RefinementControl Euler::MyEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}



