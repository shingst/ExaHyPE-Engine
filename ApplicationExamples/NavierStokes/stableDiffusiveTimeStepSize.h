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
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>

#include "tarch/la/Vector.h" 
#include "kernels/KernelUtils.h"

// Debugging
#include <iostream>

template <typename SolverType>
double stableDiffusiveTimeStepSize(
    SolverType& solver,
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx) {
  constexpr int numberOfVariables  = SolverType::NumberOfVariables;
  constexpr int numberOfParameters = SolverType::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int order              = SolverType::Order;
  constexpr int basisSize          = order+1;
  constexpr double cflFactor       = 0.7; //SolverType::CFL;

  double dt = std::numeric_limits<double>::max();
  kernels::idx4 idx_luh(basisSize, basisSize, basisSize, numberOfData);

  // Iterate over dofs.
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
	for (int dim = 0; dim < DIMENSIONS; dim++) {
	  // First compute max eigenvalues of hyperbolic and diffusive part.
	  auto hyperbolicEigenvalues = std::array<double,numberOfVariables>();
	  auto diffusiveEigenvalues = std::array<double,numberOfVariables>();
	  solver.eigenvalues(luh + idx_luh(i,j,k,0), dim, hyperbolicEigenvalues.data());
	  solver.diffusiveEigenvalues(luh + idx_luh(i,j,k,0), dim, diffusiveEigenvalues.data());

	  double maxHyperbolicEigenvalue = 0.0;
	  double maxDiffusiveEigenvalue = 0.0;

	  for (const auto eigen : hyperbolicEigenvalues) {
	    maxHyperbolicEigenvalue = std::max(maxHyperbolicEigenvalue, std::abs(eigen));
	  }
	  for (const auto eigen : diffusiveEigenvalues) {
	    maxDiffusiveEigenvalue = std::max(maxDiffusiveEigenvalue, std::abs(eigen));
	  }

	  // const double scale = cflFactor/(2 * order + 1);
	  // // TODO(Lukas) what exactly is the scale here?
	  // const double denominator = maxHyperbolicEigenvalue + 2 * maxDiffusiveEigenvalue * ((2 * order + 1)/dx[dim]);
	  // const double curDt = scale * (dx[dim]/denominator);

	  const double curDt = cflFactor * (dx[dim]/(DIMENSIONS * (2 * order + 1))) *
	    1.0/(maxHyperbolicEigenvalue + maxDiffusiveEigenvalue *
	     ((2* (2 * order + 1))/(dx[dim]))
	     );
	  // std::cout << "invDx = " << curInvDx << " curDt = " << curDt << std::endl;
	  // std::cout << "lam_h = " << maxHyperbolicEigenvalue << " lam_d = " << maxDiffusiveEigenvalue << std::endl;
	  dt = std::min(dt, curDt);
	}
      }
    }
  }
    
  return dt;
}
