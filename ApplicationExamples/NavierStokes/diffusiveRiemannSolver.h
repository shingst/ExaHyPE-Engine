#ifndef __DIFFUSIVE_RIEMANN_SOLVER_HEADER__
#define __DIFFUSIVE_RIEMANN_SOLVER_HEADER__

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

// included in ../../Kernels.h

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "kernels/KernelUtils.h"
#include "kernels/GaussLegendreQuadrature.h"

/**
* We implement a very simple Rusanov scheme with scalar dissipation
* (smax*Id).
*
#ifndef __DIFFUSIVE_RIEMANN_SOLVER_HEADER__
#define __DIFFUSIVE_RIEMANN_SOLVER_HEADER__

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

// included in ../../Kernels.h

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "kernels/KernelUtils.h"
#include "kernels/GaussLegendreQuadrature.h"

/**
* We implement a very simple Rusanov scheme with scalar dissipation
* (smax*Id).
*
* We need to consider material parameters
* in QL and QR.
* We don't need to consider material parameters
* in FL,FR.
*/
//TODO(Lukas) Accept vector for characteristicLength
template <bool useNCP, typename SolverType>
void riemannSolverNonlinear(SolverType& solver,
                          double* FL, double* FR, const double* const QL,
                          const double* const QR,
			    const tarch::la::Vector<DIMENSIONS, double>& characteristicLength,
                          const double dt,
                          const int direction) {
  constexpr int numberOfVariables  = SolverType::NumberOfVariables;
  constexpr int numberOfParameters = SolverType::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int order              = SolverType::Order;
  constexpr int basisSize          = order+1;

  using namespace kernels;
  
  // Compute the average variables and parameters from the left and the right
  double QavL[numberOfData] = {0.0};
  double QavR[numberOfData] = {0.0};
  {
    idx3 idx_QLR(basisSize, basisSize, numberOfData);
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        const double weight =
            kernels::gaussLegendreWeights[order][i] *
            kernels::gaussLegendreWeights[order][j];

        for (int k = 0; k < numberOfData; k++) {
          QavL[k] += weight * QL[idx_QLR(i, j, k)];
          QavR[k] += weight * QR[idx_QLR(i, j, k)];
        }
      }
    }
  }

  // compute fluxes (and fluctuations for non-conservative PDEs)
  double Qavg[numberOfData];
  idx2 idx_gradQ(DIMENSIONS, numberOfVariables);
  double gradQ[DIMENSIONS][numberOfVariables] = {0.0};
  double ncp[numberOfVariables]               = {0.0};
  {
    idx3 idx_FLR(basisSize, basisSize, numberOfVariables);
    idx3 idx_QLR(basisSize, basisSize, numberOfData);

    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
	double LL[numberOfVariables] = {0.0}; // do not need to store material parameters
	double LR[numberOfVariables] = {0.0};
	solver.eigenvalues(&QL[idx_QLR(i,j,0)], direction, LL);
	solver.eigenvalues(&QR[idx_QLR(i,j,0)], direction, LR);

	// skip parameters
	// hyperbolic eigenvalues
	std::transform(LL, LL + numberOfVariables, LL, std::abs<double>);
	std::transform(LR, LR + numberOfVariables, LR, std::abs<double>);
	const double smax_L = *std::max_element(LL, LL + numberOfVariables);
	const double smax_R = *std::max_element(LR, LR + numberOfVariables);
	const double maxHyperbolicEigenvalue = std::max(smax_L, smax_R);

	// diffusive eigenvalues
	solver.diffusiveEigenvalues(&QL[idx_QLR(i,j,0)], direction, LL);
	solver.diffusiveEigenvalues(&QR[idx_QLR(i,j,0)], direction, LR);
	std::transform(LL, LL + numberOfVariables, LL, std::abs<double>);
	std::transform(LR, LR + numberOfVariables, LR, std::abs<double>);
	const double smaxDiffusive_L = *std::max_element(LL, LL + numberOfVariables);
	const double smaxDiffusive_R = *std::max_element(LR, LR + numberOfVariables);
	const double maxDiffusiveEigenvalue = std::max(smaxDiffusive_L, smaxDiffusive_R);

	//constexpr double pi = std::acos(-1.0);
	//const double factor =  (2 * order + 1)/(characteristicLength + std::sqrt(0.5 * pi))
	const double factor = 2 * ((order + 1)/characteristicLength[direction]);
	const double penalty = maxHyperbolicEigenvalue + factor * maxDiffusiveEigenvalue;

	//if(useNCP) { // we don't use matrixB but the NCP call here.
        //  for(int l=0; l < numberOfVariables; l++) {
        //    gradQ[direction][l] = QR[idx_QLR(i, j, l)] - QL[idx_QLR(i, j, l)];
        //    Qavg[l] = 0.5 * (QR[idx_QLR(i, j, l)] + QL[idx_QLR(qi, j, l)]);
        //  }

        //  solver.nonConservativeProduct(Qavg, gradQ[0], ncp);
        //}

        // skip parameters
        for (int k = 0; k < numberOfVariables; k++) {
          FL[idx_FLR(i, j, k)] =
              0.5 * (FR[idx_FLR(i, j, k)] + FL[idx_FLR(i, j, k)]) -
              0.5 * penalty * (QR[idx_QLR(i, j, k)] - QL[idx_QLR(i, j, k)]);
	  assertion7(std::isfinite(FL[idx_FLR(i,j,k)]), i, j, k, factor, penalty, maxHyperbolicEigenvalue, maxDiffusiveEigenvalue);

          //if(useNCP) {
          //  FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] - 0.5 * ncp[k];
          //  FL[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] + 0.5 * ncp[k];
          //} else {
            FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)];
	//}
	}
      }
    }
  }
}

#endif
