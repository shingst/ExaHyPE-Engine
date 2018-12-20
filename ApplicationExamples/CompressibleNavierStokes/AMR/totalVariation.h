#ifndef COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#define COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#include <vector>

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"
using computeFunc = double (*)(const double* const);

double totalVariation(const double* Q, int order, int numberOfVariables,
                      int numberOfParameters,
                      const tarch::la::Vector<DIMENSIONS, double>& dx,
                      bool correctForVolume);
/*
double totalVariation(const double* Q, int order, int numberOfVariables, int
numberOfParameters, const tarch::la::Vector<DIMENSIONS,double>& dx, bool
correctForVolume, computeFunc mapObservable);
*/

template <typename F>
double totalVariation(const double* Q, int order, int numberOfVariables,
                      int numberOfParameters,
                      const tarch::la::Vector<DIMENSIONS, double>& dx,
                      bool correctForVolume, F mapObservable) {
  const auto basisSize = order + 1;
  const auto numberOfData = numberOfVariables + numberOfParameters;

  auto observable = std::vector<double>(std::pow(basisSize, DIMENSIONS));
#if DIMENSIONS == 2
  const auto idxQ = kernels::idx3(basisSize, basisSize, numberOfData);
  const auto idxObservable = kernels::idx3(basisSize, basisSize, 1);
  for (int k = 0; k < basisSize; k++) {
    for (int l = 0; l < basisSize; l++) {
      observable[idxObservable(k, l, 0)] = mapObservable(Q + idxQ(k, l, 0));
    }
  }
#else
  const auto idxQ =
      kernels::idx4(basisSize, basisSize, basisSize, numberOfData);
  const auto idxObservable = kernels::idx4(basisSize, basisSize, basisSize, 1);
  for (int j = 0; j < basisSize; j++) {
    for (int k = 0; k < basisSize; k++) {
      for (int l = 0; l < basisSize; l++) {
        observable[idxObservable(j, k, l, 0)] =
            mapObservable(Q + idxQ(j, k, l, 0));
      }
    }
  }
#endif

  return totalVariation(observable.data(), order, 1, 0, dx, correctForVolume);
}
#endif  // COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
