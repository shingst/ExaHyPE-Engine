#ifndef COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#define COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#include <vector>

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"
static_assert(DIMENSIONS, "Total variation only supported for 2D!");
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

  const auto idxQ = kernels::idx3(basisSize, basisSize, numberOfData);
  const auto idxObservable = kernels::idx3(basisSize, basisSize, 1);

  auto observable = std::vector<double>(std::pow(basisSize, DIMENSIONS));
  for (int k = 0; k < basisSize; k++) {
    for (int l = 0; l < basisSize; l++) {
      observable[idxObservable(k, l, 0)] = mapObservable(Q + idxQ(k, l, 0));
    }
  }

  return totalVariation(observable.data(), order, 1, 0, dx, correctForVolume);
}
#endif  // COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
