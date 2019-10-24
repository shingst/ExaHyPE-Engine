#ifndef COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#define COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#include <vector>

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

#include "PDE.h"
#include "kernels/aderdg/generic/Kernels.h"

using computeFunc = double (*)(const double* const);

double totalVariation(const double* Q, int order, int numberOfVariables,
                      int numberOfParameters,
                      const tarch::la::Vector<DIMENSIONS, double>& dx,
                      bool correctForVolume);

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

inline double forwardDiff(double l, double c, double h) {
  assert(std::isfinite( (1./h) * (l - c)));
  return (1./h) * (l - c);
}

inline double centralDiff(double l, double r, double h) {
  assert(std::isfinite( (1./(2*h)) * (r - l)));
  return (1./(2 * h)) * (r - l);
}

inline double backwardDiff(double c, double r, double h) {
  assert(std::isfinite( (1./h) * (c - r)));
  return (1./h) * (c - r);
}

inline double stableDiff(double l,
		  double c,
		  double r,
		  int idxC,  
		  double h,
		  size_t ghostLayerWidth,
		  size_t patchSize) {
  // TODO(Lukas) Check for off by one errors
  // TODO(Lukas) Use central differences for center?
  // TODO(Lukas) Check if order of arguments is correct.
  
  // idxC must be signed to avoid underflow!
  const auto idxL = idxC - 1;
  const auto idxR = idxC + 1;
  
  //std::cout << idxC << std::endl;
  if (idxL <= static_cast<int>(ghostLayerWidth)) {
    return backwardDiff(c, r, h);
  } else if (idxR >= static_cast<int>(patchSize + ghostLayerWidth)) {
    return forwardDiff(l, c, h);
  } else {
    return centralDiff(l, r, h);
  }
}

template <typename Solver>
void totalVariationFV(Solver* solver,
			     const double* const luh,
			     const tarch::la::Vector<DIMENSIONS, double>& cellSize,
			     double *slope) {
  constexpr unsigned int NumberOfVariables = Solver::NumberOfVariables;
  constexpr unsigned int NumberOfParameters = Solver::NumberOfParameters;
  constexpr unsigned int PatchSize = Solver::PatchSize;
  constexpr unsigned int GhostLayerWidth = Solver::GhostLayerWidth;
  const auto &ns = solver->ns;
  // Ignore efficiency for now.
  // TODO: Are derivatives correct?
  // TODO: Is ghost layer handled correctly?
  // TODO: Don't hardcore indicator variable.
  assert(DIMENSIONS == 2);
  constexpr auto numberOfData = NumberOfVariables + NumberOfParameters;
  kernels::idx3 idx(PatchSize+2*GhostLayerWidth,PatchSize+2*GhostLayerWidth,numberOfData);
  kernels::idx2 idx_obs(PatchSize+2*GhostLayerWidth,PatchSize+2*GhostLayerWidth);

  kernels::idx3 idx_slope(PatchSize+2*GhostLayerWidth,
			  PatchSize+2*GhostLayerWidth,
			  DIMENSIONS);
  const auto subcellSize = (1./PatchSize) * cellSize; 
  assert(std::isfinite(subcellSize[0]));
  assert(std::isfinite(subcellSize[1]));

  // Compute slope by finite differences
  constexpr auto variablesPerPatch = (PatchSize+2*GhostLayerWidth)*(PatchSize+2*GhostLayerWidth)*numberOfData;
  constexpr int patchBegin = GhostLayerWidth; // patchBegin cell is inside domain
  constexpr int patchEnd = patchBegin+PatchSize; // patchEnd cell is outside domain
  auto observables = std::vector<double>((PatchSize+2*GhostLayerWidth)*(PatchSize+2*GhostLayerWidth));

  auto computeIndicator = [&](const double *const Q) {
    const auto vars = ReadOnlyVariables(Q);
    const auto pressure =
    ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q),
			ns.getHeight(Q));
    const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
    return ns.evaluatePotentialTemperature(temperature, pressure);
  };

  // Compute indicator variables
  for (int j = patchBegin; j < patchEnd; j++) {
    for (int k = patchBegin; k < patchEnd; k++) {
      observables[idx_obs(j,k)] = computeIndicator(luh + idx(j,k,0));
    }
  }

  
  // Compute slopes    
  // slopex
  for (int j = patchBegin; j < patchEnd; j++) { // y
    for (int k = patchBegin; k < patchEnd; k++) { // x
      const auto left = observables[idx_obs(j, k-1)];
      const auto center = observables[idx_obs(j, k+0)];
      const auto right = observables[idx_obs(j, k+1)];

      slope[idx_slope(j, k, 0)] =
        stableDiff(left, center, right,
                   j,
                   subcellSize[0],
                   GhostLayerWidth,
                   PatchSize);
    }
  }
  // slopey
  for (int j = patchBegin; j < patchEnd; j++) { // y
    for (int k = patchBegin; k < patchEnd; k++) { // x
      const auto left = observables[idx_obs(j-1, k)];
      const auto center = observables[idx_obs(j, k)];
      const auto right = observables[idx_obs(j+1, k)];

      slope[idx_slope(j, k, 1)] =
        stableDiff(left, center, right,
                   k,
                   subcellSize[1],
                   GhostLayerWidth,
                   PatchSize);

    }
  }
}

#endif  // COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
