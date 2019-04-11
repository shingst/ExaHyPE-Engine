#ifndef COMPRESSIBLENAVIERSTOKES_CRITERION_H
#define COMPRESSIBLENAVIERSTOKES_CRITERION_H

#include <string>
#include <vector>

#include "kernels/limiter/generic/Limiter.h"

#include "AMRSettings.h"
#include "PDE.h"
#include "totalVariation.h"
#include "VarianceHelper.h"

namespace NavierStokes {

/** NEW **/
template <typename GlobalObservables>
void resetGlobalObservables(GlobalObservables& obs)  {
  if (obs.size() != 0) {
    auto *obsRaw = obs.data();
    obsRaw[0] = -1.0;
    obsRaw[1] = -1.0;
    obsRaw[2] = 0.0;
  }
}

template <
  typename GlobalObservables,
  typename ReadOnlyGlobalObservables
>
void mergeGlobalObservables(
    GlobalObservables&         obs,
    ReadOnlyGlobalObservables& other)  {
    if (obs.size() == 0) {
      return;
    }

  assertion2(obs.size() == other.size(),
             obs.size(), other.size());

  auto *reducedGlobalObservables = obs.data();
  auto const *curGlobalObservables = other.data();
  
  const auto mean0 = reducedGlobalObservables[0];
  const auto mean1 = curGlobalObservables[0];
  const auto var0 = reducedGlobalObservables[1];
  const auto var1 = curGlobalObservables[1];
  const auto count0 = reducedGlobalObservables[2];
  const auto count1 = curGlobalObservables[2];

  auto mergedMean = 0.0;
  auto mergedVariance = 0.0;

  std::tie(mergedMean, mergedVariance) =
      mergeVariance(mean0, mean1, var0, var1, count0, count1);
  reducedGlobalObservables[0] = mergedMean;
  reducedGlobalObservables[1] = mergedVariance;
  reducedGlobalObservables[2] = count0 + count1;
    
  }

template <typename GlobalObservables>
void mapGlobalObservablesDG(
    GlobalObservables& globalObservables,
    const double* const Q,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const std::string &scenarioName,
    const PDE &ns,
    const AMRSettings &amrSettings,
    int Order,
    int NumberOfVariables,
    int NumberOfParameters)  {
  
  if (globalObservables.size() == 0) {
    return;
  }

  auto indicator = amrSettings.indicator;
  auto useTV = amrSettings.useTotalVariation;

  indicator = IndicatorVariable::potentialTemperature;
  useTV = true;
  assert(indicator == IndicatorVariable::potentialTemperature);
  assert(useTV);

  auto computeIndicator = [&](const double *const Q) {
    const auto vars = ReadOnlyVariables{Q};
    if (indicator == IndicatorVariable::potentialTemperature) {
      const auto pressure =
          ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q),
                  ns.getHeight(Q));
      const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
      return ns.evaluatePotentialTemperature(temperature, pressure);
    }
    if (indicator == IndicatorVariable::Z) {
      return ns.getZ(Q)/vars.rho();
    }
    auto backgroundRho = 0.0;
    auto backgroundPressure = 0.0;
    std::tie(backgroundRho, backgroundPressure) = ns.getBackgroundState(Q);


    if (indicator == IndicatorVariable::rho) {
      return vars.rho() - backgroundRho;
    } else {
      // Pressure
      return ns.evaluatePressure(vars.E(), vars.rho(), vars.j(),
              ns.getZ(Q), ns.getHeight(Q)) -
             backgroundPressure;
    }
  };

  auto *observables = globalObservables.data();
  // TODO(Lukas) Implement global observables for 3D!

  assertion(useTV);

  if (useTV) {
    const auto tv =
        totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx,
                       amrSettings.correctForVolume, computeIndicator);
    observables[0] = tv;
    observables[1] = 0;
    observables[2] = 1;
  } else {
    throw -1;
  }
}

 template <typename Solver, typename GlobalObservables>
void mapGlobalObservablesFV(
    Solver *solver,			    
    GlobalObservables& globalObservables,
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize)  {
   if (globalObservables.size() == 0) {
     return;
   }
   constexpr unsigned int NumberOfVariables = Solver::NumberOfVariables;
   constexpr unsigned int NumberOfParameters = Solver::NumberOfParameters;
   constexpr auto numberOfData = NumberOfVariables + NumberOfParameters;
   constexpr unsigned int PatchSize = Solver::PatchSize;
   constexpr unsigned int GhostLayerWidth = Solver::GhostLayerWidth;

   constexpr auto variablesPerPatch = (PatchSize+2*GhostLayerWidth)*(PatchSize+2*GhostLayerWidth)*numberOfData;
   constexpr int patchBegin = GhostLayerWidth; // patchBegin cell is inside domain
   constexpr int patchEnd = patchBegin+PatchSize; // patchEnd cell is outside domain
   kernels::idx3 idx_slope(PatchSize+2*GhostLayerWidth,
			   PatchSize+2*GhostLayerWidth,
			   DIMENSIONS);

   double slopes[variablesPerPatch*DIMENSIONS] = {0.0};

   totalVariationFV(solver, luh, cellSize, slopes);
   
   double meanReduced = -1.0;
   double varReduced = -1.0;
   double n = 0.0;
   
   for (int i = patchBegin - 1; i < patchEnd; ++i) {
     for (int j = patchBegin-1; j < patchEnd; ++j) {
       for (int d = 0; d < DIMENSIONS; ++d) {
	 const auto curTv = slopes[idx_slope(i,j,d)];
	 std::tie(meanReduced, varReduced) = mergeVariance(meanReduced,
							   curTv,
							   varReduced,
							   0.0, n,
							   1);
	 n += 1;
       }
     }
   }



   auto *rawGlobalObservables = globalObservables.data();
   rawGlobalObservables[0] = meanReduced;
   rawGlobalObservables[1] = varReduced;
   rawGlobalObservables[2] = n;
 }

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_CRITERION_H
