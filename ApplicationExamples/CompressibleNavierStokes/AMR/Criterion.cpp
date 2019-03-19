#include "Criterion.h"

#include "VarianceHelper.h"
#include "totalVariation.h"

std::vector<double> NavierStokes::resetGlobalObservables(
    int NumberOfGlobalObservables) {
  if (NumberOfGlobalObservables == 0) {
    return {};
  }
  return {-1.0, -1.0, 0};
}

std::vector<double> NavierStokes::mapGlobalObservables(
    const double *const Q, const tarch::la::Vector<DIMENSIONS, double> &dx,
    const std::string &scenarioName, const PDE &ns,
    const AMRSettings &amrSettings, int Order, int NumberOfVariables,
    int NumberOfParameters, int NumberOfGlobalObservables) {
  if (NumberOfGlobalObservables == 0) {
    return {};
  }

  auto indicator = amrSettings.indicator;
  auto useTV = amrSettings.useTotalVariation;

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

  auto observables = resetGlobalObservables(NumberOfGlobalObservables);
  // TODO(Lukas) Implement global observables for 3D!
  // const auto idxQ = kernels::idx3(Order + 1, Order + 1,
  //                                NumberOfVariables + NumberOfParameters);
  assertion(useTV);

  if (useTV) {
    const auto tv =
        totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx,
                       amrSettings.correctForVolume, computeIndicator);
    return {tv, 0, 1};
  }
  throw -1;

  return {-1, -1, -1};  // TODO(Lukas)
}

void NavierStokes::reduceGlobalObservables(
    std::vector<double> &reducedGlobalObservables,
    const std::vector<double> &curGlobalObservables,
    int NumberOfGlobalObservables) {
  if (NumberOfGlobalObservables == 0) {
    return;
  }

  assertion2(reducedGlobalObservables.size() == curGlobalObservables.size(),
             reducedGlobalObservables.size(), curGlobalObservables.size());

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