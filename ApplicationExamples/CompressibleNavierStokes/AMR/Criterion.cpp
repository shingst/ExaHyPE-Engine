#include "Criterion.h"

#include "VarianceHelper.h"
#include "totalVariation.h"

std::vector<double> NavierStokes::resetGlobalObservables(
    int NumberOfGlobalObservables) {
  if (NumberOfGlobalObservables == 0) return {};
  return {-1.0, -1.0, 0};
}

std::vector<double> NavierStokes::mapGlobalObservables(
    const double *const Q, const tarch::la::Vector<DIMENSIONS, double> &dx,
    const std::string &scenarioName, const PDE &ns, int Order,
    int NumberOfVariables, int NumberOfParameters,
    int NumberOfGlobalObservables) {
  if (NumberOfGlobalObservables == 0) return {};

  auto observables = resetGlobalObservables(NumberOfGlobalObservables);
  // TODO(Lukas) Implement global observables for 3D!
  const auto idxQ = kernels::idx3(Order + 1, Order + 1,
                                  NumberOfVariables + NumberOfParameters);

  auto tv = 0.0;
  if (scenarioName == "two-bubbles" || scenarioName == "density-current" ||
      scenarioName == "coupling-test") {
    auto computePotT = [&](const double *const Q) {
      const auto vars = ReadOnlyVariables{Q};
      const auto pressure =
          ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
      const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
      return ns.evaluatePotentialTemperature(temperature, pressure);
    };

    tv = totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx,
                        false, computePotT);
  } else {
    auto computeIndicator = [&](const double *const Q) {
      const auto vars = ReadOnlyVariables{Q};
      const auto pressure =
          ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
      return vars.E();
    };
    tv = totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx,
                        false, computeIndicator);
  }

  return {tv, 0, 1};
}

void NavierStokes::reduceGlobalObservables(
    std::vector<double> &reducedGlobalObservables,
    const std::vector<double> &curGlobalObservables,
    int NumberOfGlobalObservables) {
  if (NumberOfGlobalObservables == 0) return;

  assertion2(reducedGlobalObservables.size() == curGlobalObservables.size(),
             reducedGlobalObservables.size(), curGlobalObservables.size());

  const auto mean0 = reducedGlobalObservables[0];
  const auto mean1 = curGlobalObservables[0];
  const auto var0 = reducedGlobalObservables[1];
  const auto var1 = curGlobalObservables[1];
  const auto count0 = static_cast<int>(reducedGlobalObservables[2]);
  const auto count1 = static_cast<int>(curGlobalObservables[2]);

  auto mergedMean = 0.0;
  auto mergedVariance = 0.0;

  std::tie(mergedMean, mergedVariance) =
      mergeVariance(mean0, mean1, var0, var1, count0, count1);
  reducedGlobalObservables[0] = mergedMean;
  reducedGlobalObservables[1] = mergedVariance;
  reducedGlobalObservables[2] = count0 + count1;
}