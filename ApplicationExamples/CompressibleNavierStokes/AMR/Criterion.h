#ifndef COMPRESSIBLENAVIERSTOKES_CRITERION_H
#define COMPRESSIBLENAVIERSTOKES_CRITERION_H

#include <string>
#include <vector>

#include "PDE.h"

namespace NavierStokes {
std::vector<double> resetGlobalObservables(int NumberOfGlobalObservables);

std::vector<double> mapGlobalObservables(
    const double *const Q, const tarch::la::Vector<DIMENSIONS, double> &dx,
    const std::string &scenarioName, const PDE &ns, int Order,
    int NumberOfVariables, int NumberOfParameters,
    int NumberOfGlobalObservables);

void reduceGlobalObservables(std::vector<double> &reducedGlobalObservables,
                             const std::vector<double> &curGlobalObservables,
                             int NumberOfGlobalObservables);
}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_CRITERION_H
