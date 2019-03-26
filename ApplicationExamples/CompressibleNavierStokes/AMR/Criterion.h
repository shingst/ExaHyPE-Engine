#ifndef COMPRESSIBLENAVIERSTOKES_CRITERION_H
#define COMPRESSIBLENAVIERSTOKES_CRITERION_H

#include <string>
#include <vector>

#include "PDE.h"

#include "AMRSettings.h"
#include "kernels/limiter/generic/Limiter.h"

namespace NavierStokes {

/** NEW **/
template <typename GlobalObservables>
void resetGlobalObservables(GlobalObservables& obs)  {
  obs.gobs(0)=-1.0;
  obs.gobs(1)=-1.0;
  obs.gobs(2)=0.0;
}

template <
  typename GlobalObservables,
  typename ReadOnlyGlobalObservables
>
void mergeGlobalObservables(
    GlobalObservables&         obs,
    ReadOnlyGlobalObservables& other)  {
  obs.gobs(0)=std::max( obs.gobs(0), other.gobs(0) ); 
  obs.gobs(1)=std::max( obs.gobs(1), other.gobs(1) );
  obs.gobs(2)=std::max( obs.gobs(2), other.gobs(2) );
}

template <typename GlobalObservables>
void mapGlobalObservablesDG(
    GlobalObservables&                          globalObservables,
    const double* const                         luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize)  {
  // TODO(Lukas): Please implement
}

template <typename GlobalObservables>
void mapGlobalObservablesFV(
    GlobalObservables&                          globalObservables,
    const double* const                         luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize)  {
  // TODO(Lukas): Please implement
}
/** NEW - END **/

std::vector<double> resetGlobalObservables(int NumberOfGlobalObservables);

std::vector<double> mapGlobalObservables(
    const double* const Q, const tarch::la::Vector<DIMENSIONS, double>& dx,
    const std::string& scenarioName, const PDE& ns,
    const AMRSettings& amrSettings, int Order, int NumberOfVariables,
    int NumberOfParameters, int NumberOfGlobalObservables);

void reduceGlobalObservables(std::vector<double>& reducedGlobalObservables,
                             const std::vector<double>& curGlobalObservables,
                             int NumberOfGlobalObservables);

template <int PatchSize, int GhostLayerWidth, int NumberOfVariables,
          int NumberOfParameters, int NumberOfGlobalObservables>
std::vector<double> mapGlobalObservablesFV(
    const double* const Q, const tarch::la::Vector<DIMENSIONS, double>& dx,
    const std::string& scenarioName, const PDE& ns,
    const AMRSettings& amrSettings) {
  // Idea: Map FV-Data to DG-Representation, and compute global variables there.
  // This way, we can use the same implementation of total variance for both
  // solvers. Additionally, we have the same TV before and after limiting. It is
  // really slow though, so it might be a good idea to approximate TV
  // differently.
  // TODO(Lukas) Skip mapping if not using TV, e.g. non-gradient based crit.
  if (NumberOfGlobalObservables == 0) {
    return {};
  }

  assert((PatchSize - 1) % 2 == 0);
  constexpr int Order = (PatchSize - 1) / 2;
  constexpr int NumberOfData = NumberOfVariables + NumberOfParameters;
  // TODO(Lukas) Stack allocate this?
  auto QDG =
      std::vector<double>(std::pow(Order + 1, DIMENSIONS) * NumberOfData);
  kernels::limiter::generic::c::projectOnDGSpace<Order + 1,
       NumberOfData, GhostLayerWidth>(Q, QDG.data());

  return ::NavierStokes::mapGlobalObservables(
      QDG.data(), dx, scenarioName, ns, amrSettings, Order, NumberOfVariables,
      NumberOfParameters, NumberOfGlobalObservables);
}

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_CRITERION_H
