// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
// ==============================================
// Please do not change the implementations below
// ==============================================
#include "NavierStokesSolver.h"

#include "kernels/limiter/generic/Limiter.h"

#include "AMR/Criterion.h"

NavierStokes::NavierStokesSolver::NavierStokesSolver(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int limiterHelperLayers,
        const int DMPObservables,
        const double DMPRelaxationParameter,
        const double DMPDifferenceScaling,
        const int iterationsToCureTroubledCell 
) :
  exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
      "NavierStokesSolver",
    new NavierStokes::NavierStokesSolver_ADERDG(
      maximumMeshSize,maximumMeshDepth,haloCells,regularisedFineGridLevels,timeStepping,limiterHelperLayers,DMPObservables),
    new NavierStokes::NavierStokesSolver_FV(
      maximumMeshSize, timeStepping),
    DMPRelaxationParameter,
    DMPDifferenceScaling,
    iterationsToCureTroubledCell) {}

void NavierStokes::NavierStokesSolver::projectOnFVLimiterSpace(const double* const luh, double* const lim) const {
  kernels::limiter::generic::c::projectOnFVLimiterSpace<Order+1,NumberOfVariables+NumberOfParameters,GhostLayerWidth>(luh, lim);
}

void NavierStokes::NavierStokesSolver::projectOnDGSpace(const double* const lim, double* const luh) const {
  kernels::limiter::generic::c::projectOnDGSpace<Order+1,NumberOfVariables+NumberOfParameters,GhostLayerWidth>(lim, luh);
}

bool NavierStokes::NavierStokesSolver::discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* boundaryMinPerVariables, double* boundaryMaxPerVariables) {
  return kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch<AbstractNavierStokesSolver_ADERDG, NumberOfDMPObservables, GhostLayerWidth>(luh, *static_cast<AbstractNavierStokesSolver_ADERDG*>(_solver.get()), _DMPMaximumRelaxationParameter, _DMPDifferenceScaling, boundaryMinPerVariables, boundaryMaxPerVariables);
}

void NavierStokes::NavierStokesSolver::findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) {
  kernels::limiter::generic::c::findCellLocalMinAndMax<AbstractNavierStokesSolver_ADERDG, NumberOfDMPObservables>(luh, *static_cast<AbstractNavierStokesSolver_ADERDG*>(_solver.get()), localMinPerVariables, localMaxPerVariable);
}
void NavierStokes::NavierStokesSolver::findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) {
  kernels::limiter::generic::c::findCellLocalLimiterMinAndMax<AbstractNavierStokesSolver_ADERDG, NumberOfDMPObservables, GhostLayerWidth>(lim, *static_cast<AbstractNavierStokesSolver_ADERDG*>(_solver.get()), localMinPerObservable,localMaxPerObservable);
}

// TODO:
std::vector<double> NavierStokes::NavierStokesSolver::mapGlobalObservables(const double* const Q, const tarch::la::Vector<DIMENSIONS, double> &dx) const {
  return {};
}

std::vector<double> NavierStokes::NavierStokesSolver::resetGlobalObservables() const {
  return {};
}

void NavierStokes::NavierStokesSolver::reduceGlobalObservables(
            std::vector<double>& reducedGlobalObservables,
            const std::vector<double>& curGlobalObservables) const { }
