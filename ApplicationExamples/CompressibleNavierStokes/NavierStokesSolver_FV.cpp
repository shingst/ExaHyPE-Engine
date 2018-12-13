#include "NavierStokesSolver_FV.h"

#include "NavierStokesSolver_FV_Variables.h"
#include "NavierStokesSolver_ADERDG_Variables.h"

#include "PDE.h"
#include "AMR/Criterion.h"

#include "Scenarios/Atmosphere.h"
#include "SetupHelper.h"

tarch::logging::Log NavierStokes::NavierStokesSolver_FV::_log( "NavierStokes::NavierStokesSolver_FV" );

void NavierStokes::NavierStokesSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  auto parsedConfig = parseConfig(cmdlineargs, constants);
  ns = std::move(parsedConfig.ns);
  scenarioName = std::move(parsedConfig.scenarioName);
  scenario = std::move(parsedConfig.scenario);
}

void NavierStokes::NavierStokesSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    ns.setHeight(Q, x[DIMENSIONS-1]);
    ns.setBackgroundState(Q, 0.0, 0.0);

    AbstractNavierStokesSolver_ADERDG::Variables vars(Q);
    scenario->initialValues(x, ns, vars);
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(Q[i]), i, Q[i]);
    }
  }

  const auto vars = ReadOnlyVariables{Q};
  const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
  assertion5(pressure >= 0, pressure, vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
}

void NavierStokes::NavierStokesSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  ns.evaluateEigenvalues(Q, dIndex, lambda);
}

void NavierStokes::NavierStokesSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // No slip, 2D
  ReadOnlyVariables varsIn(stateInside);
  Variables varsOut(stateOutside);

  // TODO(Lukas) Use gradient for analytical bcs?
  // TODO(Lukas) Support coupling!
  ns.setHeight(stateOutside, x[DIMENSIONS-1]); // TODO(Lukas) Check/refactor

#if DIMENSIONS == 2
  stateOutside[0] = stateInside[0];
  stateOutside[1] = -stateInside[1];
  stateOutside[2] = -stateInside[2];
  stateOutside[3] = stateInside[3];
  if (ns.useAdvection) {
      stateOutside[4] = stateInside[4];
  }

  if (scenario->getBoundaryType(faceIndex) == BoundaryType::movingWall) {
    const auto wallSpeed = 1.0;
    stateOutside[1] = 2 * wallSpeed - stateInside[1];
  }
#else
  stateOutside[0] = stateInside[0];
  stateOutside[1] = -stateInside[1];
  stateOutside[2] = -stateInside[2];
  stateOutside[3] = -stateInside[3];
  stateOutside[4] = stateInside[4];
  if (ns.useAdvection) {
      stateOutside[5] = stateInside[5];
  }
#endif
  if (scenarioName != "two-bubbles") {
    return;
  }

  // TODO(Lukas) Don't do this for every scenario!
  const auto g = 9.81; // TODO(Lukas) Refactor?
  const auto backgroundPotTemperature = 300;
  // TODO(Lukas) Remove
  const auto pressure = computeHydrostaticPressure(ns, g, x[DIMENSIONS-1], backgroundPotTemperature);
  // TODO(Lukas) Maybe choose T = 2WallT - Tin or sth like that
  // Then T at boundary should be 2WallT
  const auto T = potentialTToT(ns, pressure, backgroundPotTemperature);
  const auto rho = pressure / (ns.gasConstant * T);
  // TODO(Lukas) Is this also an accurate reconstruction for advection-scenarios?

  ns.setBackgroundState(stateOutside, rho, pressure); // TODO(Lukas) Is this correct?
  // TODO(Lukas) Is this correct?
  auto E = -1;
  if (ns.useGravity) {
    E = ns.evaluateEnergy(rho, pressure, varsOut.j(), ns.getZ(stateInside), x[DIMENSIONS-1]);
  } else {
    E = ns.evaluateEnergy(rho, pressure, varsOut.j(), ns.getZ(stateInside));
  }
  varsOut.E() = E;
  varsOut.rho() = rho;
}

void NavierStokes::NavierStokesSolver_FV::viscousFlux(const double* const Q,const double* const gradQ, double** F) {
  ns.evaluateFlux(Q, gradQ, F, true);
}

void NavierStokes::NavierStokesSolver_FV::viscousEigenvalues(const double* const Q, const int dIndex, double* lambda) {
  ns.evaluateDiffusiveEigenvalues(Q, dIndex, lambda);
}

//You can either implement this method or modify fusedSource
void NavierStokes::NavierStokesSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // TODO: Actually use coordinates!
  scenario->source(x, t, ns, Q, S);
}

std::vector<double> NavierStokes::NavierStokesSolver_FV::mapGlobalObservables(const double *const Q,
        const tarch::la::Vector<DIMENSIONS,double>& dx) const {
  // TODO(Lukas): Implementation is really slow but should work.
  return mapGlobalObservablesFV<PatchSize, GhostLayerWidth, NumberOfVariables,
          NumberOfVariables, NumberOfGlobalObservables>(
          Q, dx, scenarioName, ns);
}

std::vector<double> NavierStokes::NavierStokesSolver_FV::resetGlobalObservables() const {
    return ::NavierStokes::resetGlobalObservables(NumberOfGlobalObservables);
}

void NavierStokes::NavierStokesSolver_FV::reduceGlobalObservables(
        std::vector<double> &reducedGlobalObservables,
        const std::vector<double> &curGlobalObservables) const {
    ::NavierStokes::reduceGlobalObservables(reducedGlobalObservables, curGlobalObservables,
            NumberOfGlobalObservables);
}

