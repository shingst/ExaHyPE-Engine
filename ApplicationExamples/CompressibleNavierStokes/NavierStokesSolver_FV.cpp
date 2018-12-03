#include "NavierStokesSolver_FV.h"

#include "NavierStokesSolver_FV_Variables.h"
#include "NavierStokesSolver_ADERDG_Variables.h"
#include "PDE.h"

#include "Scenarios/ScenarioFactory.h"


tarch::logging::Log NavierStokes::NavierStokesSolver_FV::_log( "NavierStokes::NavierStokesSolver_FV" );

void NavierStokes::NavierStokesSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // TODO(Lukas) Refactor init!
   assert(constants.isValueValidString("scenario"));

  double referenceViscosity;
  if (constants.isValueValidString("viscosity") &&
      constants.getValueAsString("viscosity") == "default") {
   throw -1;
  } else {
    assert(constants.isValueValidDouble("viscosity"));
    referenceViscosity = constants.getValueAsDouble("viscosity");
  }

  scenarioName = constants.getValueAsString("scenario");
  scenario = ScenarioFactory::createScenario(scenarioName);

  const auto molecularDiffusionCoeff = scenario->getMolecularDiffusionCoeff();
  auto numberOfNecessaryVariables =
          1 + DIMENSIONS + 1;
  if (scenario->getUseAdvection()) {
    ++numberOfNecessaryVariables;
  }
  if (NumberOfVariables != numberOfNecessaryVariables) {
    throw -1;
  }

  ns = PDE(referenceViscosity, *scenario);
}

void NavierStokes::NavierStokesSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    AbstractNavierStokesSolver_ADERDG::Variables vars(Q);
    scenario->initialValues(x, ns, vars);
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(Q[i]), i, Q[i]);
    }
  }
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
  // TODO(Lukas) Use gradient for analytical bcs?
  // TODO(Lukas) Support coupling!
#if DIMENSIONS == 2
  stateOutside[0] = stateInside[0];
  stateOutside[1] = -stateInside[1];
  stateOutside[2] = -stateInside[2];
  stateOutside[3] = stateInside[3];

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
#endif
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


