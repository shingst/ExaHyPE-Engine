#include "SetupHelper.h"
#include "exahype/parser/ParserView.h"

NavierStokes::ScenarioConfig NavierStokes::parseConfig(
    const std::vector<std::string>& cmdlineargs,
    const exahype::parser::ParserView& constants) {
  // TODO(Lukas) Pass from solver, otherwise DG/FV have to agree!
  constexpr auto NumberOfVariables =
      AbstractNavierStokesSolver_ADERDG::NumberOfVariables;
  constexpr auto NumberOfParameters =
      AbstractNavierStokesSolver_ADERDG::NumberOfParameters;

  assert(constants.isValueValidString("scenario"));

  double referenceViscosity;
  if (constants.isValueValidString("viscosity") &&
      constants.getValueAsString("viscosity") == "default") {
    throw - 1;
  } else {
    assert(constants.isValueValidDouble("viscosity"));
    referenceViscosity = constants.getValueAsDouble("viscosity");
  }

  auto scenarioName = constants.getValueAsString("scenario");
  auto scenario = ScenarioFactory::createScenario(scenarioName);

  const bool useGravity =
      constants.getValueAsBoolOrDefault("use-gravity", false);
  const bool useBackgroundState =
      constants.getValueAsBoolOrDefault("use-background-state", false);
  assert(useGravity ||
         !useBackgroundState);  // Background state only works with gravity.

  // Make sure we have the correct number of variables/parameters:
  auto numberOfNecessaryVariables = 1 + DIMENSIONS + 1;
  if (scenario->getUseAdvection()) {
    ++numberOfNecessaryVariables;
  }
  if (NumberOfVariables != numberOfNecessaryVariables) {
    throw - 1;
  }

  auto NumberOfNecessaryParameters = 0;
  if (useGravity) {
    ++NumberOfNecessaryParameters;
  }
  if (useBackgroundState) {
    NumberOfNecessaryParameters += 2;
  }
  if (NumberOfParameters != NumberOfNecessaryParameters) {
    throw - 1;
  }

  auto ns = PDE(referenceViscosity, *scenario, useGravity, useBackgroundState);

  return {ns, scenarioName, std::move(scenario)};
}