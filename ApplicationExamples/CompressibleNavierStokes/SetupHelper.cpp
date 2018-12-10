#include "SetupHelper.h"
#include "exahype/parser/ParserView.h"

NavierStokes::ScenarioConfig NavierStokes::parseConfig(
    const std::vector<std::string>& cmdlineargs,
    const exahype::parser::ParserView& constants) {
  constexpr auto NumberOfVariables =
      AbstractNavierStokesSolver_ADERDG::NumberOfVariables;

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

  auto numberOfNecessaryVariables = 1 + DIMENSIONS + 1;
  if (scenario->getUseAdvection()) {
    ++numberOfNecessaryVariables;
  }
  if (NumberOfVariables != numberOfNecessaryVariables) {
    throw - 1;
  }

  auto ns = PDE(referenceViscosity, *scenario);

  return {ns, scenarioName, std::move(scenario)};
}