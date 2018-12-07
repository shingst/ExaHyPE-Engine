#include "Scenario.h"
#include <stdexcept>

NavierStokes::Scenario::Scenario() {}

void NavierStokes::Scenario::initialValues(const double* const x, const PDE& ns,
                                           Variables& vars) {
  // Throw away gradient here, this is not efficient but correct.
  auto gradState = std::array<double, DIMENSIONS * vars.SizeVariables>();
  analyticalSolution(x, 0.0, ns, vars, gradState.data());

  assertion2(vars.E() > 0.0, vars.E(), x)
}

void NavierStokes::Scenario::analyticalSolution(const double* const x,
                                                const double t, const PDE& ns,
                                                Variables& vars,
                                                double* gradState) {
  throw std::logic_error("Analytical solution not implemented!");
}
void NavierStokes::Scenario::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  constexpr auto NumberOfVariables = DIMENSIONS + 2;  // TODO(Lukas) generalise?
  auto source = ReadOnlyVariables(S);
  std::fill_n(S, source.SizeVariables, 0.0);
}

NavierStokes::BoundaryType NavierStokes::Scenario::getBoundaryType(int faceId) {
  return BoundaryType::wall;
}

const double NavierStokes::Scenario::getGamma() const { return gamma; }

const double NavierStokes::Scenario::getPr() const { return Pr; }

const double NavierStokes::Scenario::getC_v() const { return c_v; }

const double NavierStokes::Scenario::getC_p() const { return c_p; }

const double NavierStokes::Scenario::getGasConstant() const {
  return gasConstant;
}

const double NavierStokes::Scenario::getReferencePressure() const {
  return referencePressure;
}

double NavierStokes::Scenario::getMolecularDiffusionCoeff() const {
  return 0.0; // TODO(Lukas)
}
double NavierStokes::Scenario::getQ0() const {
  return 0.0; // TODO(Lukas)
}

bool NavierStokes::Scenario::getUseAdvection() const {
  return false;
}

