#include "Scenario.h"
#include <stdexcept>

NavierStokes::Scenario::Scenario() {}

void NavierStokes::Scenario::initialValues(const double* const x,
                                           const NavierStokes& ns,
                                           Variables& vars) {
  // Throw away gradient here, this is not efficient but correct.
  auto gradState = std::array<double, DIMENSIONS * vars.SizeVariables>();
  analyticalSolution(x, 0.0, ns, vars, gradState.data());

  if (vars[3] < 0.0) {
    std::cout << "Q[3] = " << vars[3] << " x = (" << x[0] << ", " << x[1]
              << std::endl;
  }
  return;
  std::cout << "\n";
  for (int i = 0; i < 4; ++i)
    std::cout << "Q[" << i << "] = " << vars[i] << std::endl;
}

void NavierStokes::Scenario::analyticalSolution(const double* const x,
                                                const double t,
                                                const NavierStokes& ns,
                                                Variables& vars,
                                                double* gradState) {
  throw std::logic_error("Analytical solution not implemented!");
}
void NavierStokes::Scenario::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t,
    const NavierStokes& ns, const double* const Q, double* S) {
  constexpr auto NumberOfVariables = DIMENSIONS + 2;  // TODO(Lukas) generalise?
  std::fill_n(S, NumberOfVariables, 0.0);
}

NavierStokes::BoundaryType NavierStokes::Scenario::getBoundaryType(int faceId) {
  return BoundaryType::wall;
}
