#include "ABCFlow.h"

void NavierStokes::ABCFlow::analyticalSolution(const double *const x,
                                               const double t,
                                               const NavierStokes::PDE &ns,
                                               Variables &vars,
                                               double *gradState) {
  // TODO(Lukas) Check everything here again!
  assert(DIMENSIONS == 3);

  // Variable shortcuts
  const auto rho = NavierStokesSolver_ADERDG_Variables::shortcuts::rho;
  const auto j = NavierStokesSolver_ADERDG_Variables::shortcuts::j;
  const auto E = NavierStokesSolver_ADERDG_Variables::shortcuts::E;
  const auto Z = E + 1;  // Only defined if coupling is used!

  const auto c = 100 / ns.gamma;  // TODO(Lukas) Find good c!
  const auto Ft = std::exp(-2 * ns.referenceViscosity * t);
  const auto pressure =
      (std::sin(x[1]) * std::cos(x[0]) + std::sin(x[0]) * std::cos(x[2]) +
       std::sin(x[2]) * std::cos(x[1])) *
          Ft +
      c;

  vars[rho] = 1.0;  // Solution for incompressible NS.
  vars[j + 0] = (std::sin(x[2]) + std::cos(x[1])) * Ft;
  vars[j + 1] = (std::sin(x[0]) + std::cos(x[2])) * Ft;
  vars[j + 2] = (std::sin(x[1]) + std::cos(x[0])) * Ft;
  vars[E] = ns.evaluateEnergy(vars[rho], pressure, vars.j());

  auto idxGradQ = kernels::idx2(DIMENSIONS, vars.SizeVariables);

  gradState[idxGradQ(0, rho)] = 0;
  gradState[idxGradQ(0, j + 0)] = 0;
  gradState[idxGradQ(0, j + 1)] = Ft * std::cos(x[0]);
  gradState[idxGradQ(0, j + 2)] = -Ft * std::sin(x[0]);
  gradState[idxGradQ(0, E)] =
      1.0 * (std::sin(x[0]) + std::cos(x[2])) * (Ft * Ft) * std::cos(x[0]) -
      1.0 * (std::sin(x[1]) + std::cos(x[0])) * (Ft * Ft) * std::sin(x[0]) +
      (std::sin(x[0]) * std::sin(x[1]) - std::cos(x[0]) * std::cos(x[2])) * Ft /
          (-ns.gamma + 1);

  gradState[idxGradQ(1, rho)] = 0;
  gradState[idxGradQ(1, j + 0)] = -Ft * std::sin(x[1]);
  gradState[idxGradQ(1, j + 1)] = 0;
  gradState[idxGradQ(1, j + 2)] = Ft * std::cos(x[1]);
  gradState[idxGradQ(1, E)] =
      1.0 * (std::sin(x[1]) + std::cos(x[0])) * (Ft * Ft) * std::cos(x[1]) -
      1.0 * (std::sin(x[2]) + std::cos(x[1])) * (Ft * Ft) * std::sin(x[1]) +
      (std::sin(x[1]) * std::sin(x[2]) - std::cos(x[0]) * std::cos(x[1])) * Ft /
          (-ns.gamma + 1);

  gradState[idxGradQ(2, rho)] = 0;
  gradState[idxGradQ(2, j + 0)] = Ft * std::cos(x[2]);
  gradState[idxGradQ(2, j + 1)] = -Ft * std::sin(x[2]);
  gradState[idxGradQ(2, j + 2)] = 0;
  gradState[idxGradQ(2, E)] =
      -1.0 * (std::sin(x[0]) + std::cos(x[2])) * (Ft * Ft) * std::sin(x[2]) +
      1.0 * (std::sin(x[2]) + std::cos(x[1])) * (Ft * Ft) * std::cos(x[2]) +
      (std::sin(x[0]) * std::sin(x[2]) - std::cos(x[1]) * std::cos(x[2])) * Ft /
          (-ns.gamma + 1);
}

NavierStokes::BoundaryType NavierStokes::ABCFlow::getBoundaryType(int faceId) {
  return BoundaryType::analytical;
}