#include "TaylorGreen.h"

void NavierStokes::TaylorGreen::initialValues(const double* const x,
                                              const NavierStokes& ns,
                                              Variables& vars) {
  assert(DIMENSIONS == 2);
  // 2D-Scenario
  vars.rho() = 1.0;
  vars.j(0) = 1 * std::cos(x[0]) * std::sin(x[1]);
  vars.j(1) = -1 * std::sin(x[0]) * std::cos(x[1]);
#if DIMENSIONS == 3
  vars.j(2) = 0;
#endif
  const double pressure =
      -1 * (vars.rho() / 4) * (std::cos(2 * x[0]) + std::cos(2 * x[1]));

  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
}
void NavierStokes::TaylorGreen::analyticalSolution(const double* const x,
                                                   const double t,
                                                   const NavierStokes& ns,
                                                   Variables& vars,
                                                   double* gradState) {
  kernels::idx2 idxGradQ(DIMENSIONS, vars.SizeVariables);

  const double Ft = std::exp(-2 * ns.referenceViscosity * t);

  vars.rho() = 1.0;
  vars.j(0) = 1 * std::cos(x[0]) * std::sin(x[1]) * Ft;
  vars.j(1) = -1 * std::sin(x[0]) * std::cos(x[1]) * Ft;
#if DIMENSIONS == 3
  vars.j(2) = 0;
#endif
  const auto pressure = -1 * (vars.rho() / 4) *
                        (std::cos(2 * x[0]) + std::cos(2 * x[1])) * Ft * Ft;
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

  // Assuming rho is constant.
  // j(0)
  gradState[idxGradQ(0, 1)] = -1 * std::sin(x[0]) * std::sin(x[1]) * Ft;
  gradState[idxGradQ(1, 1)] = 1 * std::cos(x[0]) * std::cos(x[1]) * Ft;

  // j(1)
  gradState[idxGradQ(0, 2)] = -1 * std::cos(x[0]) * std::cos(x[1]) * Ft;
  gradState[idxGradQ(1, 2)] = 1 * std::sin(x[0]) * std::sin(x[1]) * Ft;

  // j(2) is zero for 3d and non-existant for 2d

  // E, idx 3 for 2d and idx 4 for 3d
  // Assume that rho is 1 at boundary and constant (=zero derivative)
  // TODO(Lukas) Fix these! They are wrong!
  constexpr double e_idx = DIMENSIONS + 1;
  const double e_factor = 0.25 * Ft * Ft;
  gradState[idxGradQ(0, e_idx)] =
      e_factor *
      (1 * std::sin(2 * x[0] - 2 * x[1]) + std::sin(2 * x[0] + 2 * x[1]));
  gradState[idxGradQ(1, e_idx)] =
      e_factor *
      (-1 * std::sin(2 * x[0] - 2 * x[1]) + std::sin(2 * x[0] + 2 * x[1]));
  // gradState[idxGradQ(2, e_idx)] = 0.0;
}

NavierStokes::BoundaryType NavierStokes::TaylorGreen::getBoundaryType(
    int faceId) {
  return BoundaryType::analytical;
}