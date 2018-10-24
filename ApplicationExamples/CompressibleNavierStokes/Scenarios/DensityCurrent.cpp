#include <array>

#include "Atmosphere.h"
#include "DensityCurrent.h"

void NavierStokes::DensityCurrent::initialValues(const double* const x,
                                                 const PDE& ns,
                                                 Variables& vars) {
  // For details see:
  // TODO(Lukas) Add reference to correct paper!
  const auto posX = x[0];
  const auto posZ = (DIMENSIONS == 2) ? x[1] : x[2];

  const double backgroundPressure =
      ns.referencePressure;  // [Pa], background pressure of atmosphere at sea
                             // level
  const double backgroundT = 300;  // [K], background potential temperature

  double potentialT = backgroundT;

  const auto centerX = 0;
  const auto centerZ = 3000;
  const auto sizeX = 4000;
  const auto sizeZ = 2000;

  const auto distX = (posX - centerX) / sizeX;
  const auto distZ = (posZ - centerZ) / sizeZ;
  const auto r = std::sqrt(distX * distX + distZ * distZ);
  const auto r_c = 1.0;

  if (r <= r_c) {
    const auto pertubationSize = -15;
    const auto pi = std::acos(-1);
    potentialT += pertubationSize / 2 * (1 + std::cos(pi * r));
  }

  const double g = 9.81;  // [m/s^2]

  // Air is initially at rest.
#if DIMENSIONS == 2
  vars.j(0, 0);
#elif DIMENSIONS == 3
  vars.j(0, 0, 0);
#endif

  const auto pressure = computeHydrostaticPressure(ns, g, posZ, backgroundT);
  const auto temperature = potentialTToT(ns, pressure, potentialT);
  vars.rho() = pressure / (ns.gasConstant * temperature);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
}

void NavierStokes::DensityCurrent::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);
  const double g = -9.81;
  S[DIMENSIONS] = Q[0] * g;
  S[DIMENSIONS + 1] = Q[2] * g;
}

const double NavierStokes::DensityCurrent::getGamma() const { return gamma; }

const double NavierStokes::DensityCurrent::getPr() const { return Pr; }

const double NavierStokes::DensityCurrent::getC_v() const { return c_v; }

const double NavierStokes::DensityCurrent::getC_p() const { return c_p; }

const double NavierStokes::DensityCurrent::getGasConstant() const {
  return gasConstant;
}

const double NavierStokes::DensityCurrent::getReferencePressure() const {
  return referencePressure;
}