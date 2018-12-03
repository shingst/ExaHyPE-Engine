//
// Created by lukas on 01/12/18.
//

#include "Detonation.h"

void NavierStokes::Detonation::initialValues(const double* const x,
                                             const PDE& ns, Variables& vars) {
  // B: Burned, U: Unburned
  const auto rhoB = 1.4;
  const auto rhoU = 0.887565;
  const auto pB = 1.0;
  const auto pU = 0.191709;
  const auto ZB = 0.0;
  const auto ZU = 1.0;

  const auto alpha = std::atan2(x[1], x[0]); // Angle in polar coordinates
  const auto uB = 0.0;
  const auto uU = -0.577350 * std::cos(alpha);
  const auto vB = 0.0;
  const auto vU = -0.577350 * std::sin(alpha);

  const auto centerX = 0.5;
  const auto centerY = 0.5;
  const auto distX = x[0] - centerX;
  const auto distY = x[1] - centerY;
  const auto distToCenter = std::sqrt(
          distX * distX + distY * distY
          );
  const auto radiusBurned = 0.3;
  const auto isBurned = distToCenter < radiusBurned;

  double pressure = -1;
  if (isBurned) {
     vars.rho() = rhoB;
     pressure = pB;
     vars.j(uB * rhoB, vB * rhoB);
     ns.setZ(vars.data(), ZB);
  } else {
     vars.rho() = rhoU;
     pressure = pU;
     vars.j(uU * rhoU, vU * rhoU);
     ns.setZ(vars.data(), ZU);
  }

  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j(), ns.getZ(vars.data()));

  assert(ns.getZ(vars.data()) >= 0 && ns.getZ(vars.data()) <= 1.0);
}
double NavierStokes::Detonation::getMolecularDiffusionCoeff() const {
  return 0.0; // Euler, for now.
}

double NavierStokes::Detonation::getQ0() const {
  return 25.0;
}
void NavierStokes::Detonation::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);
  const auto vars = ReadOnlyVariables{Q};

  const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), Q[vars.size()]);
  const auto T = ns.evaluateTemperature(vars.rho(), pressure);

  const auto reactionTimeScale = 0.1; // TODO(Lukas) Try stiff ts?
  const auto ignitionT = 0.26;
  if (T > ignitionT) {
    // Critical temperature, start decay.
    ns.setZ(S, -1.0/reactionTimeScale * ns.getZ(Q));
  } else {
    ns.setZ(S, 0.0);
  }
}

bool NavierStokes::Detonation::getUseAdvection() const {
  return true;
}