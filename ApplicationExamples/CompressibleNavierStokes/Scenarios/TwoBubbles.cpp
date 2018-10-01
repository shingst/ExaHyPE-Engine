//
// Created by lukas on 28/08/18.
//

#include "TwoBubbles.h"
#include <array>

NavierStokes::TwoBubbles::Bubble::Bubble(const double tempDifference,
                                         const double size, const double decay,
                                         const double centerX,
                                         const double centerZ)
    : tempDifference(tempDifference),
      size(size),
      decay(decay),
      centerX(centerX),
      centerZ(centerZ) {}

void NavierStokes::TwoBubbles::initialValues(const double* const x,
                                             const PDE& ns, Variables& vars) {
  // For details see:
  // Robert (1993),
  // https://doi.org/10.1175/1520-0469(1993)050<1865:BCEWAS>2.0.CO;2
  const auto posX = x[0];
  const auto posZ = (DIMENSIONS == 2) ? x[1] : x[2];

  const double backgroundPressure =
      ns.referencePressure;  // [Pa], background pressure of atmosphere at sea
                             // level
  const double backgroundT = 300;  // [K], background potential temperature

  // Parameters for bubbles.
  // Large, hot bubble
  const double A_1 = 0.5;  // temp. difference
  const double a_1 = 150;  // [m]
  const double S_1 = 50;   // [m]
  const double x_1 = 500;  // [m]
  const double z_1 = 300;  // [m]

  // Small, cold bubble
  const double A_2 = -0.15;  // temp. difference
  const double a_2 = 0;      // [m]
  const double S_2 = 50;     // [m]
  const double x_2 = 560;    // [m]
  const double z_2 = 640;    // [m]

  const auto bubbles = std::array<Bubble, 2>{
      Bubble(A_1, a_1, S_1, x_1, z_1),
      Bubble(A_2, a_2, S_2, x_2, z_2),
  };

  double potentialT = backgroundT;

  for (const auto& bubble : bubbles) {
    // Check if we are in range of the bubble.
    const auto distX = posX - bubble.centerX;
    const auto distZ = posZ - bubble.centerZ;
    const auto distanceToCenter = std::sqrt(distX * distX + distZ * distZ);
    // Inside the center of the bubble apply tempDifference without decay,
    if (distanceToCenter <= bubble.size) {
      potentialT += bubble.tempDifference;
    } else {
      // Exponentially decay temperature outside of size.
      const auto d = distanceToCenter - bubble.size;
      potentialT += bubble.tempDifference *
                    std::exp(-(d * d) / (bubble.decay * bubble.decay));
    }
    break; // Only 1 bubble for now.
  }

  const double g = 9.81;  // [m/s^2]

  // Air is initially at rest.
#if DIMENSIONS == 2
  vars.j(0, 0);
#elif DIMENSIONS == 3
  vars.j(0, 0, 0);
#endif

  const double pressure = std::pow((ns.gasConstant*ns.gamma*backgroundT*std::pow(std::pow(ns.gamma, ns.gamma/(ns.gamma - 1))*ns.referencePressure, (ns.gamma - 1)/ns.gamma) - ns.gasConstant*backgroundT*std::pow(std::pow(ns.gamma, ns.gamma/(ns.gamma - 1))*ns.referencePressure, (ns.gamma - 1)/ns.gamma) - g*ns.gamma*std::pow(ns.referencePressure, ns.gasConstant/c_p)*posZ*(ns.gamma - 1) + g*std::pow(ns.referencePressure, ns.gasConstant/c_p)*posZ*(ns.gamma - 1))/(ns.gasConstant*ns.gamma*backgroundT*(ns.gamma - 1)), ns.gamma/(ns.gamma - 1));
  const double poTToT =
          std::pow((pressure / backgroundPressure), ns.gasConstant / ns.c_p);
  const double temperature = potentialT * poTToT;
  vars.rho() = pressure / (ns.gasConstant * temperature);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

  /*
   const double exner = 1 - (g * posZ)/(ns.c_p * backgroundT);
   const double pressure = ns.referencePressure * std::pow(
           exner,
           1/(ns.c_p/ns.gasConstant)
           );
  const double poTToT =
          std::pow((pressure / backgroundPressure), ns.gasConstant / ns.c_p);
  const double temperature = potentialT * poTToT;
  vars.rho() = pressure / (ns.gasConstant * temperature);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
   */
}

void NavierStokes::TwoBubbles::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);
  const double g = -9.81;
  S[DIMENSIONS] = Q[0] * g;
  S[DIMENSIONS + 1] = Q[2] * g;
}

const double NavierStokes::TwoBubbles::getGamma() const { return gamma; }

const double NavierStokes::TwoBubbles::getPr() const { return Pr; }

const double NavierStokes::TwoBubbles::getC_v() const { return c_v; }

const double NavierStokes::TwoBubbles::getC_p() const { return c_p; }

const double NavierStokes::TwoBubbles::getGasConstant() const {
  return gasConstant;
}

const double NavierStokes::TwoBubbles::getReferencePressure() const {
  return referencePressure;
}
