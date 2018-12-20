//
#include <array>

#include "Atmosphere.h"
#include "TwoBubbles.h"

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
                                             const PDE& ns,
                                             Variables& vars) {
  initialValues(x, ns, vars, 0.0);
}

void NavierStokes::TwoBubbles::initialValues(const double* const x,
                                             const PDE& ns, Variables& vars,
                                             double initialZ) {
  // For details see:
  // Robert (1993),
  // https://doi.org/10.1175/1520-0469(1993)050<1865:BCEWAS>2.0.CO;2
  const auto posX = x[0];
  const auto posZ = (DIMENSIONS == 2) ? x[1] : x[2];

  const double backgroundPotentialT = getBackgroundPotentialTemperature();  // [K], background potential temperature

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

  double potentialT = backgroundPotentialT;

  double Z = 0.0; // TODO(Lukas) Only used for coupling test!

  for (const auto& bubble : bubbles) {
    // Check if we are in range of the bubble.
    const auto distX = posX - bubble.centerX;
    const auto distZ = posZ - bubble.centerZ;
    const auto distanceToCenter = std::sqrt(distX * distX + distZ * distZ);
    // Inside the center of the bubble apply tempDifference without decay,
    if (distanceToCenter <= bubble.size) {
      potentialT += bubble.tempDifference;
      Z = initialZ;
    } else {
      // Exponentially decay temperature outside of size.
      const auto d = distanceToCenter - bubble.size;
      potentialT += bubble.tempDifference *
                    std::exp(-(d * d) / (bubble.decay * bubble.decay));
      Z = initialZ * std::exp(-(d * d) / (bubble.decay * bubble.decay));
    }
    break; // TODO(Lukas) Add second bubble.
  }

  const double g = 9.81;  // [m/s^2]

  // Air is initially at rest.
#if DIMENSIONS == 2
  vars.j(0, 0);
#elif DIMENSIONS == 3
  vars.j(0, 0, 0);
#endif

  // First compute overall state
  const auto pressure = computeHydrostaticPressure(ns, getGravity(), posZ, backgroundPotentialT);
  const auto temperature = potentialTToT(ns, pressure, potentialT);
  vars.rho() = pressure / (ns.gasConstant * temperature);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j(), Z, ns.getHeight(vars.data()));
  assertion4(std::isfinite(vars.rho()), vars.rho(), vars.E(), temperature, pressure);
  assertion4(std::isfinite(vars.E()), vars.rho(), vars.E(), temperature, pressure);
  ns.setZ(vars.data(), Z);

  if (ns.useBackgroundState) {
    // TODO(Lukas) Refactor?
    // Then compute background state (without pot.T. pertubation
    const auto backgroundTemperature = potentialTToT(ns, pressure, backgroundPotentialT);
    const auto backgroundRho = pressure / (ns.gasConstant * backgroundTemperature);
    ns.setBackgroundState(vars.data(), backgroundRho, pressure);
  }
}

void NavierStokes::TwoBubbles::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);

  double rhoPertubation = Q[0];
  if (ns.useBackgroundState) {
    auto backgroundRho = 0.0;
    auto backgroundPressure = 0.0;
    std::tie(backgroundRho, backgroundPressure) = ns.getBackgroundState(Q);

    rhoPertubation -= backgroundRho;
  }
  S[DIMENSIONS] = -1 * rhoPertubation * getGravity();

  // Only use this source term if the gravitational force is not already
  // included in the pressure.
  if (!ns.useGravity) {
    S[DIMENSIONS + 1] = -1 * Q[2] * getGravity();
  }
}

double NavierStokes::TwoBubbles::getGamma() const { return gamma; }

double NavierStokes::TwoBubbles::getPr() const { return Pr; }

double NavierStokes::TwoBubbles::getC_v() const { return c_v; }

double NavierStokes::TwoBubbles::getC_p() const { return c_p; }

double NavierStokes::TwoBubbles::getGasConstant() const {
  return gasConstant;
}

double NavierStokes::TwoBubbles::getReferencePressure() const {
  return referencePressure;
}

double NavierStokes::TwoBubbles::getGravity() const {
  return 9.81;
}

double NavierStokes::TwoBubbles::getBackgroundPotentialTemperature() const {
  return 300;
}

NavierStokes::BoundaryType NavierStokes::TwoBubbles::getBoundaryType(int faceId) {
  // For boundaries in Z direction, we need to reconstruct
  // the temperature flux coming from the boundary.
#if DIMENSIONS == 2
  const bool isZFace = faceId == 2 || faceId == 3;
#else
  const bool isZFace = faceId == 4 || faceId == 5;
#endif
  if (isZFace) {
    return BoundaryType::hydrostaticWall;
  }
  return BoundaryType::hydrostaticWall;
  // TODO(Lukas) Reconsider!
  return BoundaryType::wall;

}

