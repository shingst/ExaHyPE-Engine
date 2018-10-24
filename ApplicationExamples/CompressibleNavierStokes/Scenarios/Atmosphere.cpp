#include "Atmosphere.h"

double NavierStokes::computeHydrostaticPressure(const PDE& ns, double g,
                                                double posZ,
                                                double backgroundT) {
  return std::pow(
      (ns.gasConstant * ns.gamma * backgroundT *
           std::pow(std::pow(ns.gamma, ns.gamma / (ns.gamma - 1)) *
                        ns.referencePressure,
                    (ns.gamma - 1) / ns.gamma) -
       ns.gasConstant * backgroundT *
           std::pow(std::pow(ns.gamma, ns.gamma / (ns.gamma - 1)) *
                        ns.referencePressure,
                    (ns.gamma - 1) / ns.gamma) -
       g * ns.gamma * std::pow(ns.referencePressure, ns.gasConstant / ns.c_p) *
           posZ * (ns.gamma - 1) +
       g * std::pow(ns.referencePressure, ns.gasConstant / ns.c_p) * posZ *
           (ns.gamma - 1)) /
          (ns.gasConstant * ns.gamma * backgroundT * (ns.gamma - 1)),
      ns.gamma / (ns.gamma - 1));
}

double NavierStokes::potentialTToT(const PDE& ns, double pressure,
                                   double potentialT) {
  const double poTToT =
      std::pow((pressure / ns.referencePressure), ns.gasConstant / ns.c_p);
  return potentialT * poTToT;
}