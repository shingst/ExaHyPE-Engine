#ifndef NAVIERSTOKES_NAVIERSTOKES_H
#define NAVIERSTOKES_NAVIERSTOKES_H

#include <cmath>
#include "tarch/la/Vector.h"
#include "EulerSolver_Variables.h"
#include "kernels/KernelUtils.h"

namespace NavierStokes {
  class NavierStokes;
}

using Variables = Euler::AbstractEulerSolver::Variables;
using ReadOnlyVariables = Euler::AbstractEulerSolver::ReadOnlyVariables;
using Fluxes = Euler::AbstractEulerSolver::Fluxes;

class NavierStokes::NavierStokes {
public:
  NavierStokes();
  NavierStokes(double referenceT, double referenceViscosity, double sutherlandC);

  double evaluateEnergy(double rho, double pressure, const tarch::la::Vector<DIMENSIONS,double> &j) const;
  double evaluateTemperature(double rho, double pressure) const;
  double evaluateHeatConductionCoeff(double viscosity) const;  
  double evaluatePressure(double E, double rho, const tarch::la::Vector<DIMENSIONS,double> &j) const;
  double evaluateViscosity(double T) const;
  void evaluateEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateFlux(const double* Q, const double* gradQ, double** F) const;

  // All constants are for air.
  // TODO(Lukas) Make sure all constants are for air at same temperature.
  // TODO(Lukas) Maybe refactor those s.t. eqs can be used for different fluids?
  static constexpr double GAMMA = 1.4; // [1] Ratio of specific heats
  static constexpr double c_p = 1.005 * 1000; // [J/(kg K)] Specific heat capacity at constant pressure
  static constexpr double referencePressure = 10000; // Pa, reference pressure used by some scenarios.
  static constexpr double gasConstant = 287.058; // [m^2/(s^2 K)] = [J/(kg K)] Specific gas constant, for dry air
  static constexpr double Pr = 0.71; // [1] Prandtl number
  double referenceViscosity;
  double referenceT;
  double sutherlandC;
  double sutherlandLambda;
  
};

#endif // NAVIERSTOKES_NAVIERSTOKES_H
