#ifndef NAVIERSTOKES_NAVIERSTOKES_H
#define NAVIERSTOKES_NAVIERSTOKES_H

#include <cmath>
#include <memory>
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
  NavierStokes(double referenceViscosity, double referencePressure, double gamma, double Pr, double c_v,
          double c_p, double gasConstant);

  double evaluateEnergy(double rho, double pressure, const tarch::la::Vector<DIMENSIONS,double> &j) const;
  double evaluateTemperature(double rho, double pressure) const;
  double evaluateHeatConductionCoeff(double viscosity) const;  
  double evaluatePressure(double E, double rho, const tarch::la::Vector<DIMENSIONS,double> &j) const;
  double evaluateViscosity(double T) const;
  void evaluateEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateFlux(const double* Q, const double* gradQ, double** F, bool temperatureDiff=true) const;


  double Pr; // [1] Prandtl number
  double referencePressure; // [Pa], reference pressure used by some scenarios.

  double gamma; // [1] Ratio of specific heats, c_p/c_v

  // [m^2/(s^2 K)] = [J/(kg K)]
  double c_v; // [J/(kg K)] Specific heat capacity at constant volume
  double c_p; // [J/(kg K)] Specific heat capacity at constant pressure
  double gasConstant; // [J/(kg K)] Specific gas constant

  double referenceViscosity; // [Pa s] = [kg/(ms)] Reference dynamic viscosity


  // Unused:
  double referenceT;
  double sutherlandC;
  double sutherlandLambda;
};

#endif // NAVIERSTOKES_NAVIERSTOKES_H
