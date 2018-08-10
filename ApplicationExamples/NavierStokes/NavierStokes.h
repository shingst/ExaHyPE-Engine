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

  double evaluateTemperature(double rho, double pressure) const;
  double evaluateHeatConductionCoeff(double viscosity) const;  
  double evaluatePressure(double E, double rho, const tarch::la::Vector<3,double> &j) const;
  double evaluateViscosity(double T) const;
  void evaluateEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateFlux(const double* Q, const double* gradQ, double** F) const;

  static constexpr double GAMMA = 1.4;
  static constexpr double gasConstant = 8.3144598; // TODO
  static constexpr double Pr = 0.71; // Prandtl number, for air
  double referenceViscosity;
  double referenceT;
  double sutherlandC;
  double sutherlandLambda;
  
};
