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
  NavierStokes() :
    NavierStokes(0.001, 21.0, 0) {
  }
  
  NavierStokes(double referenceT, double referenceViscosity, double sutherlandC) :
    referenceViscosity(referenceViscosity),
    referenceT(referenceT),
    sutherlandC(sutherlandC),
    sutherlandLambda(
		     (referenceViscosity * (referenceT + sutherlandC))/
		     std::pow(referenceT, 3./2)){
  }

  double evaluateTemperature(double rho, double pressure) const {
    return pressure/(gasConstant * rho);
  }
  
  double evaluateHeatConductionCoeff(double viscosity) const {
    // Use no heat conduction for now
    return 0.0;
    //return 1./Pr * viscosity * GAMMA * (1.0/(GAMMA - 1)) * gasConstant;
  }
  
  double evaluatePressure(double E, double rho, tarch::la::Vector<3,double> j) const {
    return (GAMMA-1) * (E - 0.5 * (1.0/rho) * j * j);
  }

  double evaluateViscosity(double T) const {
    // Use constant viscosity for now
    return referenceViscosity;
    
    // Sutherland's law
    //return sutherlandLambda * (std::pow(T, 3./2) / (T + sutherlandC));
  }

  void evaluateEigenvalues(const double* const Q, const int d, double* lambda) const {
    ReadOnlyVariables vars(Q);
    Variables eigs(lambda);

    const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j());

    const double u_n = vars.j(d)/vars.rho();
    const double temperature = evaluateTemperature(vars.rho(), p);
    double c = std::sqrt(GAMMA * gasConstant * temperature);

    eigs.rho() = u_n - c;
    eigs.E() = u_n + c;
    eigs.j(u_n, u_n, u_n);
  }

  void evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const {
    ReadOnlyVariables vars(Q);
    Variables eigs(lambda);

    const double pressure = evaluatePressure(vars.E(), vars.rho(), vars.j());
    const double T = evaluateTemperature(vars.rho(), pressure);
    const double viscosity = evaluateViscosity(T);

    // max_1 = 4/3 viscosity/rho
    // max_2 = (gamma viscosity)/(Pr * rho)
    // Only set two eigenvalues.
    // Results in correct max eigenvalues.
    
    // TODO(Lukas): Need to init to zero here?
    std::fill_n(lambda, vars.variables(), 0.0);
    lambda[0] = (4./3.) * (viscosity/vars.rho());
    lambda[1] = (GAMMA * viscosity) / (Pr * vars.rho());
  }

  void evaluateFlux(const double* Q, const double* gradQ, double** F) const {
    // Dimensions                        = 3
    // Number of variables + parameters  = 5 + 0
 
    ReadOnlyVariables vars(Q);
    Fluxes f(F);

    // Identity
    tarch::la::Matrix<3,3,double> I;
    I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;

    const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j());

    // gradQ: dim, variable
    kernels::idx2 idx_gradQ(DIMENSIONS, vars.variables());

    // Velocities
    const auto vx = vars.j(0) / vars.rho();
    const auto vy = vars.j(1) / vars.rho();
    const auto vz = vars.j(2) / vars.rho();

    // TODO: What if rho is tiny? Possibly add epsilon here for stability.
    //assert(vars.rho() > 10e-6);
    const auto invRho = 1. /vars.rho();

    // Derivatives of velocities.
    const auto vx_dx = invRho * (-gradQ[idx_gradQ(0, 1)] + vx * gradQ[(idx_gradQ(0,0))]);
    const auto vy_dx = invRho * (-gradQ[idx_gradQ(0, 2)] + vy * gradQ[(idx_gradQ(0,0))]);
    const auto vz_dx = invRho * (-gradQ[idx_gradQ(0, 3)] + vz * gradQ[(idx_gradQ(0,0))]);

    const auto vx_dy = invRho * (-gradQ[idx_gradQ(1, 1)] + vx * gradQ[(idx_gradQ(1,0))]);
    const auto vy_dy = invRho * (-gradQ[idx_gradQ(1, 2)] + vy * gradQ[(idx_gradQ(1,0))]);
    const auto vz_dy = invRho * (-gradQ[idx_gradQ(1, 3)] + vz * gradQ[(idx_gradQ(1,0))]);

    const auto vx_dz = invRho * (-gradQ[idx_gradQ(2, 1)] + vx * gradQ[(idx_gradQ(2,0))]);
    const auto vy_dz = invRho * (-gradQ[idx_gradQ(2, 2)] + vy * gradQ[(idx_gradQ(2,0))]);
    const auto vz_dz = invRho * (-gradQ[idx_gradQ(2, 3)] + vz * gradQ[(idx_gradQ(2,0))]);

    // Gradient and divergence 
    tarch::la::Matrix<3,3,double> gradV;
    gradV = vx_dx, vy_dx, vz_dx,
      vx_dy, vy_dy, vz_dy,
      vx_dz, vy_dz, vz_dz;
  
    const auto div_v = vx_dx + vy_dy + vz_dz;

    // Temperature: TODO, kappa currently set to 0
    // T = (p)/(vars.rho() * Pr)
    // Tn = numerator of T, Td = denominator of T
    const double Tn = p;
    const double Td = vars.rho() * gasConstant;
    const double temperature = Tn/Td;

    const double viscosity = evaluateViscosity(temperature); // TODO(Lukas) compute visc.
    const auto stressTensor = viscosity * ((2./3.) * div_v * I - 
      gradV + tarch::la::transpose(gradV));

    // Compute derivatives of T with quotient rule:
    const double normV = vx*vx + vy*vy + vz * vz;

    // Derivative of squared velocities.
    const double vx_vx_dx = 2 * vx * vx_dx;
    const double vy_vy_dx = 2 * vy * vy_dx;
    const double vz_vz_dx = 2 * vz * vz_dx;
    const double normV_dx = vx_vx_dx + vy_vy_dx + vz_vz_dx;

    const double vx_vx_dy = 2 * vx * vy_dx;
    const double vy_vy_dy = 2 * vy * vy_dx;
    const double vz_vz_dy = 2 * vz * vy_dx;
    const double normV_dy = vx_vx_dy + vy_vy_dy + vz_vz_dy;

    const double vx_vx_dz = 2 * vx * vx_dz;
    const double vy_vy_dz = 2 * vy * vy_dz;
    const double vz_vz_dz = 2 * vz * vz_dz;
    const double normV_dz = vx_vx_dz + vy_vy_dz + vz_vz_dz;
  
    // Derivatives of numerator of T
    const double scale_Tn = (GAMMA - 1.0);
    const double Tn_dx = scale_Tn * (gradQ[idx_gradQ(0,4)] - 0.5 * (gradQ[idx_gradQ(0,0)] * normV + vars.rho() *  normV_dx));
    const double Tn_dy = scale_Tn * (gradQ[idx_gradQ(1,4)] - 0.5 * (gradQ[idx_gradQ(1,0)] * normV + vars.rho() *  normV_dy));
    const double Tn_dz = scale_Tn * (gradQ[idx_gradQ(2,4)] - 0.5 * (gradQ[idx_gradQ(2,0)] * normV + vars.rho() *  normV_dz));

    // Derivatives of denominator of T
    const double Td_dx = Pr * gradQ[(idx_gradQ(0, 0))];
    const double Td_dy = Pr * gradQ[(idx_gradQ(1, 0))]; 
    const double Td_dz = Pr * gradQ[(idx_gradQ(2, 0))]; 

    // Assemble gradient of T
    const double scale = 1./(Td * Td);
    const double T_dx = scale * (Tn_dx*Td + Tn * Td_dx);
    const double T_dy = scale * (Tn_dy*Td + Tn * Td_dy);
    const double T_dz = scale * (Tn_dz*Td + Tn * Td_dz);
  
    const tarch::la::Vector<3, double> gradT =  {T_dx, T_dy, T_dz};
    const double kappa = 0.0;
  

    f.rho(vars.j());
    f.j(invRho * outerDot(vars.j(), vars.j()) + p*I + stressTensor);
    f.E(
	((I * vars.E() + I * p + stressTensor) * (invRho * vars.j())) - kappa * gradT);
  }

  static constexpr double GAMMA = 1.4;
  static constexpr double gasConstant = 8.3144598; // TODO
  static constexpr double Pr = 0.71; // Prandtl number, for air
  double referenceViscosity;
  double referenceT;
  double sutherlandC;
  double sutherlandLambda;
  
};
