#include "NavierStokes.h"

NavierStokes::NavierStokes::NavierStokes() :
  NavierStokes(0.001, 21.0, 0) {
}
  
NavierStokes::NavierStokes::NavierStokes(double referenceT, double referenceViscosity, double sutherlandC) :
  referenceViscosity(referenceViscosity),
  referenceT(referenceT),
  sutherlandC(sutherlandC),
  sutherlandLambda(
		   (referenceViscosity * (referenceT + sutherlandC))/
		   std::pow(referenceT, 3./2)){
}



double NavierStokes::NavierStokes::evaluateTemperature(double rho, double pressure) const {
  return pressure/(gasConstant * rho);
}

double NavierStokes::NavierStokes::evaluateHeatConductionCoeff(double viscosity) const {
  // Use no heat conduction for now
  return 0.0;
  //return 1./Pr * viscosity * GAMMA * (1.0/(GAMMA - 1)) * gasConstant;
}

double NavierStokes::NavierStokes::evaluatePressure(double E, double rho, const tarch::la::Vector<3,double> &j) const {
  return (GAMMA-1) * (E - 0.5 * (1.0/rho) * j * j);
}

double NavierStokes::NavierStokes::evaluateViscosity(double T) const {
  // Use constant viscosity for now
  return referenceViscosity;

  // Sutherland's law
  //return sutherlandLambda * (std::pow(T, 3./2) / (T + sutherlandC));
}

void NavierStokes::NavierStokes::evaluateEigenvalues(const double* const Q, const int d, double* lambda) const {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j());
  assertion6(std::isfinite(p), p, vars.E(), vars.rho(), vars.j(0), vars.j(1), vars.j(2));

  const double u_n = vars.j(d)/vars.rho();
  const double temperature = evaluateTemperature(vars.rho(), p);
  //const double c = std::sqrt(GAMMA * gasConstant * temperature);
  const double c = std::sqrt(GAMMA * (p/vars.rho()));

  assertion3(std::isfinite(u_n), u_n, vars.j(d), vars.rho());
  assertion6(std::isfinite(temperature) && temperature >= 0.0, temperature, vars.rho(), vars.E(), vars.j() * vars.j(), p, u_n);
  assertion3(std::isfinite(c), c, u_n, temperature);

  eigs.rho() = u_n - c;
  eigs.E() = u_n + c;
  eigs.j(u_n, u_n, u_n);
}

void NavierStokes::NavierStokes::evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const {
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

void NavierStokes::NavierStokes::evaluateFlux(const double* Q, const double* gradQ, double** F) const {
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

  for (int dim = 0; dim < DIMENSIONS; ++dim) {
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(gradQ[idx_gradQ(dim, i)]), dim, i);
    }
  }

  // TODO: What if rho is tiny? Possibly add epsilon here for stability.
  //assert(vars.rho() > 10e-6);
  const auto invRho = 1. /vars.rho();
  assertion2(vars.rho() > 0, vars.rho(), invRho);
  assertion2(std::isfinite(invRho), vars.rho(), invRho);

  // Velocities
  const auto vx = invRho * vars.j(0);
  const auto vy = invRho * vars.j(1);
  const auto vz = invRho * vars.j(2);
  assertion2(std::isfinite(vx), invRho, vars.j(0));
  assertion2(std::isfinite(vy), invRho, vars.j(1));
  assertion2(std::isfinite(vz), invRho, vars.j(2));

  // Derivatives of velocities.
  const auto vx_dx = invRho * (gradQ[idx_gradQ(0, 1)] - vx * gradQ[(idx_gradQ(0,0))]);
  const auto vy_dx = invRho * (gradQ[idx_gradQ(0, 2)] - vy * gradQ[(idx_gradQ(0,0))]);
  const auto vz_dx = invRho * (gradQ[idx_gradQ(0, 3)] - vz * gradQ[(idx_gradQ(0,0))]);

  const auto vx_dy = invRho * (gradQ[idx_gradQ(1, 1)] - vx * gradQ[(idx_gradQ(1,0))]);
  const auto vy_dy = invRho * (gradQ[idx_gradQ(1, 2)] - vy * gradQ[(idx_gradQ(1,0))]);
  const auto vz_dy = invRho * (gradQ[idx_gradQ(1, 3)] - vz * gradQ[(idx_gradQ(1,0))]);

  const auto vx_dz = invRho * (gradQ[idx_gradQ(2, 1)] - vx * gradQ[(idx_gradQ(2,0))]);
  const auto vy_dz = invRho * (gradQ[idx_gradQ(2, 2)] - vy * gradQ[(idx_gradQ(2,0))]);
  const auto vz_dz = invRho * (gradQ[idx_gradQ(2, 3)] - vz * gradQ[(idx_gradQ(2,0))]);

  assertion3(std::isfinite(vx_dx), invRho, gradQ[idx_gradQ(0, 1)], gradQ[(idx_gradQ(0,0))]);
  assertion3(std::isfinite(vy_dx), invRho, gradQ[idx_gradQ(0, 2)], gradQ[(idx_gradQ(0,0))]);
  assertion3(std::isfinite(vz_dx), invRho, gradQ[idx_gradQ(0, 3)], gradQ[(idx_gradQ(0,0))]);
  assertion3(std::isfinite(vx_dy), invRho, gradQ[idx_gradQ(1, 1)], gradQ[(idx_gradQ(1,0))]);
  assertion3(std::isfinite(vy_dy), invRho, gradQ[idx_gradQ(1, 2)], gradQ[(idx_gradQ(1,0))]);
  assertion3(std::isfinite(vz_dy), invRho, gradQ[idx_gradQ(1, 3)], gradQ[(idx_gradQ(1,0))]);
  assertion3(std::isfinite(vx_dz), invRho, gradQ[idx_gradQ(2, 1)], gradQ[(idx_gradQ(2,0))]);
  assertion3(std::isfinite(vy_dz), invRho, gradQ[idx_gradQ(2, 2)], gradQ[(idx_gradQ(2,0))]);
  assertion3(std::isfinite(vz_dz), invRho, gradQ[idx_gradQ(2, 3)], gradQ[(idx_gradQ(2,0))]);

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
  const double Td = vars.rho();
  const double temperature = (1/gasConstant) * Tn/Td;
  assertion4(std::isfinite(temperature) && temperature >= 0.0, temperature, p, vars.rho(), vars.E());

  const double viscosity = evaluateViscosity(temperature); // TODO(Lukas) compute visc.
  const auto stressTensor = viscosity * ((2./3.) * div_v * I - 
					 (gradV + tarch::la::transpose(gradV)));
  //const auto stressTensor = 0.0 * I;

  // Compute derivatives of T with quotient rule:

  const double normJ = vars.j(0) * vars.j(0) + vars.j(1) * vars.j(1) + vars.j(2) * vars.j(2);

  // Derivative of inverse rho
  const double invRho_dx = -1 * gradQ[idx_gradQ(0,0)] / (invRho * invRho);
  const double invRho_dy = -1 * gradQ[idx_gradQ(1,0)] / (invRho * invRho);
  const double invRho_dz = -1 * gradQ[idx_gradQ(2,0)] / (invRho * invRho);
    
  // Derivative of norm of squared velocity densities times 0.5
  const double normJ_dx = gradQ[idx_gradQ(0,1)] + gradQ[idx_gradQ(0,2)] + gradQ[idx_gradQ(0,3)];
  const double normJ_dy = gradQ[idx_gradQ(1,1)] + gradQ[idx_gradQ(1,2)] + gradQ[idx_gradQ(1,3)];
  const double normJ_dz = gradQ[idx_gradQ(2,1)] + gradQ[idx_gradQ(2,2)] + gradQ[idx_gradQ(2,3)];
  
  // Derivatives of numerator of T
  const double scale_Tn = (GAMMA - 1.0);
  const double Tn_dx = scale_Tn * (gradQ[idx_gradQ(0,4)] - (invRho_dx * normJ + invRho * normJ_dx));
  const double Tn_dy = scale_Tn * (gradQ[idx_gradQ(1,4)] - (invRho_dy * normJ + invRho * normJ_dy));
  const double Tn_dz = scale_Tn * (gradQ[idx_gradQ(2,4)] - (invRho_dz * normJ + invRho * normJ_dz));

  // Derivatives of denominator of T
  const double Td_dx = gradQ[(idx_gradQ(0, 0))];
  const double Td_dy = gradQ[(idx_gradQ(1, 0))]; 
  const double Td_dz = gradQ[(idx_gradQ(2, 0))]; 

  // Assemble gradient of T
  const double scale = (1/gasConstant) *  1./(Td * Td);
  const double T_dx = scale * (Tn_dx*Td + Tn * Td_dx);
  const double T_dy = scale * (Tn_dy*Td + Tn * Td_dy);
  const double T_dz = scale * (Tn_dz*Td + Tn * Td_dz);
  assertion4(std::isfinite(T_dx), T_dx, Tn_dx, Td_dx, scale);
  assertion4(std::isfinite(T_dy), T_dy, Tn_dy, Td_dy, scale);
  assertion4(std::isfinite(T_dy), T_dz, Tn_dz, Td_dz, scale);
    
  
  //const tarch::la::Vector<3, double> gradT =  {T_dx, T_dy, T_dz};
  const tarch::la::Vector<3, double> gradT =  {0.0, 0.0, 0.0};
  const double kappa =  evaluateHeatConductionCoeff(viscosity);
  assertion2(std::isfinite(kappa), kappa, viscosity);
  

  // Full NS flux
  f.rho(vars.j());
  f.j(invRho * outerDot(vars.j(), vars.j()) + p*I + stressTensor);
  f.E(
      ((I * vars.E() + I * p + stressTensor) * (invRho * vars.j())) - kappa * gradT);
  /*
  // Only hyperbolic part, for debugging.
  f.rho ( vars.j() );
  f.j ( invRho *outerDot(vars.j(),vars.j()) + p*I );
  f.E ( invRho *(vars.E() + p) *vars.j() );
  */
 

  for (int i = 0; i < vars.variables(); ++i) {
    const auto cond = !std::isnan(Q[i]) && std::isfinite(Q[i]);
    if (!cond) {
      std::cout << i << std::endl;
      std::cout << "invRho = " << invRho << std::endl;
      for (int j = 0; j < vars.variables(); ++j) {
	std::cout << Q[j] << std::endl;
      }
    }
    assert(cond);
  }
}

