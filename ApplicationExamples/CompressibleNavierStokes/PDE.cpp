#include "PDE.h"

NavierStokes::PDE::PDE() :
  PDE(0.1, 10000, 1.4, 0.7, 1, 1.4, 0.4) {
}

NavierStokes::PDE::PDE(double referenceViscosity, double referencePressure, double gamma, double Pr,
        double c_v, double c_p, double gasConstant) :
  referenceViscosity(referenceViscosity),
  referencePressure(referencePressure),
  gamma(gamma),
  Pr(Pr),
  c_v(c_v),
  c_p(c_p),
  gasConstant(gasConstant) {
}


double NavierStokes::PDE::evaluateEnergy(double rho, double pressure, const tarch::la::Vector<DIMENSIONS,double> &j) const {
  const auto invRho = 1./rho;
  return pressure/(gamma - 1) + 0.5 * (invRho * j * j);
}

double NavierStokes::PDE::evaluateTemperature(double rho, double pressure) const {
  return pressure/(gasConstant * rho);
}

double NavierStokes::PDE::evaluateHeatConductionCoeff(double viscosity) const {
  return 1./Pr * viscosity * gamma * c_v;
}

double NavierStokes::PDE::evaluatePressure(double E, double rho, const tarch::la::Vector<DIMENSIONS,double> &j) const {
  return (gamma-1) * (E - 0.5 * (1.0/rho) * j * j);
}

double NavierStokes::PDE::evaluateViscosity(double T) const {
  return referenceViscosity;
}

void NavierStokes::PDE::evaluateEigenvalues(const double* const Q, const int d, double* lambda) const {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j());
  assertion4(std::isfinite(p), p, vars.E(), vars.rho(), vars.j());

  const double u_n = vars.j(d)/vars.rho();
  const double temperature = evaluateTemperature(vars.rho(), p);
  //const double c = std::sqrt(gamma * gasConstant * temperature);
  const double c = std::sqrt(gamma * (p/vars.rho()));

  assertion3(std::isfinite(u_n), u_n, vars.j(d), vars.rho());
  assertion6(std::isfinite(temperature) && temperature >= 0.0, temperature, vars.rho(), vars.E(), vars.j() * vars.j(), p, u_n);
  assertion3(std::isfinite(c), c, u_n, temperature);

  std::fill_n(lambda, vars.variables(), u_n);
  lambda[0] -= c;
  lambda[1] += c;
}

void NavierStokes::PDE::evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const {
  ReadOnlyVariables vars(Q);

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
  lambda[1] = (gamma * viscosity) / (Pr * vars.rho());
}

void NavierStokes::PDE::evaluateFlux(const double* Q, const double* gradQ, double** F, bool temperatureDiffusion) const {
  ReadOnlyVariables vars(Q);
  Fluxes f(F);

  // Identity
#if DIMENSIONS == 2
  tarch::la::Matrix<2,2,double> I;
  I = 1, 0,
      0, 1;
#elif DIMENSIONS == 3
  tarch::la::Matrix<3,3,double> I;
  I = 1, 0, 0,
      0, 1, 0,
      0, 0, 1;
#endif

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
  const auto invRho = 1. / (vars.rho());
  const auto invRho2 = 1./ (vars.rho() * vars.rho());
  const auto invRho3 = 1./ (vars.rho() * vars.rho() * vars.rho());

  assertion2(vars.rho() > 0, vars.rho(), invRho);
  assertion2(std::isfinite(invRho), vars.rho(), invRho);

  // Velocities
  const auto vx = invRho * vars.j(0);
  const auto vy = invRho * vars.j(1);
#if DIMENSIONS == 3
  const auto vz = invRho * vars.j(2);
#endif

  assertion2(std::isfinite(vx), invRho, vars.j(0));
  assertion2(std::isfinite(vy), invRho, vars.j(1));
#if DIMENSIONS == 3
  assertion2(std::isfinite(vz), invRho, vars.j(2));
#endif

  // Derivatives of velocities.
  const auto vx_dx = invRho * (gradQ[idx_gradQ(0, 1)] - vx * gradQ[(idx_gradQ(0,0))]);
  const auto vy_dx = invRho * (gradQ[idx_gradQ(0, 2)] - vy * gradQ[(idx_gradQ(0,0))]);

  const auto vx_dy = invRho * (gradQ[idx_gradQ(1, 1)] - vx * gradQ[(idx_gradQ(1,0))]);
  const auto vy_dy = invRho * (gradQ[idx_gradQ(1, 2)] - vy * gradQ[(idx_gradQ(1,0))]);
#if DIMENSIONS == 3
  const auto vz_dx = invRho * (gradQ[idx_gradQ(0, 3)] - vz * gradQ[(idx_gradQ(0,0))]);
  const auto vz_dy = invRho * (gradQ[idx_gradQ(1, 3)] - vz * gradQ[(idx_gradQ(1,0))]);

  const auto vx_dz = invRho * (gradQ[idx_gradQ(2, 1)] - vx * gradQ[(idx_gradQ(2,0))]);
  const auto vy_dz = invRho * (gradQ[idx_gradQ(2, 2)] - vy * gradQ[(idx_gradQ(2,0))]);
  const auto vz_dz = invRho * (gradQ[idx_gradQ(2, 3)] - vz * gradQ[(idx_gradQ(2,0))]);
#endif

  assertion3(std::isfinite(vx_dx), invRho, gradQ[idx_gradQ(0, 1)], gradQ[(idx_gradQ(0,0))]);
  assertion3(std::isfinite(vy_dx), invRho, gradQ[idx_gradQ(0, 2)], gradQ[(idx_gradQ(0,0))]);
  assertion3(std::isfinite(vx_dy), invRho, gradQ[idx_gradQ(1, 1)], gradQ[(idx_gradQ(1,0))]);
  assertion3(std::isfinite(vy_dy), invRho, gradQ[idx_gradQ(1, 2)], gradQ[(idx_gradQ(1,0))]);
#if DIMENSIONS == 3
  assertion3(std::isfinite(vz_dx), invRho, gradQ[idx_gradQ(0, 3)], gradQ[(idx_gradQ(0,0))]);
  assertion3(std::isfinite(vz_dy), invRho, gradQ[idx_gradQ(1, 3)], gradQ[(idx_gradQ(1,0))]);
  assertion3(std::isfinite(vx_dz), invRho, gradQ[idx_gradQ(2, 1)], gradQ[(idx_gradQ(2,0))]);
  assertion3(std::isfinite(vy_dz), invRho, gradQ[idx_gradQ(2, 2)], gradQ[(idx_gradQ(2,0))]);
  assertion3(std::isfinite(vz_dz), invRho, gradQ[idx_gradQ(2, 3)], gradQ[(idx_gradQ(2,0))]);
#endif

  // Gradient and divergence
#if DIMENSIONS == 2
  tarch::la::Matrix<2,2,double> gradV;
  gradV = vx_dx, vy_dx,
    vx_dy, vy_dy;

  const auto div_v = vx_dx + vy_dy;
#elif DIMENSIONS == 3
  tarch::la::Matrix<3,3,double> gradV;
  gradV = vx_dx, vy_dx, vz_dx,
    vx_dy, vy_dy, vz_dy,
    vx_dz, vy_dz, vz_dz;

  const auto div_v = vx_dx + vy_dy + vz_dz;
#endif

  // T = (p)/(vars.rho() * Pr)
  // Tn = numerator of T, Td = denominator of T
  const double Tn = p;
  const double Td = vars.rho();
  const double temperature = (1/gasConstant) * Tn/Td;
  assertion4(std::isfinite(temperature) && temperature >= 0.0, temperature, p, vars.rho(), vars.E());

  const double viscosity = evaluateViscosity(temperature); // TODO(Lukas) compute visc.
  const auto stressTensor = viscosity * ((2./3.) * div_v * I - 
					 (gradV + tarch::la::transpose(gradV)));

  const double factor = (gamma - 1)/(gasConstant);

#if DIMENSIONS == 2
  const double normJ = vars.j(0) * vars.j(0) + vars.j(1) * vars.j(1);
  const double normJ_dx = vars.j(0) * gradQ[idx_gradQ(0,1)] + vars.j(1) * gradQ[idx_gradQ(0,2)];
  const double normJ_dy = vars.j(0) * gradQ[idx_gradQ(1,1)] + vars.j(1) * gradQ[idx_gradQ(1,2)];
# elif DIMENSIONS == 3
  const double normJ = vars.j(0) * vars.j(0) + vars.j(1) * vars.j(1) + vars.j(2) * vars.j(2);
  const double normJ_dx = vars.j(0) * gradQ[idx_gradQ(0,1)] + vars.j(1) * gradQ[idx_gradQ(0,2)] + vars.j(2) * gradQ[idx_gradQ(0,3)];
  const double normJ_dy = vars.j(0) * gradQ[idx_gradQ(1,1)] + vars.j(1) * gradQ[idx_gradQ(1,2)] + vars.j(2) * gradQ[idx_gradQ(1,3)];
  const double normJ_dz = vars.j(0) * gradQ[idx_gradQ(2,1)] + vars.j(1) * gradQ[idx_gradQ(2,2)] + vars.j(2) * gradQ[idx_gradQ(2,3)];
#endif

  // Finally, assemble temperature gradient.
  // Note: idx for energy differs for 2d/3d!
#if DIMENSIONS == 2
  const double T_dx = factor * (invRho3 * gradQ[idx_gradQ(0,0)] * normJ - invRho2 * gradQ[idx_gradQ(0,0)] * vars.E() + invRho * gradQ[(idx_gradQ(0,3))] - invRho2 * normJ_dx);
  const double T_dy = factor * (invRho3 * gradQ[idx_gradQ(1,0)] * normJ - invRho2 * gradQ[idx_gradQ(1,0)] * vars.E() + invRho * gradQ[(idx_gradQ(1,3))] - invRho2 * normJ_dy);
#elif DIMENSIONS == 3
  const double T_dx = factor * (invRho3 * gradQ[idx_gradQ(0,0)] * normJ - invRho2 * gradQ[idx_gradQ(0,0)] * vars.E() + invRho * gradQ[(idx_gradQ(0,4))] - invRho2 * normJ_dx);
  const double T_dy = factor * (invRho3 * gradQ[idx_gradQ(1,0)] * normJ - invRho2 * gradQ[idx_gradQ(1,0)] * vars.E() + invRho * gradQ[(idx_gradQ(1,4))] - invRho2 * normJ_dy);
  const double T_dz = factor * (invRho3 * gradQ[idx_gradQ(2,0)] * normJ - invRho2 * gradQ[idx_gradQ(2,0)] * vars.E() + invRho * gradQ[(idx_gradQ(2,4))] - invRho2 * normJ_dz);
#endif

#if DIMENSIONS == 2
  const tarch::la::Vector<2, double> gradT =  {T_dx, T_dy};
  //const tarch::la::Vector<2, double> gradT = {0.0, 0.0};
#elif DIMENSIONS == 3
  const tarch::la::Vector<3, double> gradT =  {T_dx, T_dy, T_dz};
  //const tarch::la::Vector<3, double> gradT = {0.0, 0.0, 0.0};
#endif
  //temperatureDiffusion = false;
  const double kappa = (temperatureDiffusion) ? evaluateHeatConductionCoeff(viscosity) : 0.0;
  assertion2(std::isfinite(kappa), kappa, viscosity);

  // Full NS flux
  f.rho(vars.j());
  f.j(invRho * outerDot(vars.j(), vars.j()) + p*I + stressTensor);
//  f.E(
//      ((I * vars.E() + I * p + stressTensor) * (invRho * vars.j())) - kappa * gradT);

  // TODO(Lukas) fix for 3d
  // TODO(Lukas) names
  tarch::la::Matrix<1,2,double> t;
  t = ( invRho * vars.j(0)), (invRho * vars.j(1));

  const auto tt = t * (I * vars.E() + I * p + stressTensor);
  f.E(tt(0,0) - kappa * gradT[0], tt(0,1) - kappa * gradT[1]);


  for (int i = 0; i < vars.variables(); ++i) {
    const auto cond = std::isfinite(Q[i]);
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

