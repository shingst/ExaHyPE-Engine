#include "PDE.h"
#include "Scenarios/Scenario.h"

NavierStokes::PDE::PDE() :
  PDE(0.1, 10000, 1.4, 0.7, 1, 1.4, 0.4) {
}

NavierStokes::PDE::PDE(double referenceViscosity, NavierStokes::Scenario &scenario) :
  referenceViscosity(referenceViscosity),
  referencePressure(scenario.getReferencePressure()),
  gamma(scenario.getGamma()),
  Pr(scenario.getPr()),
  c_v(scenario.getC_v()),
  c_p(scenario.getC_p()),
  gasConstant(scenario.getGasConstant()) { }

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

#if DIMENSIONS == 2
void NavierStokes::PDE::evaluateFlux(const double* Q, const double* gradQ, double** F, bool temperatureDiffusion) const {
  ReadOnlyVariables vars(Q);
  //Fluxes f(F);

  auto idxF = kernels::idx2(vars.SizeVariables, DIMENSIONS);
  auto idxGradQ = kernels::idx2(DIMENSIONS, vars.SizeVariables);

  // Euler:
  const auto invRho = 1/vars.rho(); // Q(1)/(Q(1)*Q(1)+epsilon)
  const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j());
  //p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )

  double* f = F[0];
  double* g = F[1];
#if DIMENSIONS == 3
  double* h = F[2];
#endif

  const auto rho = NavierStokesSolverDG_Variables::shortcuts::rho;
  const auto j = NavierStokesSolverDG_Variables::shortcuts::j;
  const auto E = NavierStokesSolverDG_Variables::shortcuts::E;

  /*
  f(1) = Q(2) !- EQN%mu*gradQ(1,1)
  f(2) = irho*Q(2)*Q(2) + p
  f(3) = irho*Q(2)*Q(3)
  f(4) = irho*Q(2)*Q(4)
  f(5) = irho*Q(2)*(Q(5)+p)
   */

  f[rho] = Q[j];
  f[j] = invRho * Q[j] * Q[j] + p;
  f[j+1] = invRho * Q[j] * Q[j+1];
#if DIMENSIONS == 3
  f[j+2] = invRho * Q[j] * Q[j+2];
#endif
  f[E] = invRho * Q[j] * (Q[E] + p);

  /*
  g(1) = Q(3) !- EQN%mu*gradQ(1,2)
  g(2) = irho*Q(3)*Q(2)
  g(3) = irho*Q(3)*Q(3) + p
  g(4) = irho*Q(3)*Q(4)
  g(5) = irho*Q(3)*(Q(5)+p)
   */

  g[rho] = Q[j+1];
  g[j] = invRho * Q[j+1] * Q[j];
  g[j+1] = invRho * Q[j+1] * Q[j+1] + p;
#if DIMENSIONS == 3
  g[j+2] = invRho * Q[j+1] * Q[j+2];
#endif
  g[E] = invRho * Q[j+1] * (Q[E]+p);

  /*
  h(1) = Q(4) !- EQN%mu*gradQ(1,3)
  h(2) = irho*Q(4)*Q(2)
  h(3) = irho*Q(4)*Q(3)
  h(4) = irho*Q(4)*Q(4) + p
  h(5) = irho*Q(4)*(Q(5)+p)
  */

#if DIMENSIONS == 3
  // TODO
  h[rho] = Q[j+2];
  h[j] = invRho * Q[j+2] * Q[j];
  h[j+1] = invRho * Q[j+2] * Q[j+1] + p;
  h[j+2] = invRho * Q[j+2] * Q[j+2];
  h[E] = invRho * Q[j+2] * (Q[E]+p);
#endif

     // Viscous:
  /*
 iRho  = 1./Q(1)
 uu    = Q(2)*iRho
 vv    = Q(3)*iRho
 ww    = Q(4)*iRho
   */
  const auto uu = invRho * Q[j];
  const auto vv = invRho * Q[j+1];
#if DIMENSIONS == 3
  const auto ww = invRho * Q[j+2];
#endif

  /*
  mu    = EQN%mu
  kappa = EQN%kappa
  */
  // TODO(Lukas) Support different visc. models
  const auto mu = referenceViscosity;
  const auto kappa = evaluateHeatConductionCoeff(mu);

  /*
  uux  = iRho*( gradQ(2,1) - uu*gradQ(1,1) )
  vvx  = iRho*( gradQ(3,1) - vv*gradQ(1,1) )
  wwx  = iRho*( gradQ(4,1) - ww*gradQ(1,1) )
  uuy  = iRho*( gradQ(2,2) - uu*gradQ(1,2) )
  vvy  = iRho*( gradQ(3,2) - vv*gradQ(1,2) )
  wwy  = iRho*( gradQ(4,2) - ww*gradQ(1,2) )
  uuz  = iRho*( gradQ(2,3) - uu*gradQ(1,3) )
  vvz  = iRho*( gradQ(3,3) - vv*gradQ(1,3) )
  wwz  = iRho*( gradQ(4,3) - ww*gradQ(1,3) )
  */

  // Derivatives of velocities:
  const auto uux  = invRho * (gradQ[idxGradQ(0,j+0)] - uu * gradQ[idxGradQ(0, rho)]);
  const auto vvx  = invRho * (gradQ[idxGradQ(0,j+1)] - vv * gradQ[idxGradQ(0, rho)]);
#if DIMENSIONS == 3
  const auto wwx  = invRho * (gradQ[idxGradQ(0,j+2)] - ww * gradQ[idxGradQ(0, rho)]);
#endif

  const auto uuy  = invRho * (gradQ[idxGradQ(1,j+0)] - uu * gradQ[idxGradQ(1, rho)]);
  const auto vvy  = invRho * (gradQ[idxGradQ(1,j+1)] - vv * gradQ[idxGradQ(1, rho)]);
#if DIMENSIONS == 3
  const auto wwy  = invRho * (gradQ[idxGradQ(1,j+2)] - ww * gradQ[idxGradQ(1, rho]));
#endif

# if DIMENSIONS == 3
  const auto uuz  = invRho * (gradQ[idxGradQ(2,j+0)] - uu * gradQ[(idxGradQ(2, rho)]);
  const auto vvz  = invRho * (gradQ[idxGradQ(2,j+1)] - vv * gradQ[(idxGradQ(2, rho)]);
  const auto wwz  = invRho * (gradQ[idxGradQ(2,j+2)] - ww * gradQ[(idxGradQ(2, rho)]);
#endif

  /*
  icv   = 1./EQN%cv
  iRho2 = iRho*iRho
  iRho3 = iRho2*iRho
  */
  const auto invCv = 1/c_v;
  const auto invRho2 = invRho * invRho;
  const auto invRho3 = invRho2 * invRho;

  /*
  dTdW1 = - Q(5)*iRho2 + iRho3*( Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4) )
  dTdW2 = - Q(2)*iRho2
  dTdW3 = - Q(3)*iRho2
  dTdW4 = - Q(4)*iRho2
  */
#if DIMENSIONS == 2
  const auto dTdW1 = -1 * Q[E] * invRho2 +
          invRho3 * (Q[j]*Q[j] + Q[j+1] * Q[j+1]);
#else
  const auto dTdW1 = -1 * Q[E] * invRho2 +
          invRho3 * (Q[j]*Q[j] + Q[j+1] * Q[j+1] + Q[j+2] * Q[j+2]);
#endif
  const auto dTdW2 = -1 * Q[j] * invRho2;
  const auto dTdW3 = -1 * Q[j+1] * invRho2;
#if DIMENSIONS == 3
  const auto dTdW4 = -1 * Q[j+2] * invRho2;
#endif

  /*
  Tx = icv*( dTdW1*gradQ(1,1) + dTdW2*gradQ(2,1) + dTdW3*gradQ(3,1) + dTdW4*gradQ(4,1) + iRho*gradQ(5,1) )
  Ty = icv*( dTdW1*gradQ(1,2) + dTdW2*gradQ(2,2) + dTdW3*gradQ(3,2) + dTdW4*gradQ(4,2) + iRho*gradQ(5,2) )
  Tz = icv*( dTdW1*gradQ(1,3) + dTdW2*gradQ(2,3) + dTdW3*gradQ(3,3) + dTdW4*gradQ(4,3) + iRho*gradQ(5,3) )
  */
# if DIMENSIONS == 2
  const auto Tx = invCv * (dTdW1 * gradQ[idxGradQ(0, rho)] +
          dTdW2 * gradQ[idxGradQ(0, j)] +
          dTdW3 * gradQ[idxGradQ(0, j+1)] +
          invRho * gradQ[idxGradQ(0,E)]);

  const auto Ty = invCv * (dTdW1 * gradQ[idxGradQ(1, rho)] +
          dTdW2 * gradQ[idxGradQ(1, j)] +
          dTdW3 * gradQ[idxGradQ(1, j+1)] +
          invRho * gradQ[idxGradQ(1,E)]);
#else
  // TODO(Lukas): Untested for 3D!
   const auto Tx = invCv * (dTdW1 * gradQ[idxGradQ(0, rho)] +
          dTdW2 * gradQ[idxGradQ(0, j)] +
          dTdW3 * gradQ[idxGradQ(0, j+1)] +
          dTdW4 * gradQ[idxGradQ(0, j+2)] +
          invRho * gradQ[idxGradQ(0,E)]);

  const auto Ty = invCv * (dTdW1 * gradQ[idxGradQ(1, rho)] +
          dTdW2 * gradQ[idxGradQ(1, j)] +
          dTdW3 * gradQ[idxGradQ(1, j+1)] +
          dTdW4 * gradQ[idxGradQ(1, j+2)] +
          invRho * gradQ[idxGradQ(1,E)]);

  const auto Tz = invCv * (dTdW1 * gradQ[idxGradQ(2, rho)] +
          dTdW2 * gradQ[idxGradQ(2, j)] +
          dTdW3 * gradQ[idxGradQ(2, j+1)] +
          dTdW4 * gradQ[idxGradQ(2, j+2)] +
          invRho * gradQ[idxGradQ(2,E)]);
#endif
  /*
  divV23  = 2./3.*(uux + vvy + wwz)
  */
#if DIMENSIONS == 2
  const auto divV23 = 2./3. * (uux + vvy);
#else
  const auto divV23 = 2./3. * (uux + vvy + wwz);
#endif

  double Fv[vars.Size] = {0.0};
  double Gv[vars.Size] = {0.0};
  double Hv[vars.Size] = {0.0};
  /*
  Fv(1) = 0.
  Fv(2) = mu*( 2*uux - divV23 )
  Fv(3) = mu*(   uuy + vvx    )
  Fv(4) = mu*(   uuz + wwx    )
  Fv(5) = Fv(2)*uu + Fv(3)*vv + Fv(4)*ww + kappa*Tx
  */
  Fv[rho] = 0.0;
  Fv[j+0] = mu * (2 * uux - divV23);
  Fv[j+1] = mu * (uuy + vvx);
#if DIMENSIONS == 2
  Fv[E] = Fv[j] * uu + Fv[j+1] * vv + kappa * Tx;
#else
  Fv[j+2] = mu * (uuz + wwx);
  Fv[E] = Fv[j] * uu + Fv[j+1] * vv + Fv[j+2] * ww + kappa * Tx;
#endif

  /*
  Gv(1) = 0.
  Gv(2) = Fv(3)
  Gv(3) = mu*( 2*vvy - divV23 )
  Gv(4) = mu*(   vvz + wwy    )
  Gv(5) = Gv(2)*uu + Gv(3)*vv + Gv(4)*ww + kappa*Ty
  */
  Gv[rho] = 0.0;
  Gv[j+0] = Fv[j+1];
  Gv[j+1] = mu * (2 * vvy - divV23);
#if DIMENSIONS == 2
  Gv[E] = Gv[j] * uu + Gv[j+1] * vv + kappa * Ty;
#else
  Gv[j+2] = mu * (vvz + wwy);
  Gv[E] = Gv[j] * uu + Gv[j+1] * vv + Gv[j+2] * ww + kappa * Ty;
#endif

  // TODO(Lukas) Viscous flux for 3D!
  /*
  Hv(1) = 0.
  Hv(2) = Fv(4)
  Hv(3) = Gv(4)
  Hv(4) = mu*( 2*wwz - divV23 )
  Hv(5) = Hv(2)*uu + Hv(3)*vv + Hv(4)*ww + kappa*Tz
  */

  /*
  f = f - Fv
  g = g - Gv
  h = h - Hv
  */
  for (int i = 0; i < vars.Size; ++i) {
    f[i] -= Fv[i];
    g[i] -= Gv[i];
#if DIMENSIONS == 3
    h[i] -= Hv[i];
#endif
  }
}

#else
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

  //const double factor = (gamma - 1)/(gasConstant);
  const double factor = 1/c_v;
  // Alternatively factor = 1/cv

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
  f.j(outerDot(invRho * vars.j(), vars.j()) + p*I + stressTensor);
//  f.E(
//      ((I * vars.E() + I * p + stressTensor) * (invRho * vars.j())) - kappa * gradT);

  // TODO(Lukas) fix for 3d
  // TODO(Lukas) names
  tarch::la::Matrix<1,2,double> t;
  t = ( invRho * vars.j(0)), (invRho * vars.j(1));

  const tarch::la::Matrix<1,2,double> tt = t * (I * vars.E() + I * p + stressTensor);
  f.E(tt(0,0) - kappa * gradT[0], tt(0,1) - kappa * gradT[1]);

}
#endif
