#include "MyElasticSolver.h"

#include <memory>

#include "kernels/GaussLegendreBasis.h"
#include "peano/utils/Loop.h"

#include "kernels/aderdg/generic/Kernels.h"

Elastic::MyElasticSolver::init() {
  // This function is called inside the generated constructor.
  // @todo Please implement/augment if required
}

void Elastic::MyElasticSolver::flux(const double* const Q, double** const F) {
  // Dimensions             = 2
  // Number of variables    = 8 (#unknowns + #parameters)

  double rho; /* density */
  double lam; /* 1st Lame parameter */
  double mu;  /* 2nd Lame parameter */

  rho = Q[5];
  lam = Q[6];
  mu = Q[7];

  double* f = F[0];
  double* g = F[1];

  // @todo Please implement
  // f
  // @todo Please implement
  // f
  f[0] = 1.0 / rho * Q[2];
  f[1] = 1.0 / rho * Q[4];
  f[2] = (2.0 * mu + lam) * Q[0];
  f[3] = lam * Q[0];
  f[4] = mu * Q[1];
  f[5] = 0.0;
  f[6] = 0.0;
  f[7] = 0.0;
  // g
  // @todo Please implement
  g[0] = 1.0 / rho * Q[4];
  g[1] = 1.0 / rho * Q[3];
  g[2] = lam * Q[1];
  g[3] = (2.0 * mu + lam) * Q[1];
  g[4] = mu * Q[0];
  g[5] = 0.0;
  g[6] = 0.0;
  g[7] = 0.0;
}

void Elastic::MyElasticSolver::algebraicSource(const double* const Q, double* const S) {
  // Number of variables = 8 + 0
  // @todo Please implement
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = 0.0;
  S[6] = 0.0;
  S[7] = 0.0;
}

void Elastic::MyElasticSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 2
  // Number of variables    = 8 (#unknowns + #parameters)

  // @todo Please implement
  // fluxOut
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = 0 * fluxIn[2];
  fluxOut[3] = 0 * fluxIn[3];
  fluxOut[4] = 0 * fluxIn[4];
  fluxOut[5] = fluxIn[5];
  fluxOut[6] = fluxIn[6];
  fluxOut[7] = fluxIn[7];
  // stateOut
  // @todo Please implement
  stateOut[0] = 0 * stateIn[0];
  stateOut[1] = 0 * stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];
  stateOut[5] = stateIn[5];
  stateOut[6] = stateIn[6];
  stateOut[7] = stateIn[7];

  if (faceIndex == 0) {
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];
    fluxOut[5] = fluxIn[5];
    fluxOut[6] = fluxIn[6];
    fluxOut[7] = fluxIn[7];
    // stateOut
    // @todo Please implement
    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = stateIn[4];
    stateOut[5] = stateIn[5];
    stateOut[6] = stateIn[6];
    stateOut[7] = stateIn[7];

    fluxOut[0] = 0.0 * fluxIn[0];
    fluxOut[1] = 0.0 * fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];

    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = 0.0 * stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = 0.0 * stateIn[4];
  }
}

void Elastic::MyElasticSolver::eigenvalues(const double* const Q,
                                           const int normalNonZeroIndex,
                                           double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 8 (#unknowns + #parameters)
  // @todo Please implement
  double rho; /* density */
  double lam; /*  1st Lame parameter */
  double mu;  /* 2nd Lame parameter */

  rho = Q[5];
  lam = Q[6];
  mu = Q[7];

  // extract p-wave speed (cp) and s-wave speed (cs)
  double cp = std::sqrt((2.0 * mu + lam) / rho);
  double cs = std::sqrt(mu / rho);

  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = cp;
  lambda[3] = cs;
  lambda[4] = 0.0;
  lambda[5] = 0.0;
  lambda[6] = 0.0;
  lambda[7] = 0.0;
}

bool Elastic::MyElasticSolver::hasToAdjustSolution(
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  // @todo Please implement
  // if (tarch::la::equals(t, 0.0)) {
  return true;
  //}
  // return false;
}

void Elastic::MyElasticSolver::adjustedSolutionValues(const double* const x,
                                                      const double w,
                                                      const double t,
                                                      const double dt,
                                                      double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 8 (#unknowns + #parameters)
  // @todo Please implement
  double rho; /* density */
  double cp;  /*  p-wave speed */
  double cs;  /* s-wave speed */

  if (tarch::la::equals(t, 0.0)) {
    if (true) {
      rho = 2.6;  // gm/cm^3
      cs = 2.0;   // km/s
      cp = 4.0;   // km/s

      Q[0] =
          0.0 *
          std::exp(-((x[0] - 5) * (x[0] - 5) + (x[1] - 10.0) * (x[1] - 10.0)) /
                   (1.0 * 1.0)) *
          10.0;
      Q[1] =
          0.0 *
          std::exp(-((x[0] - 5) * (x[0] - 5) + (x[1] - 10.0) * (x[1] - 10.0)) /
                   (1.0 * 1.0)) *
          10.0;
      Q[2] = 0.0;
      Q[3] = 0.0;
      Q[4] = 0.0;
      Q[5] = rho;                            // rho
      Q[6] = rho * (cp * cp - 2 * cs * cs);  // lambda (MPa)
      Q[7] = rho * cs * cs;                  // mu (MPa)

      if (x[0] > 1.0) {
        rho = 2.7;   // gm/cm^3
        cs = 3.464;  // km/s
        cp = 6.0;    // km/s

        Q[5] = rho;                            // rho
        Q[6] = rho * (cp * cp - 2 * cs * cs);  // lambda (MPa)
        Q[7] = rho * cs * cs;                  // mu     (MPa)
      }
    }
  }
}

exahype::solvers::Solver::RefinementControl
Elastic::MyElasticSolver::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Elastic::MyElasticSolver::pointSources(
    double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt,
    int numberOfVariables, int numberOfParameters, int _order) {
  // SKIP SOURCE
  // return;

  const tarch::la::Vector<DIMENSIONS, double> _x(2.0, 5.0);
  const tarch::la::Vector<DIMENSIONS, double> offsetOfPatch = center - dx / 2.;
  const tarch::la::Vector<DIMENSIONS, double> sizeOfPatch = dx;  // / 2.;

  /** Addition of point source terms **/
  /* if (
   tarch::la::allSmaller(offsetOfPatch,_x)
   &&
   tarch::la::allGreater(offsetOfPatch+sizeOfPatch,_x)
   ) { */

  // Map coordinate vector x onto reference element.
  tarch::la::Vector<DIMENSIONS, double> xRef = _x - offsetOfPatch;
  xRef(0) /= sizeOfPatch(0);
  xRef(1) /= sizeOfPatch(1);
#ifdef Dim3
  xRef(2) /= sizeOfPatch(2);
#endif

  if (xRef(0) >= 0. && xRef(1) >= 0. && xRef(0) <= 1.0 && xRef(1) <= 1.0) {
    //
    double value = 0.0;
    double sourceInt = 0.0;
    sourceInt =
        ((1.0 - (1.0 + t / 0.1) * (std::exp(-t / 0.1))) -
         (1.0 - (1.0 + (t - dt) / 0.1) * (std::exp(-(t - dt) / 0.1))));  // 1.;
    // sourceInt = 0.5*(t/(0.1*0.1)*(std::exp(-t/0.1)) +
    // (t-dt)/(0.1*0.1)*(std::exp(-(t-dt)/0.1)));
    if (tarch::la::equals(t, 0.0)) {
      sourceInt = 0.0;
      // sourceInt = (1.0-(1.0+t/0.1)*(std::exp(-t/0.1))); //1.;
    }
    //
    // std::cout<<dx[0]<< " "<<dx[1]<<" "<<dx[2]<<std::endl;
    // exit(-1);
    for (int unknown = 0; unknown < numberOfVariables /*numberOfParameters*/;
         unknown++) {
      if (unknown == 4) {
        dfor(ii, _order + 1) {  // Gauss-Legendre node indices
          int iGauss = peano::utils::dLinearisedWithoutLookup(ii, _order + 1);
          value = kernels::basisFunctions[_order][ii(0)](xRef(0)) *
                  kernels::basisFunctions[_order][ii(1)](xRef(1)) *
#ifdef Dim3
                  kernels::basisFunctions[_order][ii(2)](xRef(2)) *
#endif
                  sourceInt;
          value /= kernels::legendre::weights[_order][ii(0)] *
                   kernels::legendre::weights[_order][ii(1)] *
#ifdef Dim3
                   kernels::legendre::weights[_order][ii(2)] *
#endif
                   1.0;
          luh[iGauss] += value;
        }
      }
    }
  }  // endif xRef inside the domain
}
