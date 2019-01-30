#include "MyElasticSolver.h"

#include <memory>

elasticwaves2d::MyElasticSolver::init() {
  // This function is called inside the generated constructor.
  // @todo Please implement/augment if required
}

void elasticwaves2d::MyElasticSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  double* f = F[0];
  double* g = F[1];

  // material parameters
  double irho = 1.0;    /* density */
  double ilambda = 1.0; /*  1st Lame parameter */
  double imu = 1.0;     /* 2nd Lame parameter */

  // @todo Please implement
  // f
  f[0] = 1.0 / irho * Q[2];
  f[1] = 1.0 / irho * Q[4];
  f[2] = (2.0 * imu + ilambda) * Q[0];
  f[3] = ilambda * Q[0];
  f[4] = imu * Q[1];
  // g
  // @todo Please implement
  g[0] = 1.0 / irho * Q[4];
  g[1] = 1.0 / irho * Q[3];
  g[2] = ilambda * Q[1];
  g[3] = (2.0 * imu + ilambda) * Q[1];
  g[4] = imu * Q[0];
}

void elasticwaves2d::MyElasticSolver::eigenvalues(const double* const Q,
                                                  const int normalNonZeroIndex,
                                                  double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement

  // material parameters
  double irho = 1.0;    /* density */
  double ilambda = 1.0; /*  1st Lame parameter */
  double imu = 1.0;     /* 2nd Lame parameter */

  // extract p-wave speed (cp) and s-wave speed (cs)
  double cp = std::sqrt((2.0 * imu + ilambda) / irho);
  double cs = std::sqrt(imu / irho);

  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = cp;
  lambda[3] = cs;
  lambda[4] = 0.0;
}

bool elasticwaves2d::MyElasticSolver::hasToAdjustSolution(
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  // @todo Please implement
  if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}

void elasticwaves2d::MyElasticSolver::adjustedSolutionValues(
    const double* const x, const double w, const double t, const double dt,
    double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  // @todo Please implement

  // initialize particle velocities with a Gaussian
  Q[0] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
                  (0.05 * 0.05)) *
         1.0e-3;
  Q[1] = std::exp(-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) /
                  (0.05 * 0.05)) *
         1.0e-3;

  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
}

exahype::solvers::Solver::RefinementControl
elasticwaves2d::MyElasticSolver::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}
