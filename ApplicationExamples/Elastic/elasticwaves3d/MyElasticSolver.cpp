#include "MyElasticSolver.h"

#include <memory>

elasticwaves3d::MyElasticSolver::init() {
  // This function is called inside the generated constructor.
  // @todo Please implement/augment if required
}

void elasticwaves3d::MyElasticSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)

  double* f = F[0];
  double* g = F[1];
  double* h = F[2];

  // @todo Please implement
  // f
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  f[3] = 0.0;
  f[4] = 0.0;
  f[5] = 0.0;
  f[6] = 0.0;
  f[7] = 0.0;
  f[8] = 0.0;
  // g
  // @todo Please implement
  g[0] = 0.0;
  g[1] = 0.0;
  g[2] = 0.0;
  g[3] = 0.0;
  g[4] = 0.0;
  g[5] = 0.0;
  g[6] = 0.0;
  g[7] = 0.0;
  g[8] = 0.0;
  // h
  // @todo Please implement
  h[0] = 0.0;
  h[1] = 0.0;
  h[2] = 0.0;
  h[3] = 0.0;
  h[4] = 0.0;
  h[5] = 0.0;
  h[6] = 0.0;
  h[7] = 0.0;
  h[8] = 0.0;
}

void elasticwaves3d::MyElasticSolver::boundaryConditions(
    const double* const x, const double t, const int faceIndex,
    const int normalNonZero, const double* const fluxIn,
    const double* const stateIn, double* const fluxOut, double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)

  // @todo Please implement
  // fluxOut
  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];
  fluxOut[5] = fluxIn[5];
  fluxOut[6] = fluxIn[6];
  fluxOut[7] = fluxIn[7];
  fluxOut[8] = fluxIn[8];
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
  stateOut[8] = stateIn[8];
}

void elasticwaves3d::MyElasticSolver::eigenvalues(const double* const Q,
                                                  const int normalNonZeroIndex,
                                                  double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement
  lambda[0] = 0.0;
  lambda[1] = 0.0;
  lambda[2] = 0.0;
  lambda[3] = 0.0;
  lambda[4] = 0.0;
  lambda[5] = 0.0;
  lambda[6] = 0.0;
  lambda[7] = 0.0;
  lambda[8] = 0.0;
}

bool elasticwaves3d::MyElasticSolver::hasToAdjustSolution(
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {
  // @todo Please implement
  return false;
}

void elasticwaves3d::MyElasticSolver::adjustedSolutionValues(
    const double* const x, const double w, const double t, const double dt,
    double* Q) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
  Q[7] = 0.0;
  Q[8] = 0.0;
}

exahype::solvers::Solver::RefinementControl
elasticwaves3d::MyElasticSolver::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void elasticwaves3d::MyElasticSolver::nonConservativeProduct(const double* const Q,
                                          const double* const gradQ,
                                          double* BgradQ) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  BgradQ[4] = 0.0;
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = 0.0;
  BgradQ[9] = 0.0;
  BgradQ[10] = 0.0;
  BgradQ[11] = 0.0;
  BgradQ[12] = 0.0;
  BgradQ[13] = 0.0;
  BgradQ[14] = 0.0;
  BgradQ[15] = 0.0;
  BgradQ[16] = 0.0;
  BgradQ[17] = 0.0;
  BgradQ[18] = 0.0;
  BgradQ[19] = 0.0;
  BgradQ[20] = 0.0;
  BgradQ[21] = 0.0;
  BgradQ[22] = 0.0;
  BgradQ[23] = 0.0;
  BgradQ[24] = 0.0;
  BgradQ[25] = 0.0;
  BgradQ[26] = 0.0;
}

void elasticwaves3d::MyElasticSolver::coefficientMatrix(const double* const Q,
                                              const int normalNonZero,
                                              double* Bn) {
  // Number of variables    = 9 (#unknowns + #parameters)
  // @todo Please implement
  Bn[0] = 0.0;
  Bn[1] = 0.0;
  Bn[2] = 0.0;
  Bn[3] = 0.0;
  Bn[4] = 0.0;
  Bn[5] = 0.0;
  Bn[6] = 0.0;
  Bn[7] = 0.0;
  Bn[8] = 0.0;
  Bn[9] = 0.0;
  Bn[10] = 0.0;
  Bn[11] = 0.0;
  Bn[12] = 0.0;
  Bn[13] = 0.0;
  Bn[14] = 0.0;
  Bn[15] = 0.0;
  Bn[16] = 0.0;
  Bn[17] = 0.0;
  Bn[18] = 0.0;
  Bn[19] = 0.0;
  Bn[20] = 0.0;
  Bn[21] = 0.0;
  Bn[22] = 0.0;
  Bn[23] = 0.0;
  Bn[24] = 0.0;
  Bn[25] = 0.0;
  Bn[26] = 0.0;
  Bn[27] = 0.0;
  Bn[28] = 0.0;
  Bn[29] = 0.0;
  Bn[30] = 0.0;
  Bn[31] = 0.0;
  Bn[32] = 0.0;
  Bn[33] = 0.0;
  Bn[34] = 0.0;
  Bn[35] = 0.0;
  Bn[36] = 0.0;
  Bn[37] = 0.0;
  Bn[38] = 0.0;
  Bn[39] = 0.0;
  Bn[40] = 0.0;
  Bn[41] = 0.0;
  Bn[42] = 0.0;
  Bn[43] = 0.0;
  Bn[44] = 0.0;
  Bn[45] = 0.0;
  Bn[46] = 0.0;
  Bn[47] = 0.0;
  Bn[48] = 0.0;
  Bn[49] = 0.0;
  Bn[50] = 0.0;
  Bn[51] = 0.0;
  Bn[52] = 0.0;
  Bn[53] = 0.0;
  Bn[54] = 0.0;
  Bn[55] = 0.0;
  Bn[56] = 0.0;
  Bn[57] = 0.0;
  Bn[58] = 0.0;
  Bn[59] = 0.0;
  Bn[60] = 0.0;
  Bn[61] = 0.0;
  Bn[62] = 0.0;
  Bn[63] = 0.0;
  Bn[64] = 0.0;
  Bn[65] = 0.0;
  Bn[66] = 0.0;
  Bn[67] = 0.0;
  Bn[68] = 0.0;
  Bn[69] = 0.0;
  Bn[70] = 0.0;
  Bn[71] = 0.0;
  Bn[72] = 0.0;
  Bn[73] = 0.0;
  Bn[74] = 0.0;
  Bn[75] = 0.0;
  Bn[76] = 0.0;
  Bn[77] = 0.0;
  Bn[78] = 0.0;
  Bn[79] = 0.0;
  Bn[80] = 0.0;
}
