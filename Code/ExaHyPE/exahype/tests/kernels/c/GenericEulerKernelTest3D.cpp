#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "kernels/aderdg/generic/Kernels.h"
#include "../testdata/generic_euler_testdata.h"

using std::cout;
using std::endl;

#if DIMENSIONS == 3

namespace exahype {
namespace tests {
namespace c {

void GenericEulerKernelTest::testFlux(const double *const Q, double *f,
                                      double *g, double *h) {
  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) *
      (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);

  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);

  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);

  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
}

void GenericEulerKernelTest::testPDEFluxes() {
  cout << "Test PDE-related functions, DIM=3" << endl;

  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5};  // pressure = 1.372
  double f[5], g[5], h[5];

  testFlux(Q, f, g, h);

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        f[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::f[i],
        i);
  }

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        g[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::g[i],
        i);
  }

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        h[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::h[i],
        i);
  }
}  // testPDEFluxes

void GenericEulerKernelTest::testVolumeIntegral() {
  cout << "Test volume integral, ORDER=3, DIM=3" << endl;

  // output:
  double *lduh = new double[320];  // intentionally left uninitialised

  // input:
  const double dx[3] = {0.05, 0.05, 0.05};  // mesh spacing
  double *lFhi = new double[960]();  // nVar * dim * nDOFx * nDOFy * nDOFz
  // lFhi = [ lFhi_x  | lFhi_y | lFhi_z ], 320 entries each
  double *lFhi_x = &lFhi[0];
  double *lFhi_y = &lFhi[320];
  double *lFhi_z = &lFhi[640];

  // seed direction
  for (int i = 0; i < 320; i += 5) {
    lFhi_x[i + 1] = 1.;
    lFhi_y[i + 2] = 1.;
    lFhi_z[i + 3] = 1.;
  }

  kernels::aderdg::generic::c::volumeIntegralNonlinear(
      lduh, lFhi, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i],
        ::exahype::tests::testdata::generic_euler::testVolumeIntegral::lduh[i],
        eps, i);
  }

  delete[] lduh;
  delete[] lFhi;
}  // testVolumeIntegral

void GenericEulerKernelTest::testSurfaceIntegral() {
  cout << "Test surface integral, ORDER=3, DIM=3" << endl;

  // inputs:
  const double dx[3] = {0.1, 0.1, 0.1};   // mesh spacing
  double *lFhbnd = new double[6 * 80]();  // 480
  double *FLeft = &lFhbnd[0];
  double *FRight = &lFhbnd[80];
  double *FFront = &lFhbnd[160];
  double *FBack = &lFhbnd[240];
  double *FBottom = &lFhbnd[320];
  double *FTop = &lFhbnd[400];

  for (int i = 0; i < 80; i += 5) {
    // in x orientation 1
    FLeft[i + 1] = 1.;
    FRight[i + 1] = 1.;
    // in y orientation 1
    FFront[i + 2] = 1.;
    FBack[i + 2] = 1.;
    // in z direction 1
    FBottom[i + 3] = 1.;
    FTop[i + 3] = 1.;
  }

  // in/out:
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i / 10.;
  }

  // lFhbnd = [ FLeft | FRight | FFront | FBack | FBottom | FTop ]
  kernels::aderdg::generic::c::surfaceIntegralNonlinear(
      lduh, lFhbnd, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i],
        ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh[i],
        eps, i);
  }

  delete[] lFhbnd;
  delete[] lduh;
}  // testSurfaceIntegral

void GenericEulerKernelTest::testRiemannSolver() {}

void GenericEulerKernelTest::testSolutionUpdate() {
  cout << "Test solution update, ORDER=3, DIM=3" << endl;

  // in/out:
  double *luh = new double[320]();
  for (int i = 0; i < 320; i += 5) {
    luh[i] = 1.0;
    luh[i + 4] = 2.5;
  }

  // inputs:
  const double dt = 1.40831757919882352703e-03;
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i;
  }

  kernels::aderdg::generic::c::solutionUpdate(luh, lduh, dt,
                                              5,  // getNumberOfVariables(),
                                              4   // getNodesPerCoordinateAxis()
                                              );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        luh[i],
        ::exahype::tests::testdata::generic_euler::testSolutionUpdate::luh[i],
        eps, i);
  }

  delete[] luh;
  delete[] lduh;
}  // testSolutionUpdate

void GenericEulerKernelTest::testSpaceTimePredictor() {
  cout << "Test space time predictor, ORDER=3, DIM=3" << endl;

  // inputs:
  // exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh[320 =
  // nVar * nDOFx * nDOFy * nDOFz]

  const double dx[3] = {0.5, 0.5, 0.5};
  const tarch::la::Vector<DIMENSIONS, double> dx_vector(0.5, 0.5, 0.5);
  const double timeStepSize = 1.267423918681417E-002;

  // local:
  double *lQi = new double[1280];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt
  double *lFi = new double[3840];  // nVar * nDOFx * nDOFy * nDOFz * nDOFt * dim

  // outputs:
  double *lQhi = new double[320];    // nVar * nDOFx * nDOFy * nDOFz
                                     // intentionally left uninitialised
  double *lFhi = new double[960];    // nVar * nDOFx * nDOFy * nDOFz * dim
  double *lQhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6
  double *lFhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6

  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<testFlux>(
      lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd,
      exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh,
      dx_vector, timeStepSize,
      5,  // numberOfVariables
      4   // basisSize
      );

  for (int i = 0; i < 1280; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQi[i], ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictor::lQi[i],
        5e-7, i);
  }

  for (int i = 0; i < 3840; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFi[i], ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictor::lFi[i],
        2.1e-5, i);
  }

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhi[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictor::lQhi[i],
        8e-8, i);
  }

  for (int i = 0; i < 960; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFhi[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictor::lFhi[i],
        3.4e-6, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictor::lQhbnd[i],
        1.1e-7, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictor::lFhbnd[i],
        4.2e-6, i);
  }

  delete[] lQi;
  delete[] lFi;
  delete[] lFhi;
  delete[] lQhi;
  delete[] lQhbnd;
  delete[] lFhbnd;
}  // testSpaceTimePredictor

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS==3
