/**
 * Initial Data for the EulerFlow.
 *
 * At the time of writing this code, there were no methods beyond environment
 * variables to put dynamical parameters into an ExaHyPE simulation. Nowadays
 * there is an API in ExaHyPE to access spec file parameters, cf. the MHD
 * application where this is done.
 *
 * For more details about the initial data implemented in ExaHyPE, see
 * https://gitlab.lrz.de/gi26det/ExaHyPE/wikis/list-of-benchmarks
 *
 **/

#include "InitialData.h"
#include "Primitives.h"
#include "MyEulerSolver.h"

#include <cstring>
#include <math.h>
#include <stdlib.h>

#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"

#include "MyEulerSolver.h"
#include "Primitives.h"

using namespace std;

// Storage for idfunc
InitialDataHandler idfunc;


/**
 * This function gives us vacuum values for Euler equations, ie. the trivial
 * solution. It fulfills the same signature as initial data functions, so it
 * can be used as a drop-in for boundary conditions.
 **/
void Vacuum(const double* x, double* Q, double t = 0.0) {
  static const double densvac = 1.0;
  static const double velvac = 0.0;
  static const double epsvac = 1.0;

  Q[0] = densvac;
  Q[1] = velvac;
  Q[2] = velvac;
  Q[3] = velvac;
  Q[4] = epsvac;
}

/**
 * Something to test the basic scheme. Advection and so.
 **/
void Polynomial(const double* const x, double* Q, double t = 0.0) {
  double V[Euler::MyEulerSolver::NumberOfVariables];
  V[0] = std::sin(x[0]/(2*M_PI));
  V[1] = 1;
  V[2] = 0;
  V[3] = 0;
  V[4] = 1;
  prim2con(Q,V);
}

/**
 * Source: MD ADERDG F90 CODE
 * ShuVortex produces primitive Variables which have to pass Prim2Con afterwards
 **/
void ShuVortex2D(const double* const x, double* Q, double t = 0.0) {
  static const double epsilon = 5.0;
  static const double pi = acos(-1.0);

  double r = sqrt(SQ(x[0] - t - 5.) + SQ(x[1] - t - 5.));
  double du = epsilon / 2. / pi * exp(0.5 * (1. - r * r)) * (5. - x[1] + t);
  double dv = epsilon / 2. / pi * exp(0.5 * (1. - r * r)) * (x[0] - 5. - t);
  double dTemp = -(eos_gamma - 1.) * SQ(epsilon) / 8. / eos_gamma / SQ(pi) *
                 exp(1. - r * r);
  double drho = pow(1. + dTemp, 1. / (eos_gamma - 1.)) - 1.;
  double dp = pow(1. + dTemp, eos_gamma / (eos_gamma - 1.)) - 1.;

  double V[Euler::MyEulerSolver::NumberOfVariables];
  V[0] = 1. + drho;
  V[1] = 1. + du;
  V[2] = 1. + dv;
  V[3] = 0.0;
  V[4] = 1. + dp;
  prim2con(Q, V);
}

/**
 * MovingGauss is a moving gaussian matter distribution where it is simple
 * to give an analytic result.
 **/
void MovingGauss(const double* const x, double* Q, double t = 0.0) {
  using namespace tarch::la;
  typedef Vector<DIMENSIONS,double> vec;

#if DIMENSIONS == 2
  vec xvec(x[0], x[1]);
  vec v0(0.5, 0);
  vec x0(0.5, 0.5);
#elif DIMENSIONS == 3
  vec xvec(x[0], x[1], x[2]);
  vec v0(0.5, 0.5, 0);
  vec x0(0.0, 0.0, 0.0);
#else
#error DIMENSIONS must be 2 or 3
#endif

  double width = 0.20;

  double V[Euler::MyEulerSolver::NumberOfVariables];
  V[0] = 0.5 +
         0.3 * exp(-norm2(xvec - x0 - v0 * t) /
                   pow(width, DIMENSIONS));
  V[1] = v0[0];
  V[2] = v0[1];
  V[3] = 0.;
  V[4] = 1.;  // pressure
  prim2con(Q, V);
}

/**
 * DiffusingGauss is not a consistent solution of Euler's equations, but instead
 * a perturbation which immediately changes rho, vx, vy and vz. It corresponds
 * to consistent initial data with roughly
 *   rho = exp( - (r - v0*t)**2)
 *   vx  = v0*t * exp(-y**2)  // actually even more complicated, it is more
 *   vy  = v0*t * exp(-x**2)  // a velocity on a kind of ring of radius v0*t
 *   E   = p/(gamma-1) + rho/2 * v**2 = p/(gamma-1) + alpha*exp(- (r-v0*t)**2)
 *   p   = 1
 * However, it is much more complicated to write the closed form solution
 *instead
 * of the perturbation approach. However, the closed form solution allows to
 * specify the solution at any time while while the initial perturbation form
 * does *not* allow to specify the solution.
 *
 * Attention: If the prim2con in the Primitives.h is used, for some reason it
 * surpresses this solution and rho=1, p=1, E=p/(gamma-1) is the resolut on the
 * whole grid. I am not sure why this happens, it seems to be an error.
 **/
void DiffusingGauss(const double* const x, double* Q, double t) {
  using namespace tarch::la;
  typedef Vector<DIMENSIONS,double> vecNd;

#if DIMENSIONS == 2
  vecNd xvec(x[0],x[1]);
  vecNd x0(0.5, 0.5);
#elif DIMENSIONS == 3
  vecNd xvec(x[0],x[1],x[2]);
  vecNd x0(0.5, 0.5, 0.5);
#else
#error DIMENSIONS must be 2 or 3
#endif
  
  
  double width = 0.25;

  Q[0] = 1.;
  Q[1] = 0.;
  Q[2] = 0.;
  Q[3] = 0.;
  Q[4] = 1. / (eos_gamma - 1) +
         exp(-norm2(xvec - x0) / pow(width, DIMENSIONS)) * 2;
}

void SmoothedSodShockTube(const double* const x, double* Q, double t) {
  // Velocities are set to zero.
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  // (Smoothed out) jump in density and pressure.
  //
  // Original distribution:
  // if(x[0] < 0.5) {
  //    Q[0] = 1.0;
  //    Q[4] = 1.0;
  //  } else {
  //    Q[0] = 0.125;
  //    Q[4] = 0.1;
  //  }

  const double steepness = 20;

  Q[0] = 0.1   + 0.5 * ( 1.0 + std::tanh(steepness*(x[0]-0.5)) ) * (1.0-0.1);
  Q[4] = 0.125 + 0.5 * ( 1.0 + std::tanh(steepness*(x[0]-0.5)) ) * (1.0-0.125);
}

