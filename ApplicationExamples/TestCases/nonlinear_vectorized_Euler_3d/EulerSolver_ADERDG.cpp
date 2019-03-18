/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include <algorithm>
#include <string>
#include <math.h>

// for the intrinsics, required by gcc7
#include <immintrin.h>
 
#include "EulerSolver_ADERDG.h"
//#include "EulerSolver_ADERDG_Variables.h"
#include "tarch/la/MatrixVectorOperations.h"

#include "kernels/GaussLegendreQuadrature.h"

tarch::logging::Log Euler::EulerSolver_ADERDG::_log("Euler::EulerSolver_ADERDG");


void Euler::EulerSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {}

// required for boundary condition
void Euler::EulerSolver_ADERDG::flux(const double* const Q, double** const F) {
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];
  #if DIMENSIONS==3
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3];
  #else
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);

  // col 1
  F[0][0] = Q[1];
  F[0][1] = irho*Q[1]*Q[1] + p;
  F[0][2] = irho*Q[2]*Q[1];
  F[0][3] = irho*Q[3]*Q[1];
  F[0][4] = irho*(Q[4]+p)*Q[1];

  // col 2
  F[1][0] = Q[2];
  F[1][1] = irho*Q[1]*Q[2];
  F[1][2] = irho*Q[2]*Q[2] + p;
  F[1][3] = irho*Q[3]*Q[2];
  F[1][4] = irho*(Q[4]+p)*Q[2];

  #if DIMENSIONS==3
  // col 3
  F[2][0] = Q[3];
  F[2][1] = irho*Q[1]*Q[3];
  F[2][2] = irho*Q[2]*Q[3];
  F[2][3] = irho*Q[3]*Q[3] + p;
  F[2][4] = irho*(Q[4]+p)*Q[3];
  #endif
}


// #define VECTSIZE 4
void Euler::EulerSolver_ADERDG::flux_vect(const double* const * const restrict Q, double* const * const * const restrict F, int size){


	__assume_aligned(F, 32);
	__assume_aligned(Q, 32);
	
	int i, j, k;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7;
	__m256d ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;
  
	double gamma[1], gamma1[1]; 
	gamma[0]   = 1.4; 
	gamma1[0]  = gamma[0] - 1.0;			// precompute gamma - 1 

	ymm0 = _mm256_load_pd(&Q[0][0]);		// load density
	ymm1 = _mm256_set1_pd(1.0);				// set unit 
	ymm2 = _mm256_div_pd(ymm1, ymm0);		// compute density^{-1} 

	ymm3 = _mm256_load_pd(&Q[1][0]);		// load x momentum 
	ymm4 = _mm256_mul_pd(ymm3, ymm2);		// compute x velocity 
	ymm5 = _mm256_load_pd(&Q[2][0]);		// load y momentum 
	ymm6 = _mm256_mul_pd(ymm5, ymm2);		// compute y velocity 
	ymm7 = _mm256_load_pd(&Q[3][0]);		// load z momentum 
	ymm8 = _mm256_mul_pd(ymm7, ymm2);		// compute z velocity 
	
	ymm9  = _mm256_load_pd(&Q[4][0]);		 // load total energy 
	ymm1  = _mm256_broadcast_sd(&gamma1[0]); // broadcast gamma-1 
	ymm15 = _mm256_set1_pd(-0.5);			 // set -1/2 
	// 
	// compute the pressure via AVX instructions (AVX2 does not make sense here, since we also need the intermediate result of the products 
	// p(:) = (EQN%gamma - 1)*(Q(:, 5) - 0.5*Q(:, 1)*(uu(:)*uu(:) + vv(:)*vv(:) + ww(:)*ww(:)))
	// 
	ymm11 = _mm256_mul_pd(ymm3, ymm4);				// compute rho*u*u 
	ymm12 = _mm256_mul_pd(ymm5, ymm6);				// compute rho*v*v 
	ymm13 = _mm256_mul_pd(ymm7, ymm8);				// compute rho*w*w 
	ymm10 = _mm256_add_pd(ymm11,ymm12);				// compute the kinetic energy * 2 
	ymm10 = _mm256_add_pd(ymm10,ymm13); 
	ymm10 = _mm256_mul_pd(ymm10, ymm15);			// divide the result by 2 and invert sign 
	ymm10 = _mm256_add_pd(ymm10, ymm9);				// subtract the kinetic energy density from the total energy density to get the internal energy density 
	ymm10 = _mm256_mul_pd(ymm10, ymm1);				// multiply the internal energy density with gamma-1 to get the pressure 
	// now compute the tensor of the nonlinear convective terms via AVX instructions 
	ymm0  = _mm256_mul_pd(ymm3, ymm6);				// compute rho*u*v  
	ymm1  = _mm256_mul_pd(ymm3, ymm8);				// compute rho*u*w  
	ymm2  = _mm256_mul_pd(ymm5, ymm8);				// compute rho*v*w 
	// add the pressure to the diagonal 
	ymm11 = _mm256_add_pd(ymm11, ymm10);			// compute rho*u*u + p 
	ymm12 = _mm256_add_pd(ymm12, ymm10);			// compute rho*v*v + p 
	ymm13 = _mm256_add_pd(ymm13, ymm10);			// compute rho*w*w + p 
	ymm14 = _mm256_add_pd(ymm9, ymm10);				// compute rhoE + p 

	// store the results in F 
	ymm15 = _mm256_mul_pd(ymm14, ymm4);				// compute u*(rhoE+p) 
	_mm256_store_pd(&F[0][0][0], ymm3);				// rho*u 
	_mm256_store_pd(&F[0][1][0], ymm11);			// rho*u*u + p   
	_mm256_store_pd(&F[0][2][0], ymm0);				// rho*u*v 
	_mm256_store_pd(&F[0][3][0], ymm1);				// rho*u*w 
	_mm256_store_pd(&F[0][4][0], ymm15);			// u*(rhoE+p)  

	ymm15 = _mm256_mul_pd(ymm14, ymm6);				// compute v*(rhoE+p) 
	_mm256_store_pd(&F[1][0][0], ymm5);				// rho*v 
	_mm256_store_pd(&F[1][1][0], ymm0);				// rho*v*u 
	_mm256_store_pd(&F[1][2][0], ymm12);			// rho*v*v + p  
	_mm256_store_pd(&F[1][3][0], ymm2);				// rho*v*w 
	_mm256_store_pd(&F[1][4][0], ymm15);			// v*(rhoE+p)  
#if DIMENSIONS==3
	ymm15 = _mm256_mul_pd(ymm14, ymm8);				// compute w*(rhoE+p) 
	_mm256_store_pd(&F[2][0][0], ymm7);				// rho*w 
	_mm256_store_pd(&F[2][1][0], ymm1);				// rho*w*u 
	_mm256_store_pd(&F[2][2][0], ymm2);				// rho*w*v   
	_mm256_store_pd(&F[2][3][0], ymm13);			// rho*w*w + p 
	_mm256_store_pd(&F[2][4][0], ymm15);			// v*(rhoE+p)  
#endif

 /*
 
  constexpr double gamma = 1.4;
  double irho[VECTSIZE] __attribute__((aligned(ALIGNMENT)));
  double j2[VECTSIZE] __attribute__((aligned(ALIGNMENT)));
  double p[VECTSIZE] __attribute__((aligned(ALIGNMENT)));
  
  #pragma vector aligned
  #pragma ivdep
  for(int i=0; i<VECTSIZE; i++){
    irho[i] = 1./Q[0][i];
    #if DIMENSIONS==3
    j2[i] = Q[1][i]*Q[1][i]+Q[2][i]*Q[2][i]+Q[3][i]*Q[3][i];
    #else
    j2[i] = Q[1][i]*Q[1][i]+Q[2][i]*Q[2][i];
    #endif
    p[i] = (gamma-1) * (Q[4][i] - 0.5*irho[i]*j2[i]);

    // col 1
    F[0][0][i] = Q[1][i];
    F[0][1][i] = irho[i]*Q[1][i]*Q[1][i] + p[i];
    F[0][2][i] = irho[i]*Q[2][i]*Q[1][i];
    F[0][3][i] = irho[i]*Q[3][i]*Q[1][i];
    F[0][4][i] = irho[i]*(Q[4][i]+p[i])*Q[1][i];

    // col 2
    F[1][0][i] = Q[2][i];
    F[1][1][i] = irho[i]*Q[1][i]*Q[2][i];
    F[1][2][i] = irho[i]*Q[2][i]*Q[2][i] + p[i];
    F[1][3][i] = irho[i]*Q[3][i]*Q[2][i];
    F[1][4][i] = irho[i]*(Q[4][i]+p[i])*Q[2][i];

    #if DIMENSIONS==3
    // col 3
    F[2][0][i] = Q[3][i];
    F[2][1][i] = irho[i]*Q[1][i]*Q[3][i];
    F[2][2][i] = irho[i]*Q[2][i]*Q[3][i];
    F[2][3][i] = irho[i]*Q[3][i]*Q[3][i] + p[i];
    F[2][4][i] = irho[i]*(Q[4][i]+p[i])*Q[3][i];
    #endif
  }
  */
}

void Euler::EulerSolver_ADERDG::eigenvalues(const double* const Q, const int direction, double* const lambda) {
  constexpr double gamma = 1.4;
  const double irho = 1./Q[0];
  #if DIMENSIONS==3
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3];
  #else
  const double j2 = Q[1]*Q[1]+Q[2]*Q[2];
  #endif
  const double p = (gamma-1) * (Q[4] - 0.5*irho*j2);

  const double u_n = Q[direction + 1] * irho;
  const double c   = std::sqrt(gamma * p * irho);

  std::fill_n(lambda,5,u_n);
  lambda[0] -= c;
  lambda[4] += c;
}

void Euler::EulerSolver_ADERDG::eigenvalues_vect(const double* const * const Q, const int direction, double** const lambda, int size) {
  
}

void Euler::EulerSolver_ADERDG::entropyWave(const double* const x,double t, double* const Q) {
  const double gamma     = 1.4;
  constexpr double width = 1.0;

  #if DIMENSIONS==2
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1]);
  tarch::la::Vector<DIMENSIONS,double> v0(4.5,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(4.5,4.5);
  #else
  tarch::la::Vector<DIMENSIONS,double> xVec(x[0],x[1],x[2]);
  tarch::la::Vector<DIMENSIONS,double> v0(4.5,0.0,0.0);
  tarch::la::Vector<DIMENSIONS,double> x0(4.5,4.5,4.5);
  #endif
  const double distance  = tarch::la::norm2( xVec - x0 - v0 * t );

  Q[0] = 0.5 + 0.3 * std::exp(-distance / std::pow(width, DIMENSIONS));
  Q[1] = Q[0] * v0[0];
  Q[2] = Q[0] * v0[1];
  Q[3] = 0.0;
  // total energy = internal energy + kinetic energy
  const double p = 1.;
  Q[4] = p / (gamma-1)  +  0.5*Q[0] * (v0[0]*v0[0]+v0[1]*v0[1]); // v*v; assumes: v0[2]=0 */
  
}


void Euler::EulerSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt, double* const Q) {
  if (tarch::la::equals(t, 0.0)) {
    entropyWave(x,0.0,Q);
  }
}

exahype::solvers::Solver::RefinementControl
Euler::EulerSolver_ADERDG::refinementCriterion(
    const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
    const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Euler::EulerSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dirichlet conditions
  double Q[NumberOfVariables]     = {0.0};
  double _F[3][NumberOfVariables] = {0.0};
  double* F[3] = {_F[0],_F[1],_F[2]};

  // initialise
  std::fill_n(stateOut, NumberOfVariables, 0.0);
  std::fill_n(fluxOut,  NumberOfVariables, 0.0);
  for (int i=0; i<Order+1; i++) {
    const double ti = t + dt * kernels::gaussLegendreNodes[Order][i];
    entropyWave(x,ti,Q);
    flux(Q,F);
    for (int v=0; v<NumberOfVariables; v++) {
      stateOut[v] += Q[v]            * kernels::gaussLegendreWeights[Order][i];
      fluxOut[v]  += F[direction][v] * kernels::gaussLegendreWeights[Order][i];
    }
  }
  
}

void Euler::EulerSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
    std::copy_n(Q, NumberOfVariables, observables);
}


bool Euler::EulerSolver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {
  if (observablesMin[0] <= 0.0) return false;
  if (observablesMin[4] < 0.0) return false;
  return true;
}
