
#include <mm_malloc.h> //g++

#include "kernels/SWE_MySWESolver_ADERDG/Quadrature.h"

// Use Gauss-Legendre quadrature
double* SWE::MySWESolver_ADERDG_kernels::aderdg::weights1;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::weights2;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::weights3;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::iweights3;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::nodes;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::uh2lob;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::dg2fv;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::fv2dg;

void SWE::MySWESolver_ADERDG_kernels::aderdg::freeQuadratureNodesAndWeights() {
  _mm_free(weights1);
  _mm_free(weights2);
  _mm_free(weights3);
  _mm_free(iweights3);
  _mm_free(nodes);
  _mm_free(uh2lob);
  _mm_free(dg2fv);
  _mm_free(fv2dg);
}

void SWE::MySWESolver_ADERDG_kernels::aderdg::initQuadratureNodesAndWeights() {
  weights1  = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT); //nDofPad
  weights2  = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT); //2D: nDofPad (==weight1), 3D: (nDof*nDof)Pad (== w1[i]*w1[j])
  weights3  = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT); //2D: (nDof*nDof)Pad (== w1[i]*w1[j]), 3D: (nDof*nDof*nDof)Pad (== w1[i]*w1[j]*w1[k])
  iweights3 = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT); //2D: (nDof*nDof)Pad (== w1[i]*w1[j]), 3D: (nDof*nDof*nDof)Pad (== w1[i]*w1[j]*w1[k])
  nodes     = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT);
  uh2lob    = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT); //nDof*nDofPad
  dg2fv     = (double *) _mm_malloc(sizeof(double)*12, ALIGNMENT); //nDof*nDofLimPad
  fv2dg     = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT); //nDofLim*nDofPad
  
  weights1[0]  = 5.000000000000000e-01;
  weights1[1]  = 5.000000000000000e-01;
  weights1[2]  = 0.000000000000000e+00;
  weights1[3]  = 0.000000000000000e+00;

  weights2[0]  = 5.000000000000000e-01;
  weights2[1]  = 5.000000000000000e-01;
  weights2[2]  = 0.000000000000000e+00;
  weights2[3]  = 0.000000000000000e+00;

  weights3[0]  = 2.500000000000000e-01;
  weights3[1]  = 2.500000000000000e-01;
  weights3[2]  = 2.500000000000000e-01;
  weights3[3]  = 2.500000000000000e-01;

  iweights3[0] = 4.000000000000000e+00;
  iweights3[1] = 4.000000000000000e+00;
  iweights3[2] = 4.000000000000000e+00;
  iweights3[3] = 4.000000000000000e+00;

  nodes[0]     = 2.113248654051871e-01;
  nodes[1]     = 7.886751345948129e-01;

  uh2lob[0]    = -3.660254037844387e-01;
  uh2lob[1]    = 1.366025403784439e+00;
  uh2lob[2]    = 0.000000000000000e+00;
  uh2lob[3]    = 0.000000000000000e+00;
  uh2lob[4]    = 1.366025403784439e+00;
  uh2lob[5]    = -3.660254037844387e-01;
  uh2lob[6]    = 0.000000000000000e+00;
  uh2lob[7]    = 0.000000000000000e+00;

  dg2fv[0]    = 1.077350269189626e+00;
  dg2fv[1]    = -7.735026918962576e-02;
  dg2fv[2]    = 0.000000000000000e+00;
  dg2fv[3]    = 0.000000000000000e+00;
  dg2fv[4]    = 5.000000000000001e-01;
  dg2fv[5]    = 4.999999999999999e-01;
  dg2fv[6]    = 0.000000000000000e+00;
  dg2fv[7]    = 0.000000000000000e+00;
  dg2fv[8]    = -7.735026918962576e-02;
  dg2fv[9]    = 1.077350269189626e+00;
  dg2fv[10]    = 0.000000000000000e+00;
  dg2fv[11]    = 0.000000000000000e+00;

  fv2dg[0]    = 7.663460352255524e-01;
  fv2dg[1]    = 3.333333333333333e-01;
  fv2dg[2]    = -9.967936855888598e-02;
  fv2dg[3]    = 0.000000000000000e+00;
  fv2dg[4]    = -9.967936855888587e-02;
  fv2dg[5]    = 3.333333333333333e-01;
  fv2dg[6]    = 7.663460352255527e-01;
  fv2dg[7]    = 0.000000000000000e+00;

}