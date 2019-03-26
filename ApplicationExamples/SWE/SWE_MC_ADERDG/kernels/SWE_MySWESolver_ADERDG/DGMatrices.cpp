
#include <mm_malloc.h> //g++

#include "kernels/SWE_MySWESolver_ADERDG/DGMatrices.h"

//DGMatrices
double* SWE::MySWESolver_ADERDG_kernels::aderdg::Kxi;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::Kxi_T;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::iK1_T;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::dudx;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::dudx_T;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::FLCoeff;
double* SWE::MySWESolver_ADERDG_kernels::aderdg::FRCoeff;
double** SWE::MySWESolver_ADERDG_kernels::aderdg::fineGridProjector1d;
double** SWE::MySWESolver_ADERDG_kernels::aderdg::fineGridProjector1d_T_weighted;


void SWE::MySWESolver_ADERDG_kernels::aderdg::freeDGMatrices() {
  _mm_free(FLCoeff);
  _mm_free(FRCoeff);
  _mm_free(dudx);
  _mm_free(dudx_T);
  _mm_free(iK1_T);
  _mm_free(Kxi);
  _mm_free(Kxi_T);
  
  _mm_free(fineGridProjector1d[0]);
  _mm_free(fineGridProjector1d[1]);
  _mm_free(fineGridProjector1d[2]);
  delete [] fineGridProjector1d;
  
  _mm_free(fineGridProjector1d_T_weighted[0]);
  _mm_free(fineGridProjector1d_T_weighted[1]);
  _mm_free(fineGridProjector1d_T_weighted[2]);
  delete [] fineGridProjector1d_T_weighted;
}


void SWE::MySWESolver_ADERDG_kernels::aderdg::initDGMatrices() {
  
  FLCoeff = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT);
  FRCoeff = (double *) _mm_malloc(sizeof(double)*4, ALIGNMENT);
  //note: FLCoeff is also F0
  
  dudx    = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);
  dudx_T  = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);
  iK1_T   = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);
  Kxi     = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);
  Kxi_T   = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);

  fineGridProjector1d            = new double* [3];
  fineGridProjector1d_T_weighted = new double* [3];
  for(int i=0; i<3; i++) {
    fineGridProjector1d[i]            = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);
    fineGridProjector1d_T_weighted[i] = (double *) _mm_malloc(sizeof(double)*8, ALIGNMENT);
  }
  
  FLCoeff[0] = 1.366025403784439e+00;
  FLCoeff[1] = -3.660254037844387e-01;
  FLCoeff[2] = 0.000000000000000e+00;
  FLCoeff[3] = 0.000000000000000e+00;

  FRCoeff[0] = -3.660254037844387e-01;
  FRCoeff[1] = 1.366025403784439e+00;
  FRCoeff[2] = 0.000000000000000e+00;
  FRCoeff[3] = 0.000000000000000e+00;

  dudx[0] = -1.732050807568877e+00;
  dudx[1] = -1.732050807568877e+00;
  dudx[2] = 0.000000000000000e+00;
  dudx[3] = 0.000000000000000e+00;
  dudx[4] = 1.732050807568877e+00;
  dudx[5] = 1.732050807568877e+00;
  dudx[6] = 0.000000000000000e+00;
  dudx[7] = 0.000000000000000e+00;

  dudx_T[0] = -1.732050807568877e+00;
  dudx_T[1] = 1.732050807568877e+00;
  dudx_T[2] = 0.000000000000000e+00;
  dudx_T[3] = 0.000000000000000e+00;
  dudx_T[4] = -1.732050807568877e+00;
  dudx_T[5] = 1.732050807568877e+00;
  dudx_T[6] = 0.000000000000000e+00;
  dudx_T[7] = 0.000000000000000e+00;

  iK1_T[0] = 6.666666666666667e-01;
  iK1_T[1] = -2.440169358562924e-01;
  iK1_T[2] = 0.000000000000000e+00;
  iK1_T[3] = 0.000000000000000e+00;
  iK1_T[4] = 9.106836025229592e-01;
  iK1_T[5] = 6.666666666666666e-01;
  iK1_T[6] = 0.000000000000000e+00;
  iK1_T[7] = 0.000000000000000e+00;

  Kxi[0] = -8.660254037844387e-01;
  Kxi[1] = 8.660254037844387e-01;
  Kxi[2] = 0.000000000000000e+00;
  Kxi[3] = 0.000000000000000e+00;
  Kxi[4] = -8.660254037844387e-01;
  Kxi[5] = 8.660254037844387e-01;
  Kxi[6] = 0.000000000000000e+00;
  Kxi[7] = 0.000000000000000e+00;

  Kxi_T[0] = -8.660254037844387e-01;
  Kxi_T[1] = -8.660254037844387e-01;
  Kxi_T[2] = 0.000000000000000e+00;
  Kxi_T[3] = 0.000000000000000e+00;
  Kxi_T[4] = 8.660254037844387e-01;
  Kxi_T[5] = 8.660254037844387e-01;
  Kxi_T[6] = 0.000000000000000e+00;
  Kxi_T[7] = 0.000000000000000e+00;

  fineGridProjector1d[0][0] = 1.244016935856292e+00;
  fineGridProjector1d[0][1] = -2.440169358562924e-01;
  fineGridProjector1d[0][2] = 0.000000000000000e+00;
  fineGridProjector1d[0][3] = 0.000000000000000e+00;
  fineGridProjector1d[0][4] = 9.106836025229591e-01;
  fineGridProjector1d[0][5] = 8.931639747704091e-02;
  fineGridProjector1d[0][6] = 0.000000000000000e+00;
  fineGridProjector1d[0][7] = 0.000000000000000e+00;
  fineGridProjector1d[1][0] = 6.666666666666666e-01;
  fineGridProjector1d[1][1] = 3.333333333333334e-01;
  fineGridProjector1d[1][2] = 0.000000000000000e+00;
  fineGridProjector1d[1][3] = 0.000000000000000e+00;
  fineGridProjector1d[1][4] = 3.333333333333333e-01;
  fineGridProjector1d[1][5] = 6.666666666666667e-01;
  fineGridProjector1d[1][6] = 0.000000000000000e+00;
  fineGridProjector1d[1][7] = 0.000000000000000e+00;
  fineGridProjector1d[2][0] = 8.931639747704107e-02;
  fineGridProjector1d[2][1] = 9.106836025229590e-01;
  fineGridProjector1d[2][2] = 0.000000000000000e+00;
  fineGridProjector1d[2][3] = 0.000000000000000e+00;
  fineGridProjector1d[2][4] = -2.440169358562926e-01;
  fineGridProjector1d[2][5] = 1.244016935856293e+00;
  fineGridProjector1d[2][6] = 0.000000000000000e+00;
  fineGridProjector1d[2][7] = 0.000000000000000e+00;

  fineGridProjector1d_T_weighted[0][0] = 4.146723119520975e-01;
  fineGridProjector1d_T_weighted[0][1] = 3.035612008409864e-01;
  fineGridProjector1d_T_weighted[0][2] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[0][3] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[0][4] = -8.133897861876413e-02;
  fineGridProjector1d_T_weighted[0][5] = 2.977213249234697e-02;
  fineGridProjector1d_T_weighted[0][6] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[0][7] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[1][0] = 2.222222222222222e-01;
  fineGridProjector1d_T_weighted[1][1] = 1.111111111111111e-01;
  fineGridProjector1d_T_weighted[1][2] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[1][3] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[1][4] = 1.111111111111111e-01;
  fineGridProjector1d_T_weighted[1][5] = 2.222222222222222e-01;
  fineGridProjector1d_T_weighted[1][6] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[1][7] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[2][0] = 2.977213249234702e-02;
  fineGridProjector1d_T_weighted[2][1] = -8.133897861876419e-02;
  fineGridProjector1d_T_weighted[2][2] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[2][3] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[2][4] = 3.035612008409863e-01;
  fineGridProjector1d_T_weighted[2][5] = 4.146723119520975e-01;
  fineGridProjector1d_T_weighted[2][6] = 0.000000000000000e+00;
  fineGridProjector1d_T_weighted[2][7] = 0.000000000000000e+00;

}