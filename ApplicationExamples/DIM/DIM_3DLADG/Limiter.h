#ifndef __LIMITER_DATA_ADAPTER_CPP_FORTRAN__
#define __LIMITER_DATA_ADAPTER_CPP_FORTRAN__


extern "C" {

// FORTRAN functions called by C
// void initialdata_(const double* x, const double* const t, double* Q);
void PDEAssurePositivity_(int* status, const double* QMin,const double* QMax);
// only initialdata_ is used, no more.

// void InitialPlaneWave_(const double* x, const double* const t, double* Q);
// void GaussianBubble_(const double* x, const double* const t, double* Q);

// Exact solutions in FORTRAN
//void alfenwave_(const double* x, double* Q, const double* /* scalar */ t);


}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
