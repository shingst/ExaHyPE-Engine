#ifndef __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__
#define __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__


extern "C" {

// FORTRAN functions called by C
void initialdata_(const double* x, const double* const t, double* Q);

// only initialdata_ is used, no more.

void ShearLayer_(const double* x, const double* const t, double* Q);
void pdelimitervalue_(int* limiter_value, const double* xx);



}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
