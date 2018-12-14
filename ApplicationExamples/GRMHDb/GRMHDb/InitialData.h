#ifndef __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__
#define __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__


extern "C" {

// FORTRAN functions called by C
void pdesetup_(const int* myrank);


void initialdata_(const double* x, const double* const t, double* Q);


void initialfield_(const double* x, const double* const t, double* Q);

/*
void ShearLayer_(const double* x, const double* const t, double* Q);
void pdelimitervalue_(int* limiter_value, const double* xx,const int* const numberOfObservables, const double* const observablesMin,const double* const observablesMax);
*/


}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
