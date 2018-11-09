#ifndef __INITIAL_DATA_ADAPTER_GRMHD__
#define __INITIAL_DATA_ADAPTER_GRMHD__

#include <string>

extern "C" {

// FORTRAN functions called by C
void initialdata_(const double* x, const double* const t, double* Q);

void initialaccretiondisc_(const double* x, const double* const t,  double* Q);
void initialaccretiondisc3d_(const double* x, const double* const t, double* Q);

}/* extern "C" */

#endif /* __INITIAL_DATA_ADAPTER_GRMHD__ */
