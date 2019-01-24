#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__


// Fortran functions:
extern "C" {
    void elementcalltecplotplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
	void elementcalltecplotaderdgplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
	void elementcalltecplotfvplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
void finishtecplotplotter_(const int* Myrank);
void initializetecplotplotter_(const double* time);
}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
