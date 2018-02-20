#ifndef __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__
#define __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__

#include "tarch/logging/Log.h"
#include <algorithm> // string lower
#include <stdlib.h> // getenv

#include "exahype/parser/ParserView.h"

extern "C" {

// FORTRAN functions called by C
void initialjet_(const double* const x, double* Q);
void initialalfenwave_(const double* const x, double* Q);
void initialblast_(const double* const x, double* Q);
void initialorsagtang_(const double* const x, double* Q);
void initialrotor_(const double* const x, double* Q);
void initialshocktube_(const double* const x, double* Q);

typedef void (*InitialDataHandler)(const double* const x, double* Q);
extern InitialDataHandler idfunc;

// Exact solutions in FORTRAN
void alfenwave_(const double* x, double* Q, const double* /* scalar */ t);
	





}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
