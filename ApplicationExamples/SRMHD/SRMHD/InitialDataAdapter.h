#ifndef __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__
#define __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__

#include "tarch/logging/Log.h"
#include <algorithm> // string lower
#include <stdlib.h> // getenv
#include "MHDSolver.h"

#include "exahype/parser/ParserView.h"

extern exahype::parser::ParserView* constants;

extern "C" {

// FORTRAN functions called by C
void initialalfenwave_(double* x, double* Q);
void initialblast_(double* x, double* Q);
void initialorsagtang_(double* x, double* Q);
void initialrotor_(double* x, double* Q);

// Exact solutions in FORTRAN
void alfenwave_(const double* x, double* Q, const double* /* scalar */ t);
	
// C functions called by FORTRAN
void initialdatabyexahypespecfile(double* x, double* Q);

}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
