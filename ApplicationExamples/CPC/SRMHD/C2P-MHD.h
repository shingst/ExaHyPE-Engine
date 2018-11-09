#ifndef __C2P_MHD_CPP__
#define __C2P_MHD_CPP__

/**
 * C definitions to use the MHD Con2Prim functions
 **/

extern "C" {

// FORTRAN functions called by C
void pdeprim2cons_(double* /* OUT */ Q, double* /* IN */ V);
void pdecons2prim_(double* /* OUT */ V, double* /* IN */ Q, int* iErr);


} /* extern "C" */
#endif /* __C2P_MHD_CPP__ */