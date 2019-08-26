// Tools.h  
// Fortran functions:
extern "C" {
//
void pdecons2primfix_(double* uPrim, const double* uCons,int* iErr);
void getnumericalsolution_(double* V,double* Q);
void getexactsolution_(double* x,double* timestep,double* V);
//
}/* extern "C" */
  
