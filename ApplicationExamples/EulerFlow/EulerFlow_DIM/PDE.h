#ifndef EULERDIM_PDE_H
#define EULERDIM_PDE_H

void initialdata(const double* const x,const double t,double* const Q);
void PDEPrim2Cons(double* Q,double* V);
void PDECons2Prim(double* V, const double* const Q);
void PDEncp(const double* const Q,const double* const gradQ,double* const BgradQ);
void PDEflux(const double* const Q,double** const F);
void PDEEigenvalues(const double* const Q,const int direction,double* const lambda);

#endif
