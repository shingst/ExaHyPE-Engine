#ifndef __INITIAL_DATA_EULERFLOW__
#define __INITIAL_DATA_EULERFLOW__

#include <cstdlib>
#include <iostream>
#include <string>

typedef void (*InitialDataHandler)(const double* const x, double* Q, double t);

extern InitialDataHandler idfunc;

void InitialData(const double* const x, double* Q, double t);

void Polynomial(const double* const x, double* Q, double t);

//ok
void ShuVortex2D(const double* const x, double* Q, double t);

//ok
void MovingGauss(const double* const x, double* V, double t);

void DiffusingGauss(const double* const x, double* Q, double /*t*/);
//Ok ->am spannendsten, Gitter muss fein genug sein
void SmoothedSodShockTube(const double* const x, double* Q, double /*t*/);

#endif /* __INITIAL_DATA_EULERFLOW__ */
