/**
 * Initial data for testing the new C++ GRMHD PDE.
 **/

#include "InitialData.h"

#include <cmath>
constexpr double pi = M_PI;

using namespace SVEC;
using namespace std;


AlfenWave::AlfenWave(const double* const x, const double t, double* Q) : VacuumInitialData(Q) {
	// Computes the AlfenWave variables (Q) at a given time t.
	// Use it ie. with t=0 for initial data
	// Use it for any other time ie. for comparison
	
	/*
	! GRID FOR ALFENWAVE:
	!     dimension const                = 2
	!     width                          = 1.0, 0.3
	!     offset                         = 0.0, 0.0
	!     end-time                       = 2.1
	!
	!  maximum-mesh-size              = 0.04
	*/

	constexpr double time_offset = 1.0;
	constexpr double gamma = SVEC::GRMHD::Parameters::gamma;

	double eta  = 1.;
	double B0   = 1.0;
	double hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0;
	double tempaa = rho0 * hh + B0*B0 * ( 1.0 + eta*eta);
	double tempab = 2.0 * eta * B0*B0 / tempaa;
	double tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab*tempab));
	double va2 = B0*B0 / ( tempaa * tempac);
	double vax = sqrt(va2);
	
	// as AlfenWave is a VacuumInitialData, all initial data are already set,
	// especially rho and press are already set to a value.

	// c2p-invariant: Magnetic field
	Bmag.up(0) = B0;
	Bmag.up(1) = eta * B0 * cos(2*pi*( x[0] - vax*(t-time_offset)));
	if(TDIM>2)
	Bmag.up(2) = eta * B0 * sin(2*pi*( x[0] - vax*(t-time_offset)));

	// primitive variables
	DFOR(i) vel.up(i) = - vax * Bmag.up(i) / B0;
	vel.up(0) = 0.0;
	
	// if ExaHyPE dimensiosn are 2. Remember: SVEC dimensions are independent and
	// these AlfenWave ID MUST be computed in 3D.
	// if(DIMENSIONS == 2 && TDIM == 3) { Bmag.up(2) = vel.up(2) = 0; } // this is WRONG
	//Prim2Cons(Q,V);
	
	// instead of doing p2c, just copy everything over in order to produce
	// primitive output variables in Q.
	
	Dens = rho;
	tau = press;
	DFOR(i) Si.lo(i) = vel.up(i);
} // AlfenWave


