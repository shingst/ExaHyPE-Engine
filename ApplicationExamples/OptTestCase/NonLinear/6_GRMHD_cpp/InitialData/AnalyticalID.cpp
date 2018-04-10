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

// The Shocktube solution:

void GRMHD_Shocktube::stateFromParameters::readParameters(const mexa::mexafile& parameters) {
	std::vector<double> vec;
	
	rho = parameters["rho"].as_double();
	press = parameters["press"].as_double();
	
	vec = parameters["vel"].vec(3).as_double();
	DFOR(i) vel.up(i) = vec[i];

	phi = 0;
	vec = parameters["Bmag"].vec(3).as_double();
	DFOR(i) Bmag.up(i) = vec[i];

	alpha = parameters["alpha"].as_double();
	vec = parameters["beta"].vec(3).as_double();
	DFOR(i) beta.up(i) = vec[i];
	
	// Metric is always unity here.
	SYMFOR(i,j) gam.lo(i,j) = SVEC::sym::delta(i,j);
}

std::string GRMHD_Shocktube::stateFromParameters::toString() const {
	// something like this should instead be part of SVEC.
	std::stringstream s2s; // state2string
	s2s << "rho = " << rho << ", ";
	DFOR(i) s2s << "vel.up[" << i << "] = " << vel.up(i) << ", ";
	s2s << "press = " << press << ", ";
	DFOR(i) s2s << "Bmag.up[" << i << "] = " << Bmag.up(i) << ", ";
	s2s << "phi = " << phi << ", ";
	s2s << "alpha = " << alpha << ",";
	DFOR(i) s2s << "beta.up[" << i << "] = " << beta.up(i) << ", ";
	SYMFOR(i,j) s2s << "gam.lo[" << i << "," << j << "] = " << gam.lo(i,j) << ", ";
	return s2s.str();
}

GRMHD_Shocktube::GRMHD_Shocktube() :
	_log("InitialData::GRMHD_Shocktube") {
	// start with reasonable defaults
}

void GRMHD_Shocktube::readParameters(const mexa::mexafile& parameters) {
	prims_right.readParameters(parameters("right"));
	prims_left.readParameters(parameters("left"));
	sep_x = parameters("sep_x").as_double();
	
	logInfo("readParameters()", "Read setup for GMRHD Shocktube at sep_x="<< sep_x);
	logInfo("readParameters()", "State at right: {" << prims_right.toString() << "}");
	logInfo("readParameters()", "State at left: {" << prims_left.toString() << "}");
}


void GRMHD_Shocktube::Interpolate(const double* x, double t, double* Q) {
	bool is_right = x[0] > sep_x;
	GRMHD_Shocktube::state prims_local(is_right ? prims_right : prims_left);
	
	SVEC::GRMHD::Prim2Cons(Q, prims_local.V).copyFullStateVector();
	
	//NVARS(i) printf("Setting ID Q[%d]=%e\n", i, Q[i]);
	
	// if you would only copy prims -> prims.
	// local_state.copy_hydro(Q); // Undefined today, but easy to implement
	// local_state.copy_magneto(Q);
	// local_state.copy_adm(Q);
}

