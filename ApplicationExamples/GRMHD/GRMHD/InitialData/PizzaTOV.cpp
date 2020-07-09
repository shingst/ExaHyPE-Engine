
#include "tarch/logging/Log.h"
#include "AbstractGRMHDSolver_ADERDG.h"
#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef PIZZATOV_AVAILABLE
/*
pizzatov::pizzatov() {
	printf("Cannot call Pizza as not compiled with -DPIZZATOV_AVAILABLE");
	abort();
}
void pizzatov::Interpolate(const double* x, double t, double* Q) {}

#else

#include "pizza_tovfront/pointwise_tov.h"

using namespace Pizza;
using namespace Pizza::TOV;

constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;


struct ExaHyPE_Pizza_Logger : public Pizza::TOV::logger {
	tarch::logging::Log _log;
	ExaHyPE_Pizza_Logger() : _log("GRMHD::PizzaTOV") {}
	void info(const std::string& msg, const std::string& codepos="") override { logInfo("", msg); }
	void warn(const std::string& msg, const std::string& codepos="") override { logWarning("", msg); }
};

pizzatov::pizzatov() {
	tov = new pointwise_tov();
	tov->log = new ExaHyPE_Pizza_Logger();
	
	// here is the place to set some TOV parameters,
	// if needed. Such as the EOS.
	
	std::string EOS_BU =
		"name      = BU (Polytropic)\n"
		"type      = polytrope\n"
		"poly_n    = 1\n"
		"poly_rmd  = 6.176e+18\n"
	;
	tov->read_eos_from_string(EOS_BU);
	
	// set mapping to GRMHD application
	// rho:1,vel:3,E:1,B:3,psi:1,lapse:1,shift:3,gij:6
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts vars;

	// GRMHD Material parameters
	tov->adm.alp = vars.lapse;
	tov->adm.shift1 = vars.shift + 0;
	tov->adm.shift2 = vars.shift + 0;
	tov->adm.shift3 = vars.shift + 0;
	tov->adm.gxx = vars.gij + 0;
	tov->adm.gxy = vars.gij + 1;
	tov->adm.gxz = vars.gij + 2;
	tov->adm.gyy = vars.gij + 3;
	tov->adm.gyz = vars.gij + 4;
	tov->adm.gzz = vars.gij + 5;
	
	// GRMHD actual computation variables
	// Note that PizzaTOV sets primitives and we have to do the cons2prim afterwards.
	tov->hydro.rho = vars.rho;
	tov->hydro.velx = vars.vel + 0;
	tov->hydro.vely = vars.vel + 1;
	tov->hydro.velz = vars.vel + 2;
	tov->hydro.press = vars.E;
	
	// do the main computation (takes no time)
	tov->compute_star();
}

void pizzatov::Interpolate(const double* x, double t, double* Q) {
	double V[nVar] = {0.0}; // IMPORTANT: Zero V before data go to PizzaTOV.
	#if DIMENSIONS == 2
	// this is only useful for lower dimensional debugging of the
	// ExaHyPE infrastructure and not suitable Initial Data for simulation.
	// If you want to do 2D, Pizza can use polar coordinates but it hasn't
	// been implemented in the pizza_tovfront.
	
	// We now abuse this to get only values on the x axis, which effectively
	// gives us a radius coordinate. This is fine for a nonrotating star.
	double x3d[3] = { x[0], 0.0, 0.0 };
	tov->initial_data(x3d, V);
	#else
	tov->initial_data(x, V);
	#endif
	
	// treatment of the atmostphere not done by PizzaTOV
	const double atmo_rho = 1e-13;
	const double atmo_press = 1e-7;
	if(V[0] < atmo_rho) {
		V[0] = atmo_rho;
	}
	if(V[4] < atmo_press) {
		V[4] = atmo_press;
	}
	
	pdeprim2cons_(Q, V);
}
*/
#endif  /* PIZZATOV_AVAILABLE */
