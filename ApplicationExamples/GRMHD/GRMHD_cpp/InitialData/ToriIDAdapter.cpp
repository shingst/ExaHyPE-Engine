#include "ToriIDAdapter.h"

// Calling Fortran...

extern "C" {
	void compute_disc_();
	void initialtorusonkerrschild_(const double* const, double* Q);

	// Members of the MODULE Parameters_Disc.
	// By the ISO C interface, they are exposed as global variables in no
	// particular order
	
	extern double
		Parameters_Disc_el0,
		Parameters_Disc_delta,
		Parameters_Disc_eos_K,
		Parameters_Disc_eos_gamma,
		Parameters_Disc_Mass_grav,
		Parameters_Disc_Mbh,
		Parameters_Disc_aom,
		Parameters_Disc_verbose,
		Parameters_Disc_debug;
}

void ToriIDAdapter::prepare() {
	logInfo("prepare()", "Setting up disc (Fortran output will follow)");
	// Parameters_Disc_verbose = true;

	Parameters_Disc_Mass_grav = 1.0;
	Parameters_Disc_el0 = 3.8;
	Parameters_Disc_delta = -0.7;
	
	compute_disc_();
}

void ToriIDAdapter::Interpolate(const double* x, double t, double* Q) {
	// ToriID provide already an conserved vector
	initialtorusonkerrschild_(x, Q);
}
