#include "AbstractGRMHDSolver_ADERDG.h"
#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#ifndef RNSID_AVAILABLE

#include <stdlib.h>
#include <stdio.h>

rnsid::rnsid() {
        printf("Cannot call RNSID as not compiled with -DRNSID_AVAILABLE");
        abort();
}

void rnsid::Interpolate(const double* x, double t, double* Q) {}

#else /* RNSID_AVAILABLE */

#include "rnsid/rnsid.h"

rnsid::rnsid() {
	id = new RNSID::rnsid();
	
	// A TOV star
	id->axes_ratio = 1.0;
	id->rnsid_rho_min = 1e-10;
	
	// mapping for quantity vector
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts var;
	id->adm_idx.gxx = var.gij + 0;
	id->adm_idx.gxy = var.gij + 1;
	id->adm_idx.gxz = var.gij + 2;
	id->adm_idx.gyy = var.gij + 3;
	id->adm_idx.gyz = var.gij + 4;
	id->adm_idx.gzz = var.gij + 5;
	
	id->adm_idx.alp    = var.lapse;
	id->adm_idx.shift1 = var.shift + 0;
	id->adm_idx.shift2 = var.shift + 1;
	id->adm_idx.shift3 = var.shift + 2;
	
	// Note these are the primitives? ...
	id->hydro_idx.rho   = var.rho;
	id->hydro_idx.velx  = var.vel + 0;
	id->hydro_idx.vely  = var.vel + 1;
	id->hydro_idx.velz  = var.vel + 2;
	id->hydro_idx.press = var.E;
	
	id->Run();
}

void get_conserved_quantities(rnsid* rnsid, const double* pos, double* Q) {
	constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;

	double V[nVar] = {0.0}; // primitive variables, as returned by rnsid

	rnsid->id->Interpolate(pos, V);
	
	// treatment of the atmostphere PROBABLY not done by RNSID
	const double atmo_rho = 1e-13;
	const double atmo_press = 1e-7;
	if(V[0] < atmo_rho) {
		V[0] = atmo_rho;
	}
	if(V[4] < atmo_press) {
		V[4] = atmo_press;
	}
	
	//NVARS(i) printf("V[%d]=%e\n", i, V[i]);
	for(int i=0;i<nVar;i++) Q[i] = 0.0;
	
	pdeprim2cons_(Q, V);
}


// A small Fortran helper routine
extern "C" void rnsid_2d_cartesian2cylindrical_(double* Cartesian, const double* const Cylindrical, const double* const xGP_cylindrical);

	
void rnsid::Interpolate(const double* pos, double t, double* Q) {
	constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;

	// Copied from GRMHD_cpp, rnsid::Interpolate:
	if(DIMENSIONS == 2) {
		// in 2D, we interpret the coordinates (x,y) as (rho,z)
		// Transfrom from cylindrical coordinates to cartesian ones:
		const double rho=pos[0], z=pos[1], phi = 0.25*M_PI;

		// The source and target coordinate systems. Mind that they have to be right-oriented
		// coordinate systems for the cross product of Bmag and Scon.
		double pos_cylindrical[3] = { rho, phi, z };
		double pos_cartesian[3]   = { rho*std::cos(phi), rho*std::sin(phi), z };
		
		if(rho==0) {
			throw std::domain_error("RNSID cannot be evaluated in 2D at the axis (radius=0) because coordinates are singular.");
		}
		
		double  Q_Cartesian[nVar] = {0.0};  // intermediate step from rnsid
		double *Q_Cylindrical = Q; // cylindrical output
		
		get_conserved_quantities(this, pos_cartesian, Q_Cartesian);
		
		// Since I don't believe this translation stuff, let me setup something on my own
		// here.
		//V_Cartesian[0] = std::sqrt( pos_cartesian[0]*pos_cartesian[0] + pos_cartesian[1]*pos_cartesian[1] + pos_cartesian[2]*pos_cartesian[2] );
		//V_Cartesian[1] = 
		

		// Transform the metric from Cartesian back to cylindrical coordinates.
		rnsid_2d_cartesian2cylindrical_(Q_Cylindrical, Q_Cartesian, pos_cylindrical);
	} else {
		get_conserved_quantities(this, pos, Q);
	}
	
	/*
	double V[nVar] = {0.0};
	#if DIMENSIONS == 2
	// this is only useful for lower dimensional debugging of the
	// ExaHyPE infrastructure and not suitable Initial Data for simulation.
	// If you want to do 2D, Pizza can use polar coordinates but it hasn't
	// been implemented in the pizza_tovfront.
	
	// We now abuse this to get only values on the x axis, which effectively
	// gives us a radius coordinate. This is fine for a nonrotating star.
	double x3d[3] = { x[0], 0.0, 0.0 };
	id->Interpolate(x3d, V);
	#else
	id->Interpolate(x, V);
	#endif
	*/
	

}

#endif /* RNSID_AVAILABLE */
