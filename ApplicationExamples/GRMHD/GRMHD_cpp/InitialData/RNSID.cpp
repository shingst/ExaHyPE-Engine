#include "AbstractGRMHDSolver_ADERDG.h"
#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "PDE/PDE.h"

using SVEC::GRMHD::Prim2Cons;
#include "PDE/tensish.cpph"
using namespace tensish;

#ifndef RNSID_AVAILABLE

#include <stdlib.h>
#include <stdio.h>

rnsid::rnsid() {
        printf("Cannot call RNSID as not compiled with -DRNSID_AVAILABLE");
        abort();
}

void rnsid::Interpolate(const double* x, double t, double* Q) {}

#else /* RNSID_AVAILABLE */

// A small Fortran helper routine
extern "C" void rnsid_2d_cartesian2cylindrical_(double*Cartesian, double*Cylindrical, double*xGP_cylindrical);

#include "rnsid/rnsid.h"
#include <array>

rnsid::rnsid() {
	id = new RNSID::rnsid();
	
	// A TOV star
	id->axes_ratio = 1.0;
	id->rnsid_rho_min = 1e-10;
	
	// mapping for quantity vector
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts var;
	
	// mind the order of Gij in the particular GRMHD application! (F vs C)
	id->adm_idx.gxx = var.gij + sym::index(0,0);
	id->adm_idx.gxy = var.gij + sym::index(0,1);
	id->adm_idx.gxz = var.gij + sym::index(0,2);
	id->adm_idx.gyy = var.gij + sym::index(1,1);
	id->adm_idx.gyz = var.gij + sym::index(1,2);
	id->adm_idx.gzz = var.gij + sym::index(2,2);
	
	id->adm_idx.alp    = var.lapse;
	id->adm_idx.shift1 = var.shift + 0;
	id->adm_idx.shift2 = var.shift + 1;
	id->adm_idx.shift3 = var.shift + 2;
	
	// Note these are the primitives? ...
	id->hydro_idx.rho   = 0;
	id->hydro_idx.velx  = var.vel + 0; // 1
	id->hydro_idx.vely  = var.vel + 1; // 2
	id->hydro_idx.velz  = var.vel + 2; // 3
	id->hydro_idx.press = 4; // we store the pressure
	
	id->Run();
}

void rnsid::readParameters(const mexa::mexafile& para) {
	/**
	 * Attention, here we don't read all parameters but only
	 * some of them. Currently, we also don't support default
	 * values but require all of them to be set. Which is not
	 * too bad if you have proper parameter files.
	 **/
	
	id->rho_center = para["rho_center"].as_double();
	id->log_enth_center = para["log_enth_center"].as_double();
	id->rho_cut = para["rho_cut"].as_double();
	id->rnsid_rho_min = para["rho_min"].as_double();
	
	id->axes_ratio = para["axes_ratio"].as_double();
	id->accuracy = para["accuracy"].as_double();
	id->perturbation = para["perturbation"].get_bool();
	
	id->rotation_type = para["rotation_type"].get_string();
	id->A_diff = para["A_diff"].get_double();
	
	id->zero_shift = para["zero_shift"].get_bool();
}

void rnsid::Interpolate(const double* pos, double t, double* Q) {
	constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
	double V[nVar] = {0.0}; // primitive variables, as returned by rnsid
	
	if(DIMENSIONS == 2) {
		// in 2D, we interpret the coordinates (x,y) as (r,z)
		// Transfrom from cylindrical coordinates to cartesian ones:
		const double r=pos[0], z=pos[1], phi = 0;

		double pos_cylindrical[3] = { r, z, phi };
		double pos_cartesian[3]   = { r*std::cos(phi), r*std::sin(phi), z  };
		
		double  V_Cartesian[nVar];  // intermediate step from rnsid
		double *V_Cylindrical = V; // cylindrical output
		
		id->Interpolate(pos_cartesian, V_Cartesian);

		// Transform the metric from Cartesian back to cylindrical coordinates.
		// ALERT: This is suitable for the GRMHD part but does not cover the curvature
		//        part (CCZ4). It is ignored here anyway.

		rnsid_2d_cartesian2cylindrical_(V_Cartesian, V_Cylindrical, pos_cartesian);
	} else {
		id->Interpolate(pos, V);
	}
	
	// treatment of the atmostphere PROBABLY not done by RNSID
	const double atmo_rho = SVEC::GRMHD::Parameters::atmo_rho;
	const double atmo_press = SVEC::GRMHD::Parameters::atmo_press;
	if(V[0] < atmo_rho) {
		V[0] = atmo_rho;
	}
	if(V[4] < atmo_press) {
		V[4] = atmo_press;
	}
	
	//NVARS(i) printf("V[%d]=%e\n", i, V[i]);
	for(int i=0;i<nVar;i++) Q[i] = 0.0;
	
	Prim2Cons p2c(Q, V);
	//NVARS(i) printf("Q[%d]=%e\n", i, Q[i]);
	p2c.copyFullStateVector();
}

#endif /* RNSID_AVAILABLE */
