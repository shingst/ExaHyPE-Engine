#include "InitialData.h"

/* GUARD to compile only if the TwoPunctures code is available */
#ifdef TWOPUNCTURES_AVAILABLE

// TP Black Holes Initial Data code
#include "libtwopunctures/TwoPunctures.h"

// as always...
#include "AbstractCCZ4Solver_ADERDG.h"
typedef CCZ4::AbstractCCZ4Solver_ADERDG::Variables Variables;

#include "kernels/KernelUtils.h" // idx classes
#include "kernels/GaussLegendreBasis.h" // gaussLegendreNodes
using namespace kernels;

#include "tarch/la/Vector.h"
using namespace tarch::la;
typedef Vector<DIMENSIONS,double> dvec;
typedef Vector<ImportedTwoPunctures::numberOfVariables,double> Qvec;

#include "CCZ4Solver_ADERDG_Variables.h" // subvectors and so

#include "kernels/aderdg/generic/c/computeGradients.cpph" // computeGradQ
#include "Helpers/tensors.h" // Tensor math

#include <algorithm>  // std::fill_n

#include "Helpers/printonce.h" // debugging

#include <cmath> // std::pow etc.

// Tie TwoPunctures output to Tarch Logging
struct ExaHyPE_TwoPunctures_Logger : public TP::logger {
	tarch::logging::Log _log;
	ExaHyPE_TwoPunctures_Logger() : _log("GRMHD::PizzaTOV") {}

	void log(const std::string& msg) override { logInfo("", msg); }
	void error(const std::string& msg) override { logError("", msg); }
	void info(const std::string& msg) override { logInfo("", msg); }
	void warn(const std::string& msg) override { logWarning("", msg); }
};

ImportedTwoPunctures::ImportedTwoPunctures() {
	tp = new TP::TwoPunctures(); // Pointer for lesser compil dependency
	tp->log = new ExaHyPE_TwoPunctures_Logger;

	// parameters similar to CCZ4-Vakuum-ET TwoPunctures
	// cf LOEWE /home/astro/koeppel/parfiles/CCZ4-Vakuum-ET
	tp->par_b             =  1;

	const dvec BHPOS(0,0,0);
	
	tp->target_M_plus     =  1.0;
	tp->par_P_plus[0]     =  BHPOS[0];
	tp->par_P_plus[1]     =  BHPOS[1];
	tp->par_P_plus[2]     =  BHPOS[2];
	tp->par_S_plus[0]     =  0.0;
	tp->par_S_plus[1]     =  0.0;
	tp->par_S_plus[2]     =  0.0;
	tp->center_offset[0] =  -1;

	tp->target_M_minus    =  0.0;
	tp->par_P_minus[0]    =  0.0;
	tp->par_P_minus[1]    =  0.0;
	tp->par_P_minus[2]    =  0.0;
	tp->par_S_minus[0]    =  0.0;
	tp->par_S_minus[1]    =  0.0;
	tp->par_S_minus[2]    =  0.0;

	tp->grid_setup_method = "evaluation";
	// smoothen out the infinities at punctures
	tp->TP_epsilon = 1e-6;
}

void ImportedTwoPunctures::prepare() {
	static printonce msg("ImportedTwoPunctures", "Preperation run started");
	tp->Run();
}

void ImportedTwoPunctures::Interpolate(const double* const x, double* Q) {
	Variables V(Q);
	std::fill_n(Q, numberOfVariables, 0.0);
	
	// TwoPunctures sets us only gij, kij, alpha.
	tp->Interpolate(x, Q);
	
	// do a check if we are at the position of the BH
	// if( abs(x - BHPOS) < 1e-9 ) { Q() = 0; }
	
	// An arbitary but well defined initial value for beta = 0
	// Z, Theta, b all zero.
	// dLapse, dShift, dG all zero.
	
	// Compute trace of K = K_{ij} gamma^{ij}
	using namespace TP::Z4VectorShortcuts;
	const Tensors::metric m(Q + g11);
	const Tensors::mats_u K(Q + K11);
	
	V.traceK() = m.trace(K);
	V.K0() = V.traceK();
	V.phi() = std::pow(Tensors::det(m.lo), (-1./6.));
	double phi2 = V.phi()*V.phi();
	
	// ABORT ABORT ABORT HERE
	
	/*
	// set K_{ij} := \tilde A_{ij}, the conformal traceless part of the extr. curvature
	for(int ij=0; ij<6; ij++) {
		V.K(ij) = phi2*( V.K(ij) - 1./3.*V.traceK()*V.G(ij) );
	}
	
	// set g_{ij} := \tilde g_{ij} = phi^2 * g_{ij}
	for(int ij=0; ij<6; ij++) {
		V.G(ij) = phi2*V.G(ij);
	}
	
	// Z_i = \tilde Gamma^i = \tilde gamma^{ij} \tilde Gamma^i_{jk}
	const Tensors::metric gtilde(Q + g11);
	kernels::dshadow DD( G......, 3, 3, 3);
	for(int i=0; i<3; i++) {
	for(int j=0; j<3; i++) {
	for(int k=0; k<3; i++) {
	for(int l=0; l<3; i++) {
		V.Z(i) = V.Z(i) + 1./phi2 * (gtilde.lo(i,j)*gtilde.up(k,l)*(2 * DD(l,j,k) + 2*V.P(l) * gtilde.lo(j,k) ) );
	}}}}
	*/


	// in-place Prim2Cons:
	V.lapse() = std::exp(V.lapse());
	V.phi()   = std::exp(V.phi());
}

void ImportedTwoPunctures::ComputeAuxillaries(const double* const x, double* Q, double *Qx, double *Qy, double *Qz) {
	// Compute pure derivatives to auxilliary variables
	
	// We now have G, K and the lapse
	// We have to compute the auxiliary variables and the shiftmake
	// We also must ensure that NaNs are cured which may appear
	// when the interpolation point is exactly at a BH position.

	Variables V(Q), Vx(Qx), Vy(Qy), Vz(Qz);
	
	// Ai = (\partial_i alpha) / alpha
	V.dLapse(0) = Vx.lapse() / V.lapse();
	V.dLapse(1) = Vy.lapse() / V.lapse();
	V.dLapse(2) = Vz.lapse() / V.lapse();
	
	// Bki = \partial_k beta^i
	V.dxShift( Vx.shift() );
	V.dyShift( Vy.shift() );
	V.dzShift( Vz.shift() );
	
	// D in Z4: D_kij = 0.5 \partial_k g_ij
	// D in CCZ4::  D_kij = 0.5 \partial i \tilde \gamma_jk 
	//      = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
	// this is wrong TODO:
	V.dxG( 0.5 * Vx.G() );
	V.dyG( 0.5 * Vy.G() );
	V.dzG( 0.5 * Vz.G() );
	
	// P in CCZ4: P_i = \partial_i \phi
	V.P(0) = Vx.phi();
	V.P(1) = Vx.phi();
	V.P(2) = Vx.phi();
}

void ImportedTwoPunctures::dInterpolate(const double* const x, double *Qx, double *Qy, double *Qz) {
	const double epsilon = 1e-7, eps4 = 1e-4;
	double Qp1[numberOfVariables], Qm1[numberOfVariables], Qp2[numberOfVariables], Qm2[numberOfVariables];
	//Qvec Qp1, Qm1, Qp2, Qm2 // would be so nice if there was a modifiable Qd.data()
	dvec xc(x[0],x[1],x[2]); // fucking workaround for missing tarch::la::vector(const double* const) constructor.
	
	// Metric derivative computed with a fourth order central finite difference 
	// as done by Michael D
	double *Qd;
	for(int d=0; d<3; d++) {
		dvec xp1(xc), xm1(xc), xp2(xc), xm2(xc);
		xp1(d) += eps4;
		xm1(d) -= eps4;
		xp2(d) += 2*eps4;
		xm2(d) -= 2*eps4;
		
		Interpolate(xp1.data(), Qp1);
		Interpolate(xm1.data(), Qm1);
		Interpolate(xp2.data(), Qp2);
		Interpolate(xm2.data(), Qm2);
		
		if(d==0) Qd = Qx;
		if(d==1) Qd = Qy;
		if(d==2) Qd = Qz;
		
		for(int i=0; i<numberOfVariables; i++) {
			Qd[i] = ( 8.0*Qp1[i] - 8.0*Qm1[i]  + Qm2[i]   - Qp2[i]  )/(12.0*eps4);
		}
	}
}

void ImportedTwoPunctures::adjustPointSolution(const double* const x,const double t,double* Q) {
	// Compute the TwoPunctures with Finite Difference aux variables
	static printonce msg("ImportedTwoPunctures", "Pointwise adjustPointSolution");
	
	Interpolate(x, Q);
	
	// determine gradients with Finite Difference
	double Qx[numberOfVariables], Qy[numberOfVariables], Qz[numberOfVariables];
	dInterpolate(x, Qx, Qy, Qz);
	
	ComputeAuxillaries(x, Q, Qx, Qy, Qz);
}

void ImportedTwoPunctures::adjustPatchSolution(
		const tarch::la::Vector<DIMENSIONS, double>& center,
		const tarch::la::Vector<DIMENSIONS, double>& dx,
		const double t, const double dt, double* luh) {
	static printonce msg("ImportedTwoPunctures", "Patchwise solutionAdjustment");
	
	dvec r, q;
	double *Q, *Qx, *Qy, *Qz;
	double gradQ[basisSize3 * DIMENSIONS * numberOfVariables];
	idx5 idx_gradQ(basisSize, basisSize, basisSize, DIMENSIONS, numberOfVariables);
	idx4 idx_luh(basisSize, basisSize, basisSize, numberOfVariables);
	
	// First, set the incomplete ID on all DOF
	for(int iz = 0; iz < basisSize; iz++) {
		q[2] = kernels::legendre::nodes[order][iz];
		for(int iy = 0; iy < basisSize; iy++) {
			q[1] = kernels::legendre::nodes[order][iy];
			for(int ix = 0; ix < basisSize; ix++) {	
				q[0] = kernels::legendre::nodes[order][ix];
				
				r  = center + dx * ( q - 0.5 );
				Q  = &luh[idx_luh(iz,iy,ix,0)];
				Interpolate(r.data(), Q);
			} // x
		} // y
	} // z
	
	// Second, for the time being, just compute *all* gradients.
	// Initial Data are computed only once. So it can be slow.
        kernels::aderdg::generic::c::computeGradQ<CCZ4::AbstractCCZ4Solver_ADERDG>(gradQ, luh, dx);
	
	// Third, compute the gradients
	for(int iz = 0; iz < basisSize; iz++) {
		q[2] = kernels::legendre::nodes[order][iz];
		for(int iy = 0; iy < basisSize; iy++) {
			q[1] = kernels::legendre::nodes[order][iy];
			for(int ix = 0; ix < basisSize; ix++) {
				q[0] = kernels::legendre::nodes[order][ix];
				
				r  = center + dx * ( q - 0.5 );
				Q  = &luh[idx_luh(iz,iy,ix,0)];
				// gradients in x, y, z direction
				Qx = &gradQ[idx_gradQ(iz,iy,ix,0,0)];
				Qy = &gradQ[idx_gradQ(iz,iy,ix,1,0)];
				Qz = &gradQ[idx_gradQ(iz,iy,ix,2,0)];
				ComputeAuxillaries(r.data(), Q, Qx, Qy, Qz);
			} // x
		} // y
	} // z
}

// For the fortran interface, abusing actually only the Constructor from
// ImportedTwoPunctures in order to have a single place where paremters are specified...
ImportedTwoPunctures* ImportedTP;
extern "C" {
void twopunctures_interpolate_(const double* x, double* Q) {
	if(!ImportedTP) {
		ImportedTP = new ImportedTwoPunctures();
		ImportedTP->prepare();
	}
	ImportedTP->tp->Interpolate(x, Q);
}
} // extern C


#endif /* Guard TWOPUNCTURES_AVAILABLE */

