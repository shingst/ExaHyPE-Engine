/**
 * This is the plotter abused to compute global spacetime integrals.
 * This is actually used to compare with the exact solution for
 * convergence tests.
 **/

#include "PDE/PDE.h"
#include "Writers/IntegralsWriter.h"
#include "InitialData/InitialData.h"
using SVEC::GRMHD::Cons2Prim;

#include "kernels/GaussLegendreBasis.h"
#include <cmath>

#include "kernels/aderdg/generic/c/sizes.cpph"

GRMHD::IntegralsWriter::IntegralsWriter(exahype::solvers::LimitingADERDGSolver&  solver)
	: IntegralsWriter() { plotForADERSolver = true; }

GRMHD::IntegralsWriter::IntegralsWriter(GRMHD::GRMHDSolver_ADERDG&  solver)
	: IntegralsWriter() { plotForADERSolver = true; }

GRMHD::IntegralsWriter::IntegralsWriter(GRMHD::GRMHDSolver_FV&  solver)
	: IntegralsWriter() { plotForADERSolver = false; }

GRMHD::IntegralsWriter::IntegralsWriter() :
	conserved("output/cons-"),
	primitives("output/prim-"),
	adm("output/const-adm-"), // for proof
	errors("output/error-"),
	statistics("output/volform")
{
	conserved.add(0, "dens");
	conserved.add(1, "sconx");
	conserved.add(2, "scony");
	conserved.add(3, "sconz");
	conserved.add(4, "ener");
	conserved.add(5, "bx");
	conserved.add(6, "by");
	conserved.add(7, "bz");
	conserved.add(8, "psi");
		
	primitives.add(0, "rho"); // V[0]=Q[0]
	primitives.add(1, "velx");
	primitives.add(2, "vely");
	primitives.add(3, "velz");
	primitives.add(4, "press");
	primitives.add(5, "bx");
	primitives.add(6, "by");
	primitives.add(7, "bz");
	primitives.add(8, "psi");
	
	// should all be static:
	adm.add(9, "lapse");
	adm.add(10, "shiftx");
	adm.add(11, "shifty");
	adm.add(12, "shiftz");
	adm.add(13, "gxx"); 	// gij ordering: Tensish (C)
	adm.add(14, "gxy");
	adm.add(15, "gyy");
	adm.add(16, "gxz");
	adm.add(17, "gyz");
	adm.add(18, "gzz");

	errors.add(0, "rho");
	errors.add(1, "velx");
	errors.add(2, "vely");
	errors.add(3, "velz");
	errors.add(4, "press");
	errors.add(5, "bx");
	errors.add(6, "by");
	errors.add(7, "bz");
	errors.add(8, "psi");
	// for proof
	adm.add(9, "lapse");
	adm.add(10, "shiftx");
	adm.add(11, "shifty");
	adm.add(12, "shiftz");
	adm.add(13, "gxx"); 	// gij ordering: Tensish (C)
	adm.add(14, "gxy");
	adm.add(15, "gyy");
	adm.add(16, "gxz");
	adm.add(17, "gyz");
	adm.add(18, "gzz");
}


GRMHD::IntegralsWriter::~IntegralsWriter() {
	// delete all the TimeSeriesReductions. Not that important
	// as this is a kind of a Singleton object.
}


void GRMHD::IntegralsWriter::startPlotting(double time) {
	conserved.startRow(time);
	primitives.startRow(time);
	errors.startRow(time);
	statistics.startRow(time);
}


void GRMHD::IntegralsWriter::finishPlotting() {
	conserved.finishRow();
	primitives.finishRow();
	errors.finishRow();
	statistics.finishRow();
}


void GRMHD::IntegralsWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
	// make sure this plotter has no output associated
	assertion( outputQuantities == nullptr );


	double dV;
	if(plotForADERSolver) {
		const int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
		dV = kernels::ADERDGVolume(order, sizeOfPatch, pos);
	} else {
		const int patchSize = GRMHD::AbstractGRMHDSolver_FV::PatchSize;
		dV = tarch::la::volume(sizeOfPatch)/patchSize; // correct is probably (patchSize+1)
	}

	statistics.addValue(dV, 1);


	// reduce the conserved quantities
	conserved.addValue(Q, dV);
	adm.addValue(Q, dV);

	// reduce the primitive quantities
	double V[nVar];
	constexpr bool force_crash = false;
	Cons2Prim(V, Q, force_crash).copyFullStateVector();
	primitives.addValue(V, dV);

	// now do the convergence test, as we have exact initial data
	
	// Quick ugly fix (SVEN): do not do the convergence test at the axis,
	//   for the 2D RNSID.
	if(pos[0]==0 || pos[1]==0)
		return;
	
	double ExactCons[nVar], ExactPrim[nVar];
	
	//id->Interpolate(xpos, timeStamp, ExactCons);
	InitialData(x.data(),timeStamp,ExactCons);
	Cons2Prim(ExactPrim, ExactCons).copyFullStateVector();
	
	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
	}
	
	errors.addValue(localError, dV);
}


