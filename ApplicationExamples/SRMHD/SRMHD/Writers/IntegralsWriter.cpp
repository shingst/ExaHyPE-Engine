/**
 * This is the plotter abused to compute global spacetime integrals.
 * This is actually used to compare with the exact solution for
 * convergence tests.
 **/

#include "Writers/IntegralsWriter.h"
#include "Writers/TimeSeriesReductions.h"
#include "C2P-MHD.h"
#include "InitialDataAdapter.h"

#include "kernels/GaussLegendreBasis.h"
#include <cmath>

SRMHD::IntegralsWriter::IntegralsWriter(MHDSolver&  solver) :
	conserved("output/conserved-"),
	primitives("output/primitive-"),
	errors("output/error-"),
	statistics("output/volume.asc")
{
	// open all the reductions

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
	
	errors.add(0, "rho");
	errors.add(1, "velx");
	errors.add(2, "vely");
	errors.add(3, "velz");
	errors.add(4, "press");
	errors.add(5, "bx");
	errors.add(6, "by");
	errors.add(7, "bz");
	errors.add(8, "psi");

	// here, we plot all variables.
	assert( conserved.size() <= nVar );
	assert( primitives.size() <= nVar );
	assert( errors.size() <= nVar );
}


SRMHD::IntegralsWriter::~IntegralsWriter() {
}


void SRMHD::IntegralsWriter::startPlotting(double time) {
	conserved.startRow(time);
	primitives.startRow(time);
	errors.startRow(time);
	statistics.startRow(time);
}


void SRMHD::IntegralsWriter::finishPlotting() {
	conserved.finishRow();
	primitives.finishRow();
	errors.finishRow();
	statistics.finishRow();
}


void SRMHD::IntegralsWriter::mapQuantities(
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

	const double NumberOfLagrangePointsPerAxis = SRMHD::AbstractMHDSolver::Order + 1;
	//const double NumberOfUnknownsPerGridPoint = nVar;
	
	// Gauss-Legendre weights from pos argument
	double wx = kernels::legendre::weights[SRMHD::AbstractMHDSolver::Order][pos[0]];
	double wy = kernels::legendre::weights[SRMHD::AbstractMHDSolver::Order][pos[1]];
	double wz = 1;
	#ifdef Dim3
	wz = kernels::legendre::weights[SRMHD::AbstractMHDSolver::Order][pos[2]];
	#endif

	// volume form for integration
	double scaling = wx*wy*wz*tarch::la::volume(sizeOfPatch);
	statistics.addValue(scaling, 1);

	// reduce the conserved quantities
	conserved.addValue(Q, scaling);

	// reduce the primitive quantities
	double V[nVar];
	int err;
	pdecons2prim_(V, Q, &err);
	primitives.addValue(V, scaling);

	// now do the convergence test, as we have exact initial data
	double ExactCons[nVar];
	double ExactPrim[nVar];
	const double *xpos = x.data();
	
	//InitialData(xpos, Exact, time);
	// we compare against AlfenWave:
	alfenwave_(xpos, ExactCons, &timeStamp);
	pdecons2prim_(ExactPrim, ExactCons, &err);
	
	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
		//localError[i] = abs(Q[i] - ExactCons[i]);
	}

	errors.addValue(localError, scaling);
		
		/*if(i==0) {
			printf("RHO; t=%e x=(%e,%e) num=%14e exact=%14e error=%e\n", timeStamp, x[0],x[1], V[0], ExactPrim[0], localError[0]);
		}*/
    // Uncomment for debugging reasons
//    std::cout << "V_ana["<<i<<"]="<<ExactPrim[i]<<",V_h["<<i<<"]="<<V[i]<< std::endl;
//    std::cout << "Q_ana["<<i<<"]="<<ExactCons[i]<<",Q_h["<<i<<"]="<<Q[i] <<std::endl;
}


