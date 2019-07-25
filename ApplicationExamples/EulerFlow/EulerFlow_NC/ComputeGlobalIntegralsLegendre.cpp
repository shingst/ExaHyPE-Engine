#include "ComputeGlobalIntegralsLegendre.h"
#include "Primitives.h"
#include "InitialData.h"

#include "kernels/GaussLegendreBasis.h"

Euler::ComputeGlobalIntegralsLegendre::ComputeGlobalIntegralsLegendre(MyEulerSolver&  solver) :
	conserved("output/conserved-"),
	primitives("output/primitive-"),
	errors("output/error-"),
	statistics("output/volume.asc")	
{
	conserved.add(0, "dens");
	conserved.add(1, "sconx");
	conserved.add(2, "scony");
	conserved.add(3, "sconz");
	conserved.add(4, "ener");
	
	primitives.add(0, "rho"); // V[0]=Q[0]
	primitives.add(1, "velx");
	primitives.add(2, "vely");
	primitives.add(3, "velz");
	primitives.add(4, "presss");
	
	errors.add(0, "rho");
	errors.add(1, "velx");
	errors.add(2, "vely");
	errors.add(3, "velz");
	errors.add(4, "presss");

	// here, we plot all variables.
	assert( conserved.size() <= NumberOfVariables );
	assert( primitives.size() <= NumberOfVariables );
	assert( errors.size() <= NumberOfVariables );
}


Euler::ComputeGlobalIntegralsLegendre::~ComputeGlobalIntegralsLegendre() {
}


void Euler::ComputeGlobalIntegralsLegendre::startPlotting(double time) {
	conserved.startRow(time);
	primitives.startRow(time);
	errors.startRow(time);
	statistics.startRow(time);
}


void Euler::ComputeGlobalIntegralsLegendre::finishPlotting() {
	conserved.finishRow();
	primitives.finishRow();
	errors.finishRow();
	statistics.finishRow();
}


void Euler::ComputeGlobalIntegralsLegendre::mapQuantities(
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

	// volume form for integration
	double scaling = tarch::la::volume(sizeOfPatch);
	statistics.addValue(scaling, 1);

	// Gauss-Legendre weights from pos argument
	double wx = kernels::legendre::weights[Euler::MyEulerSolver::Order][pos[0]];
	double wy = kernels::legendre::weights[Euler::MyEulerSolver::Order][pos[1]];
	double wz = 1;
	#ifdef Dim3
	wz = kernels::legendre::weights[Euler::MyEulerSolver::Order][pos[2]];
	#endif
	scaling *= wx*wy*wz;

	// reduce the conserved quantities
	conserved.addValue(Q, scaling);

	// reduce the primitive quantities
	double V[NumberOfVariables];
	cons2prim(V, Q);
	primitives.addValue(V, scaling);

	// now do the convergence test, as we have exact initial data
	double Exact[NumberOfVariables];
	const double *xpos = x.data();
	
	idfunc(xpos, Exact, timeStamp);
	double ExactPrim[NumberOfVariables];
	cons2prim(ExactPrim, Exact);

	// Uncomment for debugging reasons
//	std::cout << "x="<<x.toString()<<"J="<<scaling << std::endl;
//  std::cout << "wx="<<wx<<",wy="<<wy<<",wz="<<wz << std::endl;

	double localError[NumberOfVariables];
	for(int i=0; i<NumberOfVariables; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
	}
	errors.addValue(localError, scaling);

		// Uncomment for debugging reasons
//		std::cout << "V_ana["<<i<<"]="<<ExactPrim[i]<<",V_h["<<i<<"]="<<V[i]<< std::endl;
//		std::cout << "Q_ana["<<i<<"]="<<Exact[i]<<",Q_h["<<i<<"]="<<Q[i] <<std::endl;
}


