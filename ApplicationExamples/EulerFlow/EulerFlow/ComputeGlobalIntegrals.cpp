#include "ComputeGlobalIntegrals.h"
#include "Primitives.h"
#include "InitialData.h"

Euler::ComputeGlobalIntegrals::ComputeGlobalIntegrals(MyEulerSolver&  solver) :
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
	primitives.add(4, "press");
	
	errors.add(0, "rho");
	errors.add(1, "velx");
	errors.add(2, "vely");
	errors.add(3, "velz");
	errors.add(4, "press");

	// here, we plot all variables.
	assert( conserved.size() <= NumberOfVariables );
	assert( primitives.size() <= NumberOfVariables );
	assert( errors.size() <= NumberOfVariables );
}


Euler::ComputeGlobalIntegrals::~ComputeGlobalIntegrals() {
}


void Euler::ComputeGlobalIntegrals::startPlotting(double time) {
	conserved.startRow(time);
	primitives.startRow(time);
	errors.startRow(time);
	statistics.startRow(time);
}


void Euler::ComputeGlobalIntegrals::finishPlotting() {
	conserved.finishRow();
	primitives.finishRow();
	errors.finishRow();
	statistics.finishRow();
}


void Euler::ComputeGlobalIntegrals::mapQuantities(
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

	const double NumberOfLagrangePointsPerAxis = Euler::MyEulerSolver::Order + 1;
	//const double NumberOfUnknownsPerGridPoint = NumberOfVariables;

	// volume form for integration
	double scaling = tarch::la::volume(sizeOfPatch* (1.0/NumberOfLagrangePointsPerAxis));
	statistics.addValue(scaling, 1);
	
	// NOTE: This is *not* the correct scaling for the ADERDG scheme.
	// Instead, you should use the GlobalIntegralsLegendre plotter.

	// reduce the conserved quantities
	conserved.addValue(Q, scaling);

	// reduce the primitive quantities
	double V[NumberOfVariables];
	cons2prim(V, Q);
	primitives.addValue(V, scaling);

	// now do the convergence test, as we have exact initial data
	double ExactCons[NumberOfVariables];
	double ExactPrim[NumberOfVariables];
	const double *xpos = x.data();
	
	idfunc(xpos, ExactCons, timeStamp);
	cons2prim(ExactPrim, ExactCons);

	double localError[NumberOfVariables];
	for(int i=0; i<NumberOfVariables; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
	}
	errors.addValue(localError, scaling);
}


