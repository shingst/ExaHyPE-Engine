#include "MyEulerSolver_Plotter1.h"

/**
 * This is the plotter abused to compute global spacetime integrals.
 * This is actually used to compare with the exact solution for
 * convergence tests.
 * 
 * 
 **/

#include "TimeSeriesReductions.h"
#include "Primitives.h"
#include "InitialData.h"

Euler::MyEulerSolver_Plotter1::MyEulerSolver_Plotter1(MyEulerSolver& solver) {
	// open all the reductions
	assert( 5 == nVars );
	
	conserved[0] = new TimeSeriesReductions("output/dens.asc");
	conserved[1] = new TimeSeriesReductions("output/sconx.asc");
	conserved[2] = new TimeSeriesReductions("output/scony.asc");
	conserved[3] = new TimeSeriesReductions("output/sconz.asc");
	conserved[4] = new TimeSeriesReductions("output/ener.asc");
	
	primitives[0] = new TimeSeriesReductions("output/rho.asc"); // V[0]=Q[0]
	primitives[1] = new TimeSeriesReductions("output/velx.asc");
	primitives[2] = new TimeSeriesReductions("output/vely.asc");
	primitives[3] = new TimeSeriesReductions("output/velz.asc");
	primitives[4] = new TimeSeriesReductions("output/press.asc");
	
	errors[0] = new TimeSeriesReductions("output/error-rho.asc");
	errors[1] = new TimeSeriesReductions("output/error-velx.asc");
	errors[2] = new TimeSeriesReductions("output/error-vely.asc");
	errors[3] = new TimeSeriesReductions("output/error-velz.asc");
	errors[4] = new TimeSeriesReductions("output/error-press.asc");
	
	statistics = new TimeSeriesReductions("output/volform.asc");
}


Euler::MyEulerSolver_Plotter1::~MyEulerSolver_Plotter1() {
	// delete all reductions.
}


void Euler::MyEulerSolver_Plotter1::startPlotting(double time) {
	this->time = time;
	for(int i=0; i<nVars; i++) {
		conserved[i]->initRow(time);
		primitives[i]->initRow(time);
		errors[i]->initRow(time);
	}
	statistics->initRow(time);
}


void Euler::MyEulerSolver_Plotter1::finishPlotting() {
	for(int i=0; i<nVars; i++) {
		conserved[i]->writeRow();
		primitives[i]->writeRow();
		errors[i]->writeRow();
	}
	statistics->writeRow();
}

void Euler::MyEulerSolver_Plotter1::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
	// make sure this plotter has no output associated
	assertion( outputQuantities == nullptr );

	const double NumberOfLagrangePointsPerAxis = MY_POLYNOMIAL_DEGREE + 1;
	const double NumberOfUnknownsPerGridPoint = MY_NUMBER_OF_VARIABLES;

	// volume form for integration
	double scaling = tarch::la::volume(sizeOfPatch* (1.0/NumberOfLagrangePointsPerAxis));
	statistics->addValue(scaling, 1);

	// reduce the conserved quantities
	for (int i=0; i<nVars; i++)
		conserved[i]->addValue( Q[i], scaling );

	// reduce the primitive quantities
	double V[nVars];
	cons2prim(V, Q);
	for(int i=0; i<nVars; i++)
		primitives[i]->addValue( V[i], scaling );

	// now do the convergence test, as we have exact initial data
	double Exact[nVars];
	// TODO: Need a way to access _data in tarch::la::Vector.
	double xpos[DIMENSIONS];
	for(int i=0; i<DIMENSIONS; i++) xpos[i] = x[i];
	
//	ShuVortex2D(xpos, Exact, time); // TODO(Dominic): Bug: Here, you used the Shu-Vortex as reference solution for every problem. You compared your (Euler) lin. adv. against it.
	InitialData(xpos, Exact, time);

	double localError[nVars];
	for(int i=0; i<nVars; i++) {
		localError[i] = abs(V[i] - Exact[i]);
		errors[i]->addValue( localError[i], scaling );
	}
	
}


