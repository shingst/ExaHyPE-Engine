#include "PrimitiveWriter.h"
#include "Fortran/PDE.h"
#include "GRMHDSolver_ADERDG_Variables.h"


GRMHD::PrimitiveWriter::PrimitiveWriter(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}



GRMHD::PrimitiveWriter::PrimitiveWriter(GRMHDSolver_ADERDG&  solver) {
  // @todo Please insert your code here
}


GRMHD::PrimitiveWriter::PrimitiveWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

GRMHD::PrimitiveWriter::~PrimitiveWriter() {
  // @todo Please insert your code here
}

void GRMHD::PrimitiveWriter::writtenQuantitiesNames(char** outputfileNames) {
	// also works in principle, but not as convenient
	// GRMHD::AbstractGRMHDSolver_ADERDG::VariableNames specfileNames;
	// for(int i=0; i<outputQuantities; i++) outputfileNames[i] = specfileNames[i];
	outputfileNames[0] = "rho";
	outputfileNames[1] = "velx";
	outputfileNames[2] = "vely";
	outputfileNames[3] = "velz";
	outputfileNames[4] = "press";
	outputfileNames[5] = "bx";
	outputfileNames[6] = "by";
	outputfileNames[7] = "bz";
	outputfileNames[8] = "psi";
	outputfileNames[9] = "lapse";
	outputfileNames[10] = "shiftx";
	outputfileNames[11] = "shifty";
	outputfileNames[12] = "shiftz";
	// fixme: Check ordering of matrix. In GRMHD, we have Fortran ordering (row major)
	// not tensish or C++ ordering (col major). I'm currently too sick to do that.
	outputfileNames[13] = "gxx";
	outputfileNames[14] = "gxy";
	outputfileNames[15] = "gyy";
	outputfileNames[16] = "gxz";
	outputfileNames[17] = "gyz";
	outputfileNames[18] = "gzz";
	// Info: This notation spills errors, cf https://stackoverflow.com/q/20944784/1656042
	// for solution. Again, I'm too lazy. It works.
}

void GRMHD::PrimitiveWriter::startPlotting(double time) {
  numberOfC2PFailures = 0;
  allConversions = 0;
}


void GRMHD::PrimitiveWriter::finishPlotting() {
  if(numberOfC2PFailures > 0) {
    static tarch::logging::Log _log("GRMHD::PrimitiveWriter");
    logInfo("finishPlotting", "Counted " << numberOfC2PFailures << " Cons2Prim failures during plotting (" << (allConversions/numberOfC2PFailures*100) << "%)");
  }
  startPlotting(0);
}


void GRMHD::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
    int err;
    pdecons2prim_(outputQuantities, Q, &err);
    if(err != 0) {
	    //printf("Cons2Prim Failure in PrimitiveWriter!!!");
	    numberOfC2PFailures++;

	    constexpr int c2pFailureIndicatingMagicNumberForDensity = -1;
	    Q[0] = c2pFailureIndicatingMagicNumberForDensity;
    }
    allConversions++;
}


