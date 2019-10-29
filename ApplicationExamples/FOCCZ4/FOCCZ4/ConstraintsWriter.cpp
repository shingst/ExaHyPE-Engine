#include "ConstraintsWriter.h"

#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h"

#include "PDE.h" // ADMConstraints()

#include "AbstractFOCCZ4Solver_ADERDG.h"
#include "kernels/aderdg/generic/c/computeGradients.cpph"
#include "kernels/GaussLegendreBasis.h"

#include "tarch/logging/Log.h"
static tarch::logging::Log _log("FOCCZ4::ConstraintsWriter::");


FOCCZ4::ConstraintsWriter::ConstraintsWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

FOCCZ4::ConstraintsWriter::ConstraintsWriter(FOCCZ4Solver_FV&  solver) {
  // @todo Please insert your code here
}

FOCCZ4::ConstraintsWriter::ConstraintsWriter(FOCCZ4Solver_ADERDG&  solver) {
  // @todo Please insert your code here
}

FOCCZ4::ConstraintsWriter::~ConstraintsWriter() {
  // @todo Please insert your code here
}


void FOCCZ4::ConstraintsWriter::startPlotting(double time) {
  // @todo Please insert your code here
  timeLastWarned = -1;
}


void FOCCZ4::ConstraintsWriter::finishPlotting() {
  // @todo Please insert your code here
}

bool FOCCZ4::ConstraintsWriter::mapWithDerivatives() {
  return true;
}

void FOCCZ4::ConstraintsWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q, double* gradQ,
    double* outputQuantities,
    double timeStamp
) {
	// this should be a UserOnTheFlyPostProcessing constant,
	// allowing to ensure we write out 6 unknowns.
	static constexpr int writtenUnknowns = 6;

	for(int i=0; i<59; i++) { if(Q[i]!=Q[i]) std::abort(); }
	
	admconstraints_(outputQuantities, Q, gradQ);

	for(int i=0; i<6; i++) { if(outputQuantities[i]!=outputQuantities[i]) std::abort(); }
}


void FOCCZ4::ConstraintsWriter::mapQuantities(
        const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, int>&    pos,
        double* Q,
        double* outputQuantities,
        double timeStamp) {
	// this should be a UserOnTheFlyPostProcessing constant,
	// allowing to ensure we write out 6 unknowns.
	static constexpr int writtenUnknowns = 6;

	//for(int i=0; i<59; i++) { if(Q[i]!=Q[i]) std::abort(); }

    double gradQ[basisSize3 * DIMENSIONS * numberOfVariables];
	
	//kernels::aderdg::generic::c::computeGradQ<FOCCZ4::AbstractFOCCZ4Solver_ADERDG>(gradQ, Q, sizeOfPatch);

	admconstraints_(outputQuantities, Q, gradQ);

	for(int i=0; i<6; i++) { if(outputQuantities[i]!=outputQuantities[i]) std::abort(); }

}
	
