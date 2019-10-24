#include "ConstraintsWriter.h"

#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h"

#include "Fortran/PDE.h" // ADMConstraints()

#include "AbstractCCZ4Solver_ADERDG.h"


CCZ4::ConstraintsWriter::ConstraintsWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

CCZ4::ConstraintsWriter::ConstraintsWriter(CCZ4Solver_FV&  solver) {
  // @todo Please insert your code here
}

CCZ4::ConstraintsWriter::ConstraintsWriter(CCZ4Solver_ADERDG&  solver) {
  // @todo Please insert your code here
}

CCZ4::ConstraintsWriter::~ConstraintsWriter() {
  // @todo Please insert your code here
}


void CCZ4::ConstraintsWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void CCZ4::ConstraintsWriter::finishPlotting() {
  // @todo Please insert your code here
}

bool CCZ4::ConstraintsWriter::mapWithDerivatives() {
  return true;
}

void CCZ4::ConstraintsWriter::mapQuantities(
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
	
	CCZ4Fortran::DeriveADMConstraints(outputQuantities, Q, gradQ);
}

