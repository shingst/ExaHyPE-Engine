#include "PrimitiveWriter.h"
#include "AbstractGRMHDSolver_ADERDG.h"
#include "PDE/PDE-GRMHD-ExaHyPE.h"
using SVEC::GRMHD::Cons2Prim;

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


void GRMHD::PrimitiveWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void GRMHD::PrimitiveWriter::finishPlotting() {
  // @todo Please insert your code here
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
	for (int i=0; i<23; i++) outputQuantities[i] = -42;
	constexpr bool force_crash = false;
	Cons2Prim c2p(outputQuantities, Q, force_crash);
	c2p.copyFullStateVector();
	
	// abuse the last output quantity to inform about the C2P status
	constexpr int c2pstatus_field = 22;
	outputQuantities[c2pstatus_field] = c2p.failed ? 1 : 0;
}

