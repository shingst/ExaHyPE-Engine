#include "PrimitiveWriter.h"
#include "CCZ4Solver_ADERDG.h"
#include "CCZ4Solver_FV.h"

#include "Fortran/PDE.h"

CCZ4::PrimitiveWriter::PrimitiveWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

CCZ4::PrimitiveWriter::PrimitiveWriter(CCZ4Solver_FV&  solver) {
  // @todo Please insert your code here
}

CCZ4::PrimitiveWriter::PrimitiveWriter(CCZ4Solver_ADERDG&  solver) {
  // @todo Please insert your code here
}

CCZ4::PrimitiveWriter::~PrimitiveWriter() {
  // @todo Please insert your code here
}


void CCZ4::PrimitiveWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void CCZ4::PrimitiveWriter::finishPlotting() {
  // @todo Please insert your code here
}


void CCZ4::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {

  // log(alpha) and all that 
  CCZ4Fortran::Cons2Prim(outputQuantities, Q);
}
