// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "DifferenceWithIDWriter.h"

#include "InitialData/InitialData.h"
using SVEC::GRMHD::Cons2Prim;


GRMHD::DifferenceWithIDWriter::DifferenceWithIDWriter(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}



GRMHD::DifferenceWithIDWriter::DifferenceWithIDWriter(GRMHDSolver_ADERDG&  solver) {
  // @todo Please insert your code here
}


GRMHD::DifferenceWithIDWriter::DifferenceWithIDWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

GRMHD::DifferenceWithIDWriter::~DifferenceWithIDWriter() {
}

void GRMHD::DifferenceWithIDWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GRMHD::DifferenceWithIDWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void GRMHD::DifferenceWithIDWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 23;
  
  double evolvedPrims[writtenUnknowns]={0.}, exactPrims[writtenUnknowns]={0.};
  Cons2Prim(evolvedPrims, Q).copyFullStateVector();

  double exactCons[writtenUnknowns]={0.};
  InitialData(x.data(),timeStamp,exactCons);
  Cons2Prim(exactPrims, exactCons).copyFullStateVector();

  double *localError = outputQuantities;
  for(int i=0; i<writtenUnknowns; i++) {
    localError[i] = evolvedPrims[i] - exactPrims[i];
  }
}