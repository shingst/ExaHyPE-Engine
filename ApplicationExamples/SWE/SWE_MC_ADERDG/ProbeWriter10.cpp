// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter10.h"

SWE::ProbeWriter10::ProbeWriter10(SWE::MySWESolver& solver) {
  // @TODO Please insert your code here.
}

SWE::ProbeWriter10::~ProbeWriter10() {
}

void SWE::ProbeWriter10::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void SWE::ProbeWriter10::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ProbeWriter10::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
        outputQuantities[0] = Q[3];
    outputQuantities[1] = Q[3] + Q[4];
}
