// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "TecplotWriter.h"
#include "TECPLOTinterface.h"
#include "tarch/parallel/Node.h"
#include <stdio.h>

DIM::TecplotWriter::TecplotWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
  // @TODO Please insert your code here.
}


void DIM::TecplotWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  // @TODO Please insert your code here.
  //plotForADERSolver=1;
  elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&plotForADERSolver);
}


void DIM::TecplotWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
  initializetecplotplotter_(&time);
}


void DIM::TecplotWriter::finishPlotting() {
  // @TODO Please insert your code here.
  int mpirank = tarch::parallel::Node::getInstance().getRank();
  finishtecplotplotter_(&mpirank);
}