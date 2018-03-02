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


GRGPR::TecplotWriter::TecplotWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
  // @TODO Please insert your code here.
}

GRGPR::TecplotWriter::TecplotWriter(GRGPR::GPRSolver_ADERDG& solver) : TecplotWriter(){
	plotForADERSolver = 1;
}

GRGPR::TecplotWriter::TecplotWriter(GRGPR::GPRSolver_FV& solver) : TecplotWriter(){
	plotForADERSolver = 0;
}

void GRGPR::TecplotWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  // @TODO Please insert your code here.
  elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&plotForADERSolver);

}


void GRGPR::TecplotWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
      int mpirank = tarch::parallel::Node::getInstance().getRank();
    //printf("rank %d\n", mpirank);
    initializetecplotplotter_(&time);
}


void GRGPR::TecplotWriter::finishPlotting() {
  // @TODO Please insert your code here.
  int mpirank = tarch::parallel::Node::getInstance().getRank();
  finishtecplotplotter_(&mpirank);
}
