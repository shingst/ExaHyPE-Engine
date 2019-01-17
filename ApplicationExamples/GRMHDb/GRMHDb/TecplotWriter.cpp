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
#include "tarch/parallel/NodePool.h"
#include <stdio.h>



GRMHDb::TecplotWriter::TecplotWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(){
  // @TODO Please insert your code here.
}

GRMHDb::TecplotWriter::TecplotWriter(GRMHDb::GRMHDbSolver_ADERDG& solver) : TecplotWriter(){
	plotForADERSolver = 1;
}



void GRMHDb::TecplotWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  // @TODO Please insert your code here.
  //counterloc++;
  //if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) { 
		elementcalltecplotplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0],&plotForADERSolver);
//		}
    //elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&counterloc2);
}


void GRMHDb::TecplotWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
  mpirank = tarch::parallel::Node::getInstance().getRank();
  //counterloc = 0;
  //counterloc2 = 0;
  //if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) { 
		initializetecplotplotter_(&time);
//		}
}


void GRMHDb::TecplotWriter::finishPlotting() {
  // @TODO Please insert your code here.
  //printf("1: finishtecplotplotter: DONE! %d ", mpirank);
  //fflush(stdout);
  //printf("  counterloc=%d", counterloc);
  //printf("  counterloc2=%d", counterloc2);

  //if (!tarch::parallel::NodePool::getInstance().isIdleNode(mpirank)) { 
		finishtecplotplotter_(&mpirank);
//		}
  //printf("2: finishtecplotplotter: DONE! %d ", mpirank);
  //fflush(stdout);
}