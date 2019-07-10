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
#include "GPRDRSolver_FV.h"
#include "GPRDRSolver_ADERDG.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorOperations.h"

#include <algorithm>

#include <iomanip>
GPRDR::TecplotWriter::TecplotWriter() : exahype::plotters::LimitingADERDG2UserDefined::LimitingADERDG2UserDefined(){
  // @TODO Please insert your code here.
}

void GPRDR::TecplotWriter::plotADERDGPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
    double timeStamp) {
  // @TODO Please insert your code here.
  plotForADERSolver = 0;
  elementcalltecplotplotter_(u,&offsetOfPatch[0],&sizeOfPatch[0],&plotForADERSolver);
}

void GPRDR::TecplotWriter::plotFiniteVolumesPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* const u,
    double timeStamp) {
  // @TODO Please insert your code here.
  plotForADERSolver = 1;
  elementcalltecplotfvplotter_(u, &offsetOfPatch[0], &sizeOfPatch[0], &plotForADERSolver);
}


void GPRDR::TecplotWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
  mpirank = tarch::parallel::Node::getInstance().getRank();
  initializetecplotplotter_(&time);
}


void GPRDR::TecplotWriter::finishPlotting() {
  // @TODO Please insert your code here.
  finishtecplotplotter_(&mpirank);
}