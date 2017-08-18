// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorWriter.h"

#include "EulerSolver_FV.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"

#include "peano/utils/Loop.h"

#include "tarch/la/VectorOperations.h"

#include <algorithm>

#include <iomanip>

Euler::ErrorWriter::ErrorWriter() : exahype::plotters::FiniteVolumes2UserDefined::FiniteVolumes2UserDefined(){
  // @TODO Please insert your code here.
}


void Euler::ErrorWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {
  constexpr int numberOfVariables = AbstractEulerSolver_FV::NumberOfVariables;
  constexpr int basisSize         = AbstractEulerSolver_FV::PatchSize;
  constexpr int ghostLayerWidth   = AbstractEulerSolver_FV::GhostLayerWidth;

  double x[DIMENSIONS];

  kernels::idx4 idx(basisSize+2*ghostLayerWidth,basisSize+2*ghostLayerWidth,basisSize+2*ghostLayerWidth,numberOfVariables);
  dfor(i,basisSize) {
     double w_dV = 1.0;
     for (int d=0; d<DIMENSIONS; d++) {
       const double cellSize = sizeOfPatch[d] / basisSize;
       x[d]  = offsetOfPatch[d] + cellSize * (i(d)+0.5);
       w_dV *= cellSize;
     }

     double uAna[numberOfVariables];
     EulerSolver_FV::referenceSolution(x,timeStamp,uAna);

     const double* uNum = u +
         idx ( (DIMENSIONS==3) ? i(2)+ghostLayerWidth : 0, i(1)+ghostLayerWidth, i(0)+ghostLayerWidth, 0);

     for (int v=0; v<numberOfVariables; v++) {
        const double uDiff = std::abs(uNum[v]-uAna[v]);
        errorL2[v]   += uDiff*uDiff * w_dV;
        errorL1[v]   += uDiff * w_dV;
        errorLInf[v]  = std::max( errorLInf[v], uDiff );

        normL1Ana[v]  += std::abs(uAna[v]) * w_dV;
        normL2Ana[v]  += uAna[v] * uAna[v] * w_dV;
        normLInfAna[v] = std::max( normLInfAna[v], std::abs(uAna[v]) );
     }
  }
}

void Euler::ErrorWriter::startPlotting( double time) {
  _timeStamp = time;

  std::fill_n(errorL1,  AbstractEulerSolver_FV::NumberOfVariables, 0.0);
  std::fill_n(errorL2,  AbstractEulerSolver_FV::NumberOfVariables, 0.0);
  std::fill_n(errorLInf,AbstractEulerSolver_FV::NumberOfVariables, 0.0);
  
  std::fill_n(normL1Ana,  AbstractEulerSolver_FV::NumberOfVariables, 0.0);
  std::fill_n(normL2Ana,  AbstractEulerSolver_FV::NumberOfVariables, 0.0);
  std::fill_n(normLInfAna,AbstractEulerSolver_FV::NumberOfVariables, 0.0);
}

void Euler::ErrorWriter::finishPlotting() {
  constexpr int numberOfVariables = AbstractEulerSolver_FV::NumberOfVariables;

  for (int v=0; v<numberOfVariables; v++) {
    errorL2[v]   = sqrt(errorL2[v]);
    normL2Ana[v] = sqrt(normL2Ana[v]);
  }

  std::cout << "**Errors for FV solver with patch size="<<AbstractEulerSolver_FV::PatchSize<<"**" << std::endl;
  std::cout << "t_eval : "<<_timeStamp << std::endl;
  std::cout << "variable     : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << v << ", ";
  }
  std::cout << std::endl;

  std::cout << "absErrorL1   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL1[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "absErrorL2   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL2[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "absErrorLInf : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorLInf[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "relErrorL1   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL1[v]/normL1Ana[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "relErrorL2   : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorL2[v]/normL2Ana[v] << ", ";
  }
  std::cout << std::endl;

  std::cout << "relErrorLInf : ";
  for (int v=0; v<numberOfVariables; v++) {
    std::cout << std::setprecision(2) << errorLInf[v]/normLInfAna[v] << ", ";
  }
  std::cout << std::endl;
}
