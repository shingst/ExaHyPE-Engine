// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "MHDSolver_FV.h"
#include "kernels/finitevolumes/godunov/c/2d/godunov.cpph"
#include "kernels/finitevolumes/godunov/c/3d/godunov.cpph"

MHD::MHDSolver_FV::MHDSolver_FV(int cellsPerCoordinateAxis,double maximumMeshSize,exahype::solvers::Solver::TimeStepping timeStepping, exahype::Parser::ParserView constants):
  exahype::solvers::FiniteVolumesSolver("MHDSolver_FV",9,0,cellsPerCoordinateAxis,1 /* ghost layer width */,maximumMeshSize,timeStepping) {
  // @todo Please implement/augment if required
}

double MHD::MHDSolver_FV::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double maxAdmissibleDt = kernels::finitevolumes::godunov::c::stableTimeStepSize<eigenvalues>(luh,tempEigenvalues,dx,getNumberOfVariables(),getNodesPerCoordinateAxis());
  return maxAdmissibleDt;
}

void MHD::MHDSolver_FV::solutionUpdate(double* luhNew,const double* luh,double** tempStateSizedArrays,double** tempUnknowns,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt,double& maxAdmissibleDt) {
  maxAdmissibleDt = kernels::finitevolumes::godunov::c::solutionUpdate<flux,source,eigenvalues>(luhNew,luh,tempStateSizedArrays,tempUnknowns,dx,dt,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

void MHD::MHDSolver_FV::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  kernels::finitevolumes::godunov::c::solutionAdjustment<adjustedSolutionValues>(luh,center,dx,t,dt,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

void MHD::MHDSolver_FV::ghostLayerFilling(double* luh,const double* luhNeighbour,const tarch::la::Vector<DIMENSIONS,int>& neighbourPosition) {
  kernels::finitevolumes::godunov::c::ghostLayerFilling(luh,luhNeighbour,neighbourPosition,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

void MHD::MHDSolver_FV::ghostLayerFillingAtBoundary(double* luh,const double* luhbnd,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) {
  kernels::finitevolumes::godunov::c::ghostLayerFillingAtBoundary(luh,luhbnd,boundaryPosition,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

void MHD::MHDSolver_FV::boundaryLayerExtraction(double* luhbnd,const double* luh,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) {
  kernels::finitevolumes::godunov::c::boundaryLayerExtraction(luhbnd,luh,boundaryPosition,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

void MHD::MHDSolver_FV::boundaryConditions(double* luhbndOutside,const double* const luhbndInside,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {
  kernels::finitevolumes::godunov::c::boundaryConditions(*this,luhbndOutside,luhbndInside,cellCentre,cellSize,t,dt,faceIndex,normalNonZero);
}