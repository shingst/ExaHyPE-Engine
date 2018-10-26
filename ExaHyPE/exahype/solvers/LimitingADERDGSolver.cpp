/**onte
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 *
 * \author Dominic E. Charrier, Tobias Weinzierl
 **/

#include <algorithm> //copy_n

#include "LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"
#include "exahype/amr/AdaptiveMeshRefinement.h"

namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    exahype::solvers::ADERDGSolver* solver,
    exahype::solvers::FiniteVolumesSolver* limiter,
    const double DMPRelaxationParameter,
    const double DMPDifferenceScaling,
    const int iterationsToCureTroubledCell)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNodesPerCoordinateAxis(), solver->getMaximumMeshSize(),
          solver->getMaximumAdaptiveMeshDepth(),
          solver->getTimeStepping()),
          _solver(std::move(solver)),
          _limiter(std::move(limiter)),
          _DMPMaximumRelaxationParameter(DMPRelaxationParameter),
          _DMPDifferenceScaling(DMPDifferenceScaling),
          _iterationsToCureTroubledCell(iterationsToCureTroubledCell)
{
  solver->disableCheckForNaNs();

  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());

  #ifdef Parallel
  // TODO(WORKAROUND): Not sure for what anymore
  const int numberOfObservables = _solver->getDMPObservables();
  _invalidObservables.resize(numberOfObservables);
  std::fill_n(_invalidObservables.data(),_invalidObservables.size(),-1);

  _receivedMax.resize(numberOfObservables);
  _receivedMin.resize(numberOfObservables);
  assertion( numberOfObservables==0 || !_receivedMax.empty());
  assertion( numberOfObservables==0 || !_receivedMin.empty());
  #endif
}
/** Wire through to limiting ADER-DG solver */
exahype::solvers::Solver::MeshUpdateEvent exahype::solvers::LimitingADERDGSolver::getNextMeshUpdateEvent() const {
  return _solver->getNextMeshUpdateEvent();
}
void exahype::solvers::LimitingADERDGSolver::updateNextMeshUpdateEvent(exahype::solvers::Solver::MeshUpdateEvent meshUpdateEvent) {
  return _solver->updateNextMeshUpdateEvent(meshUpdateEvent);
}
void exahype::solvers::LimitingADERDGSolver::setNextMeshUpdateEvent() {
  _solver->setNextMeshUpdateEvent();
}
exahype::solvers::Solver::MeshUpdateEvent exahype::solvers::LimitingADERDGSolver::getMeshUpdateEvent() const {
  return _solver->getMeshUpdateEvent();
}
void exahype::solvers::LimitingADERDGSolver::overwriteMeshUpdateEvent(MeshUpdateEvent newMeshUpdateEvent) {
   _solver->overwriteMeshUpdateEvent(newMeshUpdateEvent);
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStamp() const {
  return _solver->getMinTimeStamp();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStepSize() const {
  return _solver->getMinTimeStepSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinNextTimeStepSize() const {
  return _solver->getMinNextTimeStepSize();
}

void exahype::solvers::LimitingADERDGSolver::updateMinNextTimeStepSize(double value) {
  _solver->updateMinNextTimeStepSize(value);
}

void exahype::solvers::LimitingADERDGSolver::initSolver(
    const double timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
    const std::vector<std::string>& cmdlineargs,
    const exahype::parser::ParserView& parserView) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  std::pair<double,int> coarsestMeshInfo =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(_maximumMeshSize,boundingBoxSize[0]);
  _coarsestMeshSize  = coarsestMeshInfo.first; // TODO(Dominic): Wire through as well
  _coarsestMeshLevel = coarsestMeshInfo.second;

  updateNextMeshUpdateEvent(MeshUpdateEvent::InitialRefinementRequested);
  setNextMeshUpdateEvent();

  _solver->initSolver(timeStamp, domainOffset, domainSize, boundingBoxSize, cmdlineargs, parserView);
  _limiter->initSolver(timeStamp, domainOffset, domainSize, boundingBoxSize, cmdlineargs, parserView);
}

bool exahype::solvers::LimitingADERDGSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  return _solver->isPerformingPrediction(section);
}

bool exahype::solvers::LimitingADERDGSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  bool isMergingMetadata = false;

  switch (section) {
    case exahype::State::AlgorithmSection::RefinementStatusSpreading:
      isMergingMetadata = getMeshUpdateEvent()==MeshUpdateEvent::IrregularLimiterDomainChange ||
                          getMeshUpdateEvent()==MeshUpdateEvent::RefinementRequested;
      break;
    case exahype::State::AlgorithmSection::MeshRefinement:
      isMergingMetadata = hasRequestedMeshRefinement();
      break;
    default:
      break;
  }

  return isMergingMetadata;
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  _solver->synchroniseTimeStepping(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex) const {
  _solver->synchroniseTimeStepping(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellDescriptionsIndex);
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    SolverPatch& solverPatch,
    FiniteVolumesSolver::Heap::HeapEntries& limiterPatches,
    const int limiterElement) const {
  _solver->synchroniseTimeStepping(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,limiterPatches,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  _solver->startNewTimeStep();
  ensureLimiterTimeStepDataIsConsistent();

  logDebug("startNewTimeStep()","getMeshUpdateEvent()="<<Solver::toString(getMeshUpdateEvent())<<
           ",getNextMeshUpdateEvent()="<<Solver::toString(getNextMeshUpdateEvent()));
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  _solver->startNewTimeStepFused(isFirstIterationOfBatch,isLastIterationOfBatch);
  ensureLimiterTimeStepDataIsConsistent();

  logDebug("startNewTimeStep()","getMeshUpdateEvent()="<<Solver::toString(getMeshUpdateEvent())<<
           ",getNextMeshUpdateEvent()="<<Solver::toString(getNextMeshUpdateEvent()));
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterTimeStepDataIsConsistent() const {
  _limiter->_minTimeStamp            = _solver->_minCorrectorTimeStamp;
  _limiter->_minTimeStepSize         = _solver->_minCorrectorTimeStepSize;
  _limiter->_previousMinTimeStamp    = _solver->_previousMinCorrectorTimeStamp;
  _limiter->_previousMinTimeStepSize = _solver->_previousMinCorrectorTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::updateTimeStepSizesFused()  {
  _solver->updateTimeStepSizesFused();
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::updateTimeStepSizes()  {
  _solver->updateTimeStepSizes();
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes() {
  _solver->zeroTimeStepSizes();
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStepFused() {
  _solver->rollbackToPreviousTimeStepFused();
  ensureLimiterTimeStepDataIsConsistent();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMaxLevel(int maxLevel) {
  _solver->updateNextMaxLevel(maxLevel);
}

int exahype::solvers::LimitingADERDGSolver::getNextMaxLevel() const {
  return _solver->getNextMaxLevel();
}

int exahype::solvers::LimitingADERDGSolver::getMaxLevel() const {
  return _solver->getMaxLevel();
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////

exahype::solvers::Solver::MeshUpdateEvent exahype::solvers::LimitingADERDGSolver::updateRefinementStatusDuringRefinementStatusSpreading(
    SolverPatch& solverPatch) const {
  _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
  if ( 
      solverPatch.getType()==SolverPatch::Type::Descendant &&
      solverPatch.getRefinementStatus() > 0 &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()
  ) {
    solverPatch.setRefinementFlag(true);
    return MeshUpdateEvent::RefinementRequested;
  } else {
    return MeshUpdateEvent::None;
  }
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber,
    const bool stillInRefiningMode) {
 return
     _solver->progressMeshRefinementInEnterCell(
         fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
         coarseGridCell,coarseGridVerticesEnumerator,
         solverNumber,stillInRefiningMode);
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber,
    const bool stillInRefiningMode) {
  return
      _solver->progressMeshRefinementInLeaveCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,fineGridPositionOfCell,solverNumber,stillInRefiningMode);
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::LimitingADERDGSolver::eraseOrRefineAdjacentVertices(
        const int cellDescriptionsIndex,
        const int solverNumber,
        const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize,
        const bool checkThoroughly) const {
  return _solver->eraseOrRefineAdjacentVertices(
             cellDescriptionsIndex,solverNumber,cellOffset,cellSize,checkThoroughly);
}

bool exahype::solvers::LimitingADERDGSolver::attainedStableState(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int solverNumber) const {
  return
      _solver->attainedStableState(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) {
  _solver->finaliseStateUpdates(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,solverNumber);

  const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
  const int solverElement = _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  if ( solverElement!=exahype::solvers::Solver::NotFound ) {
    SolverPatch& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,solverElement);
    if (
        solverPatch.getType()==SolverPatch::Type::Cell &&
        solverPatch.getRefinementStatus()<-1
    ) {
      logError("determineMinAndMax(...)","solverPatch.getRefinementStatus()<-1 for cell="<<solverPatch.toString());
      std::abort();
    }

   // only done when doing initial refinement
   const bool newLimiterPatchAllocated =
       getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested &&
       ensureRequiredLimiterPatchIsAllocated(
           solverPatch,cellDescriptionsIndex,solverPatch.getRefinementStatus());
   if ( newLimiterPatchAllocated ) {
     assertion1(tarch::la::equals(solverPatch.getCorrectorTimeStamp(),0.0),solverPatch.toString());
     const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverNumber);
     LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
     adjustLimiterSolution(solverPatch,limiterPatch);
   }
    ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellDescriptionsIndex);
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellDescriptionsIndex);
  }
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex) {
  double admissibleTimeStepSize =
      _solver->startNewTimeStep(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellDescriptionsIndex);

  return admissibleTimeStepSize;
}

double exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch)  {
  double admissibleTimeStepSize =
      _solver->startNewTimeStepFused(solverPatch,isFirstIterationOfBatch,isLastIterationOfBatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellDescriptionsIndex);

  return admissibleTimeStepSize;
}

double exahype::solvers::LimitingADERDGSolver::updateTimeStepSizesFused(
      const int cellDescriptionsIndex,
      const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const double admissibleTimeStepSize =
        _solver->updateTimeStepSizesFused(cellDescriptionsIndex,solverElement);

    ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

double exahype::solvers::LimitingADERDGSolver::updateTimeStepSizes(
      const int cellDescriptionsIndex,
      const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const double admissibleTimeStepSize =
        _solver->updateTimeStepSizes(cellDescriptionsIndex,solverElement);

    ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex) const {
  _solver->zeroTimeStepSizes(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellDescriptionsIndex);
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  synchroniseTimeStepping(cellDescriptionsIndex, solverElement); // TODO(Dominic): have version with solver patch
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  _solver->rollbackToPreviousTimeStep(solverPatch); // TODO(Dominic): Too many element lookups
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStepFused(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  synchroniseTimeStepping(cellDescriptionsIndex, solverElement); // TODO(Dominic): have version with solver patch
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  _solver->rollbackToPreviousTimeStepFused(solverPatch); // TODO(Dominic): Too many element lookups
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::adjustSolutionDuringMeshRefinementBody(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex,
    const bool isInitialMeshRefinement) {
  zeroTimeStepSizes(solverPatch,cellDescriptionsIndex);      // TODO(Dominic): Still necessary?
  synchroniseTimeStepping(solverPatch,cellDescriptionsIndex);

  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    if (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::Prolongating) {
      _solver->prolongateVolumeData(solverPatch,isInitialMeshRefinement);
      solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
    }
    _solver->adjustSolution(solverPatch);

    determineSolverMinAndMax(solverPatch,false);
    if ( !evaluatePhysicalAdmissibilityCriterion(solverPatch) ) {
       solverPatch.setRefinementStatus(_solver->getMinimumRefinementStatusForTroubledCell());
       solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
    } else {
      _solver->markForRefinement(solverPatch);
    }
  }
}

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatchForSolverPatch(
    const SolverPatch& solverPatch,const int cellDescriptionsIndex) const {
  assertion1(solverPatch.getRefinementStatus()>=0,solverPatch.toString());
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  assertion1(limiterElement!=Solver::NotFound,solverPatch.toString());
  LimiterPatch& limiterPatch =
      FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  // Ensure time stamps and step sizes are consistent
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return limiterPatch;
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
  const SolverPatch& solverPatch,const int cellDescriptionsIndex) const {
  const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if ( limiterElement!=Solver::NotFound ) {
    #ifdef Asserts
    LimiterPatch& limiterPatch = FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
    #endif
    assertion2(solverPatch.getPreviousRefinementStatus() >=_solver->_minimumRefinementStatusForPassiveFVPatch ||
               solverPatch.getRefinementStatus() >=_solver->_minimumRefinementStatusForPassiveFVPatch,solverPatch.toString(),limiterPatch.toString());
    copyTimeStepDataFromSolverPatch(solverPatch,cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if ( limiterElement!=Solver::NotFound ) {
    #ifdef Asserts
    LimiterPatch& limiterPatch = FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
    #endif
    assertion2(solverPatch.getPreviousRefinementStatus() >=_solver->_minimumRefinementStatusForPassiveFVPatch ||
               solverPatch.getRefinementStatus() >=_solver->_minimumRefinementStatusForPassiveFVPatch,solverPatch.toString(),limiterPatch.toString());
    copyTimeStepDataFromSolverPatch(solverPatch,cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
    const SolverPatch& solverPatch,
    FiniteVolumesSolver::Heap::HeapEntries& limiterPatches,
    const int limiterElement) const {
  if ( limiterElement!=Solver::NotFound ) {
    LimiterPatch& limiterPatch = limiterPatches[limiterElement];
    assertion2(solverPatch.getPreviousRefinementStatus() >=_solver->_minimumRefinementStatusForPassiveFVPatch ||
               solverPatch.getRefinementStatus() >=_solver->_minimumRefinementStatusForPassiveFVPatch,solverPatch.toString(),limiterPatch.toString());
    copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  }
}

void exahype::solvers::LimitingADERDGSolver::copyTimeStepDataFromSolverPatch(
    const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement) {
  LimiterPatch& limiterPatch =
      FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
}

void exahype::solvers::LimitingADERDGSolver::copyTimeStepDataFromSolverPatch(
    const SolverPatch& solverPatch, LimiterPatch& limiterPatch) {
  limiterPatch.setPreviousTimeStamp(solverPatch.getPreviousCorrectorTimeStamp());
  limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
  limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
  limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStepOrRestriction(
    const int  solverNumber,
    CellInfo&  cellInfo,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isAtRemoteBoundary) {
  const int element        = indexOfCellDescription(cellInfo._ADERDGCellDescriptions,solverNumber);
  SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[element];

  // Write the previous limiter status back onto the patch for all cell description types
  solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    const bool isAMRSkeletonCell     = ADERDGSolver::belongsToAMRSkeleton(solverPatch,isAtRemoteBoundary);
    const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
    const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

    if (
        !SpawnPredictionAsBackgroundJob ||
        isFirstIterationOfBatch ||
        isLastIterationOfBatch  ||
        mustBeDoneImmediately
    ) {
      return fusedTimeStepBody(
          solverPatch, cellInfo._cellDescriptionsIndex, element,
          isFirstIterationOfBatch,isLastIterationOfBatch,
          isSkeletonCell,
          mustBeDoneImmediately,
          solverPatch.getNeighbourMergePerformed());
    } else {
      solverPatch.setHasCompletedTimeStep(false); // done here in order to skip lookup of cell description in job constructor
      FusedTimeStepJob fusedTimeStepJob( *this, solverPatch, cellInfo._cellDescriptionsIndex, element,
          solverPatch.getNeighbourMergePerformed(), isSkeletonCell );
      Solver::submitJob(fusedTimeStepJob,isSkeletonCell);
      return UpdateResult();
    }
  }
  else {
    UpdateResult result;

    if (
        solverPatch.getType()==SolverPatch::Type::Descendant &&
        solverPatch.getCommunicationStatus()>=ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication
    ) {
      _solver->restrictToTopMostParent(solverPatch);
    }

    _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
    result._meshUpdateEvent =
        _solver->evaluateRefinementCriteriaAfterSolutionUpdate(
            solverPatch,solverPatch.getNeighbourMergePerformed()); // must be done by all cell types TODO(Dominic): Clean up
    return result;
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStepBody(
    SolverPatch& solverPatch,
    const int    cellDescriptionsIndex,
    const int    element,
    const bool   isFirstIterationOfBatch,
    const bool   isLastIterationOfBatch,
    const bool   isSkeletonJob,
    const bool   mustBeDoneImmediately,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed) {
  // synchroniseTimeStepping(cellDescriptionsIndex,element); // assumes this was done in neighbour merge
  updateSolution(solverPatch,cellDescriptionsIndex,neighbourMergePerformed,isFirstIterationOfBatch);

  UpdateResult result;
  result._timeStepSize = startNewTimeStepFused(
      solverPatch,cellDescriptionsIndex,isFirstIterationOfBatch,isLastIterationOfBatch);
  result._meshUpdateEvent =
      updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
          solverPatch,cellDescriptionsIndex,neighbourMergePerformed);

  if ( solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForTroubledCell() ) {
    _solver->performPredictionAndVolumeIntegral(
        solverPatch,cellDescriptionsIndex,element,
        solverPatch.getCorrectorTimeStamp(), // corrector time step data is correct; see docu
        solverPatch.getCorrectorTimeStepSize(),
        false/*already uncompressed*/,isSkeletonJob);
  }

  return result;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::update(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary){
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  
  // Write the previous limiter status back onto the patch for all cell description types
  solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

  UpdateResult result;
  switch ( solverPatch.getType() ) {
  case SolverPatch::Type::Cell:
    if (CompressionAccuracy>0.0) {
      uncompress(solverPatch,cellDescriptionsIndex);
    }

    // the actual computations
    updateSolution(solverPatch,cellDescriptionsIndex,solverPatch.getNeighbourMergePerformed(),true);
    result._timeStepSize    = startNewTimeStep(solverPatch,cellDescriptionsIndex);
    result._meshUpdateEvent = updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
        solverPatch,cellDescriptionsIndex,solverPatch.getNeighbourMergePerformed());

    if (CompressionAccuracy>0.0) {
      compress(solverPatch,cellDescriptionsIndex,isAtRemoteBoundary);
    }
    break;
  default:
    _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
    result._meshUpdateEvent =
      _solver->evaluateRefinementCriteriaAfterSolutionUpdate(
          solverPatch,solverPatch.getNeighbourMergePerformed());  // must be done by all cell types TODO(Dominic): Clean up
    break;
  }
  return result;
}

void exahype::solvers::LimitingADERDGSolver::uncompress(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex) const {
  _solver->uncompress(solverPatch);
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
    _limiter->uncompress(limiterPatch);
  }
}

void exahype::solvers::LimitingADERDGSolver::compress(
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  compress(solverPatch,cellDescriptionsIndex,isAtRemoteBoundary);
}

void exahype::solvers::LimitingADERDGSolver::compress(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex,
    const bool isAtRemoteBoundary) const {
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    _solver->compress(solverPatch,isAtRemoteBoundary);
    const int limiterElement =
        _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      _limiter->compress(limiterPatch,isAtRemoteBoundary);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustSolutionDuringMeshRefinement(
    const int cellDescriptionsIndex,
    const int element) {
  const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
    AdjustSolutionDuringMeshRefinementJob job(*this,solverPatch,cellDescriptionsIndex,isInitialMeshRefinement);
    peano::datatraversal::TaskSet spawnedSet( job, peano::datatraversal::TaskSet::TaskType::Background  );
  } else {
    adjustSolutionDuringMeshRefinementBody(solverPatch,cellDescriptionsIndex,isInitialMeshRefinement);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateSolution(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
    const bool backupPreviousSolution) {
  // 1. Erase old cells; now it's safe (TODO(Dominic): Add to docu)
  ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellDescriptionsIndex);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellDescriptionsIndex);

  // 2. Update the solution in the cells
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      assertion(solverPatch.getRefinementStatus()>=-1);
      if (solverPatch.getRefinementStatus()       < _solver->_minimumRefinementStatusForPassiveFVPatch) {
        _solver->updateSolution(solverPatch,neighbourMergePerformed,backupPreviousSolution);
      }
      else if ( solverPatch.getRefinementStatus() < _solver->_minimumRefinementStatusForActiveFVPatch ) {
        _solver->updateSolution(solverPatch,neighbourMergePerformed,backupPreviousSolution);

        LimiterPatch& limiterPatch =
            getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
      else { // solverPatch.getRefinementStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
        _limiter->updateSolution(limiterPatch,cellDescriptionsIndex,backupPreviousSolution);
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        projectFVSolutionOnDGSpace(solverPatch,limiterPatch); // TODO(Dominic): Required for healing
      }
    } else {
      _solver->updateSolution(solverPatch,neighbourMergePerformed,backupPreviousSolution);
    }
  }
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::LimitingADERDGSolver::updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed
) {
  MeshUpdateEvent meshUpdateEvent = MeshUpdateEvent::None;
  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    if (
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch()
    ) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
    } // else: Keep the previously computed min and max values
    
    meshUpdateEvent =
        determineRefinementStatusAfterSolutionUpdate(solverPatch,neighbourMergePerformed);

    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellDescriptionsIndex);
    bool limiterPatchAllocated =
        meshUpdateEvent==MeshUpdateEvent::None &&
        ensureRequiredLimiterPatchIsAllocated(
            solverPatch,cellDescriptionsIndex,solverPatch.getRefinementStatus());
    if ( limiterPatchAllocated ) { // TODO(Dominic): Unsafe
      assertion(solverPatch.getLevel()==getMaximumAdaptiveMeshLevel());
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }
  } else {
    _solver->updateRefinementStatus(solverPatch,neighbourMergePerformed);
    ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellDescriptionsIndex);
  }
  return meshUpdateEvent;
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::LimitingADERDGSolver::determineRefinementStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed) {
  assertion1(solverPatch.getType()==SolverPatch::Type::Cell,solverPatch.toString());

  MeshUpdateEvent meshUpdateEvent = MeshUpdateEvent::None;
  bool dmpViolated = !evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch);
  bool padViolated = !evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found
  bool isTroubled = dmpViolated || padViolated;
  if ( isTroubled ) {
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
    if ( solverPatch.getRefinementStatus() > _solver->_minimumRefinementStatusForActiveFVPatch ) {
      meshUpdateEvent = MeshUpdateEvent::None;
      solverPatch.setRefinementStatus(_solver->_minimumRefinementStatusForTroubledCell);
    } else {
      meshUpdateEvent = MeshUpdateEvent::IrregularLimiterDomainChange;
      solverPatch.setRefinementStatus(_solver->_minimumRefinementStatusForTroubledCell);
      //logInfo("determineLimiterStatusAfterSolutionUpdate()","irregular for x="<<solverPatch.getOffset() << ", level="<<solverPatch.getLevel() << "status="<<solverPatch.getRefinementStatus()<<","<<solverPatch.getPreviousLimiterStatus()<<","<<solverPatch.getExternalLimiterStatus()<<",max status="<<max status );
    }
    if (solverPatch.getLevel()<getMaximumAdaptiveMeshLevel()) {
      //logInfo("determineRefinementStatusAfterSolutionUpdate(...)","troubled on coarse grid. dmpViolated="<<dmpViolated<<". padViolated="<<padViolated<<". cell="<<solverPatch.toString());
      meshUpdateEvent = MeshUpdateEvent::RefinementRequested;
    }
  } else { // we cool the troubled cells down slowly
    if ( solverPatch.getPreviousRefinementStatus() >= _solver->getMinimumRefinementStatusForTroubledCell() ) {
      solverPatch.setRefinementStatus(_solver->getMinimumRefinementStatusForTroubledCell());
      solverPatch.setIterationsToCureTroubledCell(
          solverPatch.getIterationsToCureTroubledCell()-1);
      if (solverPatch.getIterationsToCureTroubledCell()==0) {
        solverPatch.setRefinementStatus(_solver->getMinimumRefinementStatusForTroubledCell()-1);
        solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
      }
    } else {
      _solver->updateRefinementStatus(solverPatch,neighbourMergePerformed); // update the limiter status first
      if ( solverPatch.getRefinementStatus() < _solver->_minimumRefinementStatusForPassiveFVPatch ) { // if below threshold evaluate the refinement criterion
        meshUpdateEvent = Solver::mergeMeshUpdateEvents(
            meshUpdateEvent,
            _solver->evaluateRefinementCriteriaAfterSolutionUpdate(
                solverPatch,solverPatch.getNeighbourMergePerformed()));
      }
      if ( solverPatch.getRefinementStatus() - solverPatch.getPreviousRefinementStatus() > _solver->_limiterHelperLayers ) {
        meshUpdateEvent = Solver::mergeMeshUpdateEvents(
            meshUpdateEvent, MeshUpdateEvent::IrregularLimiterDomainChange );
        //logInfo("determineLimiterStatusAfterSolutionUpdate()","irregular for x="<<solverPatch.getOffset() << ", level="<<solverPatch.getLevel() << "status="<<solverPatch.getRefinementStatus()<<","<<solverPatch.getPreviousLimiterStatus()<<","<<solverPatch.getExternalLimiterStatus()<<",max status="<<max status );
      }
    }
  }

  assertion2(
      !isTroubled ||
      solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForTroubledCell(),
      isTroubled,
      solverPatch.getRefinementStatus());
  return meshUpdateEvent;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = getDataHeapArray(solverPatch.getSolution());

  const int numberOfObservables = _solver->getDMPObservables();
  if ( numberOfObservables>0 ) {
    double* observablesMin = getDataHeapArray(solverPatch.getSolutionMin());
    double* observablesMax = getDataHeapArray(solverPatch.getSolutionMax());

    // 1. Check if the DMP is satisfied and search for the min and max
    // Write the new min and max to the storage reserved for face 0
    bool dmpIsSatisfied = discreteMaximumPrincipleAndMinAndMaxSearch(solution, observablesMin,observablesMax);

    // 2. Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    return dmpIsSatisfied;
  } else {
    return true;
  }
}

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch) {
  double* observablesMin = nullptr;
  double* observablesMax = nullptr;

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables > 0) {
    observablesMin = getDataHeapArray(solverPatch.getSolutionMin());
    observablesMax = getDataHeapArray(solverPatch.getSolutionMax());
  }

  const double* const solution = getDataHeapArray(solverPatch.getSolution());

  return _solver->isPhysicallyAdmissible(
      solution,
      observablesMin,observablesMax,
      solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForTroubledCell,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      solverPatch.getCorrectorTimeStamp(),solverPatch.getCorrectorTimeStepSize());
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if ( solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion1( solverPatch.getRefinementStatus()>=-1,solverPatch.getRefinementStatus() );

      if ( solverPatch.getRefinementStatus()<-1 ) {
        logError("determineMinAndMax(...)","solverPatch.getRefinementStatus()<-1 for cell="<<solverPatch.toString());
        std::abort();
      }

      if (solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch()) {
        determineSolverMinAndMax(solverPatch,true);
      } else { // solverPatch.getRefinementStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
    } else {
      determineSolverMinAndMax(solverPatch,true);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(
    SolverPatch& solverPatch, const bool validate) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),
            solverPatch.toString());
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()),
            solverPatch.toString());

    const double* const solution = getDataHeapArray(solverPatch.getSolution());

    double* observablesMin = getDataHeapArray(solverPatch.getSolutionMin());
    double* observablesMax = getDataHeapArray(solverPatch.getSolutionMax());

    // Write the result to the face with index "0"
    findCellLocalMinAndMax(solution, observablesMin, observablesMax);

    // Copy the result on the other faces as well
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    #ifdef Asserts
    if ( validate ) {
      for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
        assertion2(*(observablesMin+i)<std::numeric_limits<double>::max(),i,solverPatch.toString());
        assertion2(*(observablesMax+i)>-std::numeric_limits<double>::max(),i,solverPatch.toString());
      } // Dead code elimination will get rid of this loop
    }
    #endif
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* limiterSolution = getDataHeapArray(limiterPatch.getSolution());

    double* observablesMin = getDataHeapArray(solverPatch.getSolutionMin());
    double* observablesMax = getDataHeapArray(solverPatch.getSolutionMax());

    // Write the result to the face with index "0"
    findCellLocalLimiterMinAndMax(limiterSolution, observablesMin, observablesMax);

    // Copy the result on the other faces as well
    const int numberOfObservables = _solver->getDMPObservables();
    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
      std::copy_n(
          observablesMin,numberOfObservables, // past-the-end element
          observablesMin+i*numberOfObservables);
      std::copy_n(
          observablesMax,numberOfObservables, // past-the-end element
          observablesMax+i*numberOfObservables);
    }

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion(*(observablesMin+i)<std::numeric_limits<double>::max());
      assertion(*(observablesMax+i)>-std::numeric_limits<double>::max());
    } // Dead code elimination will get rid of this loop
  }
}

void exahype::solvers::LimitingADERDGSolver::deallocateLimiterPatch(
    const SolverPatch& solverPatch,
    const int cellDescriptionsIndex) const {
  const int limiterElement =
        _solver->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  assertion(limiterElement!=Solver::NotFound);
  LimiterPatch& limiterPatch = FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  limiterPatch.setType(LimiterPatch::Type::Erased);
  _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex).erase(
      FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex).begin()+limiterElement);
  lock.free();
}

void exahype::solvers::LimitingADERDGSolver::ensureNoLimiterPatchIsAllocatedOnHelperCell(
    const SolverPatch& solverPatch,
    const int cellDescriptionsIndex) const {
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (
      limiterElement!=exahype::solvers::Solver::NotFound &&
      solverPatch.getType()!=SolverPatch::Type::Cell
  ) {
    deallocateLimiterPatch(solverPatch,cellDescriptionsIndex);
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
    const SolverPatch& solverPatch,
    const int cellDescriptionsIndex) const {
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (
      limiterElement!=Solver::NotFound      &&
      solverPatch.getRefinementStatus()         < _solver->_minimumRefinementStatusForPassiveFVPatch &&
      solverPatch.getPreviousRefinementStatus() < _solver->_minimumRefinementStatusForPassiveFVPatch
  ) {
    deallocateLimiterPatch(solverPatch,cellDescriptionsIndex);
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustLimiterSolution(
    SolverPatch& solverPatch,
    LimiterPatch& limiterPatch) const {
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  _limiter->adjustSolution(limiterPatch);
}

void exahype::solvers::LimitingADERDGSolver::allocateLimiterPatch(
        const SolverPatch& solverPatch,
        const int cellDescriptionsIndex) const {
  #if defined(Asserts)
  const int previousLimiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  #endif
  assertion(previousLimiterElement==Solver::NotFound);

  exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
      cellDescriptionsIndex,
      solverPatch.getSolverNumber(),
      LimiterPatch::Type::Cell,
      LimiterPatch::RefinementEvent::None,
      solverPatch.getLevel(),
      solverPatch.getParentIndex(),
      solverPatch.getSize(),
      solverPatch.getOffset());

  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolution()),solverPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),solverPatch.toString());

  LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
  _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);
}

bool exahype::solvers::LimitingADERDGSolver::ensureRequiredLimiterPatchIsAllocated(
        const SolverPatch& solverPatch,
        const int cellDescriptionsIndex,
        const int limiterStatus) const {
  const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      limiterElement==Solver::NotFound                      &&
      solverPatch.getType()==SolverPatch::Type::Cell        &&
      limiterStatus>=_solver->_minimumRefinementStatusForPassiveFVPatch
  ) {
    assertion1(solverPatch.getPreviousRefinementStatus()<_solver->_minimumRefinementStatusForPassiveFVPatch,
               solverPatch.toString());
    allocateLimiterPatch(solverPatch,cellDescriptionsIndex);
    return true;
  }
  return false;
}


void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = getDataHeapArray(solverPatch.getSolution());
  double*       limiterSolution = getDataHeapArray(limiterPatch.getSolution());

  projectOnFVLimiterSpace(solverSolution, limiterSolution);
}

// TODO(Dominic): Check that we have rolled back in time as well
void exahype::solvers::LimitingADERDGSolver::rollbackSolutionGlobally(
    const int cellDescriptionsIndex, const int solverElement,
    const bool fusedTimeStepping) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  if ( fusedTimeStepping ) { // TODO merge with synchronisation
    rollbackToPreviousTimeStepFused(cellDescriptionsIndex,solverElement);
  } else {
    rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
  }

  // 1. Ensure limiter patch is allocated (based on previous limiter status
  ensureRequiredLimiterPatchIsAllocated(
      solverPatch,cellDescriptionsIndex,
      solverPatch.getPreviousRefinementStatus());

  // 2. Rollback solution to previous time step
  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    if ( solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

      _solver->swapSolutionAndPreviousSolution(solverPatch);   // roll back solver

      if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch ) {
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch); // roll back limiter (must exist!)
      } else {
        const int limiterElement =
            tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
        if (limiterElement!=Solver::NotFound) {
          LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
          copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
          projectDGSolutionOnFVSpace(solverPatch,limiterPatch); // project DG solution on patch
        }
      }
    }
    else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      _solver->swapSolutionAndPreviousSolution(solverPatch);
    }
  }

  // 3. Reset the previous refinement status.
  solverPatch.setRefinementStatus(solverPatch.getPreviousRefinementStatus());
}

void exahype::solvers::LimitingADERDGSolver::rollbackSolutionLocally(
    const int cellDescriptionsIndex,
    const int solverElement,
    const bool fusedTimeStepping) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  if (fusedTimeStepping) { // TODO merge with synchronisation
    rollbackToPreviousTimeStepFused(cellDescriptionsIndex,solverElement);
  } else {
    rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
  }

  // 1. Ensure limiter patch is allocated (based on current limiter status)
  ensureRequiredLimiterPatchIsAllocated(
      solverPatch,cellDescriptionsIndex,
      solverPatch.getRefinementStatus());

  // 2. Now roll back to the last valid solution
  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForPassiveFVPatch
  ) { // this is one of the important differences to the global recomputation where we rollback also cells with limiter status == 0
    assertion(solverPatch.getType()==SolverPatch::Type::Cell);
    assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);
    if ( solverPatch.getType()!=SolverPatch::Type::Cell ) {
      logError("rollbackSolutionLocally(..)","type is not cell but is=" << solverPatch.toString());
      std::abort();
    }

    if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch ) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
    }
    else { // We need to project limiter data for the previous stage
      _solver->swapSolutionAndPreviousSolution(solverPatch);
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }

    // 2.1 Reset the iterationsToCure on all troubled cells to maximum value if cell is troubled
    if (solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForTroubledCell()) {
      solverPatch.setIterationsToCureTroubledCell(1+_iterationsToCureTroubledCell);
    }
  }

  // 3. Only after the reinitialisation, it is safe to deallocate the limiter patch
  ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellDescriptionsIndex);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellDescriptionsIndex);
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = getDataHeapArray(limiterPatch.getSolution());
  double*       solverSolution  = getDataHeapArray(solverPatch.getSolution());

  projectOnDGSpace(limiterSolution, solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
    SolverPatch& solverPatch,
    const int cellDescriptionsIndex) {
  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForPassiveFVPatch
  ) {
    if (
        solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch()
    ) { // these guys are recomputing with the limiter
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(
          solverPatch,cellDescriptionsIndex);

      _limiter->updateSolution(limiterPatch,cellDescriptionsIndex,true);
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    }
    else { // these guys are just swapping and projecting
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex);

      if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch ) { // these did just do a swap
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      }
      else { // this one has a new FV patch
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
    }
  }
  // limiterStatus==0 or cell is not on finest level
#if defined(Asserts)
  const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  assertion1(
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() || limiterElement==Solver::NotFound,
      solverPatch.toString());
#endif
}

double exahype::solvers::LimitingADERDGSolver::recomputeSolutionLocally(
        const int cellDescriptionsIndex, const int solverElement) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  recomputeSolution(solverPatch,cellDescriptionsIndex);

  return startNewTimeStep(solverPatch,cellDescriptionsIndex);
}

double exahype::solvers::LimitingADERDGSolver::recomputeSolutionLocallyFused(
        const int cellDescriptionsIndex,
        const int element,
        const bool isAtRemoteBoundary) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  recomputeSolution(solverPatch,cellDescriptionsIndex);

  double admissibleTimeStepSize =
      startNewTimeStepFused(solverPatch,cellDescriptionsIndex,true,true);

  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getRefinementStatus() < _solver->getMinimumRefinementStatusForTroubledCell()  &&
      solverPatch.getPredictorTimeStepSize() > 0
  ) {
    if (
        solverPatch.getRefinementStatus() >= _solver->getMinimumRefinementStatusForActiveFVPatch()
        || // TODO(Dominic): Reassess what this is doing
        (solverPatch.getRefinementStatus() < _solver->getMinimumRefinementStatusForActiveFVPatch() &&
            solverPatch.getPreviousRefinementStatus()>=_solver->getMinimumRefinementStatusForTroubledCell())
    ) {
      _solver->performPredictionAndVolumeIntegral(
          solverPatch,cellDescriptionsIndex,element,
          solverPatch.getCorrectorTimeStamp(), // corrector time step data is correct; see docu
          solverPatch.getCorrectorTimeStepSize(),
          false/*already uncompressed*/,isAtRemoteBoundary);
    }
  }

  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::prolongateFaceData(
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) {
  _solver->prolongateFaceData(cellDescriptionsIndex,element,isAtRemoteBoundary);
}

void exahype::solvers::LimitingADERDGSolver::restriction(
        const int cellDescriptionsIndex,
        const int element) {
  _solver->restriction(cellDescriptionsIndex,element);
}


///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursData(
    const int                                  solverNumber,
    Solver::CellInfo&                       context1,
    Solver::CellInfo&                       context2,
    const tarch::la::Vector<DIMENSIONS, int>&  pos1,
    const tarch::la::Vector<DIMENSIONS, int>&  pos2,
    const bool                                 isRecomputation) const {
  const int solverElement1 = indexOfCellDescription(context1._ADERDGCellDescriptions,solverNumber);
  const int solverElement2 = indexOfCellDescription(context2._ADERDGCellDescriptions,solverNumber);

  if ( solverElement1!=Solver::NotFound && solverElement2!=Solver::NotFound) {
    SolverPatch& solverPatch1 = context1._ADERDGCellDescriptions[solverElement1];
    SolverPatch& solverPatch2 = context2._ADERDGCellDescriptions[solverElement2];

    const int limiterElement1 = indexOfCellDescription(context1._FiniteVolumesCellDescriptions,solverNumber); // might be NotFound
    const int limiterElement2 = indexOfCellDescription(context2._FiniteVolumesCellDescriptions,solverNumber);

    synchroniseTimeStepping(solverPatch1,context1._FiniteVolumesCellDescriptions,limiterElement1);
    synchroniseTimeStepping(solverPatch2,context2._FiniteVolumesCellDescriptions,limiterElement2);

    //
    // 1. Solve the Riemann problems/copy boundary layers
    //

    // We only limit on the finest mesh level
    if ( solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion2(solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());
      // 1.1. Merge solver solution or limiter solution values in
      // non-overlapping parts of solver and limiter domain:
      if (solverPatch1.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch() &&
          solverPatch2.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch()) {
        if (!isRecomputation) {
          _solver->mergeNeighboursData(solverNumber,context1,context2,pos1,pos2);
           // Which one is left and right is checked internally again.
        }
      }
      else if (solverPatch1.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch() &&
               solverPatch2.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch()) {
        _limiter->mergeNeighboursData(solverNumber,context1,context2,pos1,pos2);
      }
      // 1.2. Merge limiter solution values in overlapping part
      // of solver and limiter domain:
      else if (
          (solverPatch1.getRefinementStatus() >= _solver->getMinimumRefinementStatusForActiveFVPatch() &&
           solverPatch2.getRefinementStatus() <  _solver->getMinimumRefinementStatusForActiveFVPatch() &&
           solverPatch2.getRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch)
          ||
          (solverPatch1.getRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch    &&
           solverPatch1.getRefinementStatus() <  _solver->getMinimumRefinementStatusForActiveFVPatch() &&
           solverPatch2.getRefinementStatus() >= _solver->getMinimumRefinementStatusForActiveFVPatch())
      ) {
        assertion2(limiterElement1!=Solver::NotFound,solverPatch1.toString(),solverPatch2.toString());
        assertion2(limiterElement2!=Solver::NotFound,solverPatch2.toString(),solverPatch1.toString());
        _limiter->mergeNeighboursData(solverNumber,context1,context2,pos1,pos2);
        if (!isRecomputation) {
          _solver->mergeNeighboursData(solverNumber,context1,context2,pos1,pos2);
             // Which one is left and right is checked internally again.
        }
      }
      else {
        logError("mergeNeighboursBasedOnLimiterStatus(...)","Neighbours cannot communicate. " <<
            std::endl << "cell1=" << solverPatch1.toString() <<
            std::endl << ".cell2=" << solverPatch2.toString());
        std::terminate();
      }

    // On the other levels, we work with the ADER-DG solver only
    } else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      if ( !isRecomputation ) {
        _solver->mergeNeighboursData(solverNumber,context1,context2,pos1,pos2);
            // Which one is left and right is checked internally again.
      }
    }

    // 2. Merge the min and max of both cell description's solver's
    // solution value. This is not done during recomputations.

    if (
        !isRecomputation &&
        _solver->getDMPObservables() > 0
    ) {
      Solver::InterfaceInfo face(pos1,pos2);
      mergeSolutionMinMaxOnFace(solverPatch1,solverPatch2,face);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  Solver::InterfaceInfo& face) const {
  if (
      (solverPatch1.getCommunicationStatus()==ADERDGSolver::CellCommunicationStatus &&
      solverPatch1.getFacewiseCommunicationStatus(face._faceIndex1) >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
      solverPatch1.getFacewiseAugmentationStatus(face._faceIndex1)  <  ADERDGSolver::MaximumAugmentationStatus) // excludes Ancestors
      ||
      (solverPatch2.getCommunicationStatus()==ADERDGSolver::CellCommunicationStatus &&
      solverPatch2.getFacewiseCommunicationStatus(face._faceIndex2) >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
      solverPatch2.getFacewiseAugmentationStatus(face._faceIndex2)  <  ADERDGSolver::MaximumAugmentationStatus) // excludes Ancestors
  ) {
    assertion( solverPatch1.getSolverNumber() == solverPatch2.getSolverNumber() );
    const int numberOfObservables = _solver->getDMPObservables();
    double* min1 = getDataHeapArray(solverPatch1.getSolutionMin()) + face._faceIndex1 * numberOfObservables;
    double* min2 = getDataHeapArray(solverPatch2.getSolutionMin()) + face._faceIndex2 * numberOfObservables;
    double* max1 = getDataHeapArray(solverPatch1.getSolutionMax()) + face._faceIndex1 * numberOfObservables;
    double* max2 = getDataHeapArray(solverPatch2.getSolutionMax()) + face._faceIndex2 * numberOfObservables;

    for (int i=0; i<numberOfObservables; i++) {
      const double min = std::min(
          *(min1+i),
          *(min2+i)
      );
      const double max = std::max(
          *(max1+i),
          *(max2+i)
      );

      *(min1+i) = min;
      *(min2+i) = min;

      *(max1+i) = max;
      *(max2+i) = max;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
    const int                                 solverNumber,
    Solver::CellInfo&                         context,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
    const bool                                isRecomputation) {
  assertion2(tarch::la::countEqualEntries(posCell,posBoundary)==(DIMENSIONS-1),posCell.toString(),posBoundary.toString());
  Solver::BoundaryFaceInfo face(posCell,posBoundary);

  const int solverElement = indexOfCellDescription(context._ADERDGCellDescriptions,solverNumber);
  if ( solverElement != Solver::NotFound ) {
    SolverPatch& solverPatch = context._ADERDGCellDescriptions[solverElement];

    if ( !solverPatch.getNeighbourMergePerformed(face._faceIndex) ) { // just check. no reset necessary as it is done by sub solvers
      const int limiterElement = indexOfCellDescription(context._FiniteVolumesCellDescriptions,solverNumber);

      synchroniseTimeStepping(solverPatch,context._FiniteVolumesCellDescriptions,limiterElement);

      if (solverPatch.getType()==SolverPatch::Type::Cell) {
        if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
          if (solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch()) {
            assertion(solverPatch.getRefinementStatus()<=_solver->_minimumRefinementStatusForPassiveFVPatch || limiterElement!=Solver::NotFound);
            if (!isRecomputation) {
              _solver->mergeWithBoundaryData(solverNumber,context,posCell,posBoundary);
            }
          }
          else {
            assertion(limiterElement!=Solver::NotFound);
            _limiter->mergeWithBoundaryData(solverNumber,context,posCell,posBoundary);
          }
        }
        else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
          if (!isRecomputation) {
            _solver->mergeWithBoundaryData(solverNumber,context,posCell,posBoundary);
          }
        }
      }
    }
  }
}

#ifdef Parallel
///////////////////////////////////
// NEIGHBOUR - Mesh refinement
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  _solver->appendNeighbourCommunicationMetadata(
      metadata,src,dest,
      cellDescriptionsIndex,solverNumber);
}

///////////////////////////////////
// NEIGHBOUR - Time marching
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbour(
    const int                                    toRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const int element = indexOfCellDescription(cellInfo._ADERDGCellDescriptions,solverNumber);
  if ( element != NotFound ) {
    sendMinAndMaxToNeighbour(toRank,cellInfo._ADERDGCellDescriptions[element],src,dest,x,level);

    sendDataToNeighbourBasedOnLimiterStatus(
        toRank,solverNumber,cellInfo,src,dest,
        false,/* isRecomputation */
        x,level);
  } else {
    sendEmptyDataToNeighbour(toRank,x,level);
  }

  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
}

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                    toRank,
    const SolverPatch&                           solverPatch,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  if ( numberOfObservables>0 ) {
    if (ADERDGSolver::holdsFaceData(solverPatch)) {
      Solver::BoundaryFaceInfo& face(src,dest);
      const double* observablesMin = getDataHeapArrayFacePart(solverPatch.getSolutionMin(),numberOfObservables,face._faceIndex);
      const double* observablesMax = getDataHeapArrayFacePart(solverPatch.getSolutionMax(),numberOfObservables,face._faceIndex);

      DataHeap::getInstance().sendData(
          observablesMin, numberOfObservables, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().sendData(
          observablesMax, numberOfObservables, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
    } else {
      for(int sends=0; sends<2; ++sends) {
        #if defined(UsePeanosSymmetricBoundaryExchanger)
        DataHeap::getInstance().sendData(
            _invalidObservables.data(), _invalidObservables.size(),toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        #else
        DataHeap::getInstance().sendData(
            exahype::EmptyDataHeapMessage, toRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
        #endif
      }
    }
  }
}


void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    solverNumber,
        CellInfo&                                    cellInfo,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const {
  const int element = indexOfCellDescription(cellInfo._ADERDGCellDescriptions,solverNumber);
  if ( element != NotFound ) {
    if ( level==getMaximumAdaptiveMeshLevel() ) {
      logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " to rank="<<toRank<<",x="<<x<<",level="<<level);

      // solver sends
      if (
          solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForTroubledCell() &&
          isRecomputation==false
      ) {
        _solver->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
      } else {
        _solver->sendEmptyDataToNeighbour(toRank,x,level);
      }

      // limiter sends (receive order must be inverted)
      if (solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForPassiveFVPatch) {
        const int limiterElement = indexOfCellDescription(cellInfo._FiniteVolumesCellDescriptions,solverNumber);
        if ( limiterElement!=NotFound ) {
          _limiter->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
        } else { // if the limiter status of a cell changes dramatically, a limiter patch might not been allocated
          // at the time data is sent to neighbouring ranks if fused time stepping is used.
          assertion1(Solver::FuseADERDGPhases,solverPatch.toString());
          _limiter->sendEmptyDataToNeighbour(toRank,x,level);
        }
      } else {
        _limiter->sendEmptyDataToNeighbour(toRank,x,level);
      }

    } else {
      if ( isRecomputation==false ) {
        _solver->sendDataToNeighbour(toRank,solverNumber,cellInfo,element,src,dest,x,level);
      } else {
        _solver->sendEmptyDataToNeighbour(toRank,x,level);
      }
    }
  } else { // solver element not found
    _solver->sendEmptyDataToNeighbour(toRank,x,level);
    if ( level==getMaximumAdaptiveMeshLevel() ) {
      _limiter->sendEmptyDataToNeighbour(toRank,x,level);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // send an empty minAndMax message
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables > 0) {
    for(int sends=0; sends<2; ++sends) {
      #if defined(UsePeanosSymmetricBoundaryExchanger)
      DataHeap::getInstance().sendData(
          _invalidObservables.data(), _invalidObservables.size(), toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      #else
      DataHeap::getInstance().sendData(
          exahype::EmptyDataHeapMessage, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      #endif
    }
  }

  _solver->sendEmptyDataToNeighbour(toRank,x,level);
  if (level==getMaximumAdaptiveMeshLevel()) {
    _limiter->sendEmptyDataToNeighbour(toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
          fromRank << " at vertex x=" << x << ", level=" << level <<
          ", source=" << src << ", destination=" << dest);

  mergeWithNeighbourDataBasedOnLimiterStatus(
      fromRank,cellDescriptionsIndex,element,src,dest,
      false,/*isRecomputation*/x,level);

  mergeWithNeighbourMinAndMax(fromRank,cellDescriptionsIndex,element,src,dest,x,level);

  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataBasedOnLimiterStatus(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const bool                                   isRecomputation,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  synchroniseTimeStepping(cellDescriptionsIndex,element);

  if ( level==getMaximumAdaptiveMeshLevel() ) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

    logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
        fromRank << " at vertex x=" << x << ", level=" << level <<
        ", source=" << src << ", destination=" << dest << ",limiterStatus=" << solverPatch.getRefinementStatus());
    assertion1(solverPatch.getRefinementStatus()>=ADERDGSolver::Pending,solverPatch.toString());

    if ( solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch() ) {
      _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      if (!isRecomputation) {
        _solver->mergeWithNeighbourData(
            fromRank,cellDescriptionsIndex,element,
            src,dest,x,level);
      } else {
        _solver->dropNeighbourData(fromRank,src,dest,x,level);
      }
    } else { // solverPatch.getRefinementStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=Solver::NotFound);
      _limiter->mergeWithNeighbourData(
          fromRank,cellDescriptionsIndex,limiterElement,
          src,dest,x,level);
      _solver->dropNeighbourData(fromRank,src,dest,x,level);
    }

  } else {

    logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
        fromRank << " at vertex x=" << x << ", level=" << level <<
        ", source=" << src << ", destination=" << dest);

    if (!isRecomputation) {
      _solver->mergeWithNeighbourData(
          fromRank,cellDescriptionsIndex,element,
          src,dest,x,level);
    } else {
      _solver->dropNeighbourData(fromRank,src,dest,x,level);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMinAndMax(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

    Solver::BoundaryFaceInfo face(dest,src);
    if(
        (solverPatch.getCommunicationStatus()                       == ADERDGSolver::CellCommunicationStatus &&
        solverPatch.getFacewiseCommunicationStatus(face._faceIndex) >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
        solverPatch.getFacewiseAugmentationStatus(face._faceIndex)  <  ADERDGSolver::MaximumAugmentationStatus)
        ||
        (solverPatch.getFacewiseCommunicationStatus(face._faceIndex) == ADERDGSolver::CellCommunicationStatus &&
        solverPatch.getCommunicationStatus()                         >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
        solverPatch.getAugmentationStatus()                          <  ADERDGSolver::MaximumAugmentationStatus)
      ){
      // Inverted send-receive order: TODO(Dominic): Add to docu
      // Send order:    min,max
      // Receive order; max,min
      DataHeap::getInstance().receiveData(
          const_cast<double*>(_receivedMax.data()), numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().receiveData(
          const_cast<double*>(_receivedMin.data()), numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);

      mergeSolutionMinMaxOnFace(
          solverPatch,face,_receivedMin.data(),_receivedMax.data());
    } else {
      for(int receives=0; receives<2; ++receives)
        DataHeap::getInstance().receiveData(
            fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  solverPatch,
  Solver::BoundaryFaceInfo& face,
  const double* const min,
  const double* const max) const {
  assertion1(ADERDGSolver::holdsFaceData(solverPatch),solverPatch.toString());

  double* solutionMin = getDataHeapArray(solverPatch.getSolutionMin());
  double* solutionMax = getDataHeapArray(solverPatch.getSolutionMax());

  const int numberOfObservables = _solver->getDMPObservables();
  for (int i=0; i<numberOfObservables; i++) {
    solutionMin[i+face._faceIndex*numberOfObservables]  = std::min( solutionMin[i+face._faceIndex*numberOfObservables], min[i] );
    solutionMax[i+face._faceIndex*numberOfObservables]  = std::max( solutionMax[i+face._faceIndex*numberOfObservables], max[i] );
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
  if ( level==getMaximumAdaptiveMeshLevel() ) {
    _limiter->dropNeighbourData(fromRank,src,dest,x,level);
  }
  _solver->dropNeighbourData(fromRank,src,dest,x,level);

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    for(int receives=0; receives<2; ++receives)
      DataHeap::getInstance().receiveData(
          fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
  }
}

///////////////////////////////////////////////////////////
// NEIGHBOUR - Solution recomputation
///////////////////////////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendEmptySolverAndLimiterDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToNeighbour(toRank,x,level);
  if ( level==getMaximumAdaptiveMeshLevel() ) {
    _limiter->sendEmptyDataToNeighbour(toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourSolverAndLimiterData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  logDebug("dropNeighbourSolverAndLimiterData(...)", "drop data for solver " << _identifier << " from rank " <<
            fromRank << " at vertex x=" << x << ", level=" << level <<
            ", src=" << src << ", dest=" << dest);

  if ( level==getMaximumAdaptiveMeshLevel() ) {
    _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  }
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

/////////////////////////////////////
// MASTER<=>WORKER
/////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  _solver->progressMeshRefinementInPrepareSendToWorker(
      workerRank, fineGridCell, fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell, coarseGridVerticesEnumerator,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerIfProlongating(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToWorkerIfProlongating(
      workerRank,cellDescriptionsIndex,element,x,level);
}

void exahype::solvers::LimitingADERDGSolver::receiveDataFromMasterIfProlongating(
    const int masterRank,
    const int receivedCellDescriptionsIndex,
    const int receivedElement,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  _solver->receiveDataFromMasterIfProlongating(
      masterRank,receivedCellDescriptionsIndex,receivedElement,x,level);
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement) {
 return _solver->progressMeshRefinementInMergeWithWorker(
      localCellDescriptionsIndex,
      receivedCellDescriptionsIndex,receivedElement);
}

void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInPrepareSendToMaster(
    const int masterRank,
    const int cellDescriptionsIndex, const int element,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  _solver->progressMeshRefinementInPrepareSendToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int  level,
    const bool stillInRefiningMode) {
  return _solver->progressMeshRefinementInMergeWithMaster(
      worker,localCellDescriptionsIndex,localElement,coarseGridCellDescriptionsIndex,x,level,stillInRefiningMode);
}

void exahype::solvers::LimitingADERDGSolver::appendMasterWorkerCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  _solver->appendMasterWorkerCommunicationMetadata(
      metadata,cellDescriptionsIndex,solverNumber);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToWorkerOrMasterDueToForkOrJoin(
      toRank,cellDescriptionsIndex,element,messageType,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=Solver::NotFound) {
    _limiter->sendDataToWorkerOrMasterDueToForkOrJoin(
        toRank,cellDescriptionsIndex,limiterElement,messageType,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
      fromRank,cellDescriptionsIndex,element,messageType,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=Solver::NotFound) {
    _limiter->mergeWithWorkerOrMasterDataDueToForkOrJoin(
        fromRank,cellDescriptionsIndex,limiterElement,messageType,x,level);
  }
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  _solver->sendDataToMaster(masterRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  _solver->sendDataToWorker(workerRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithMasterData(masterRank,x,level);
}
#endif

std::string exahype::solvers::LimitingADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::LimitingADERDGSolver::toString (std::ostream& out) const {
  out << getIdentifier() << "{_ADERDG: ";
  out << _solver->toString() << "}\n";
  out << getIdentifier() << "{_FV: ";
  out << _limiter->toString() << "}";
}


exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::FusedTimeStepJob(
  LimitingADERDGSolver& solver,
  SolverPatch&          solverPatch,
  const int             cellDescriptionsIndex,
  const int             element,
  const bool            isSkeletonJob):
  _solver(solver),
  _solverPatch(solverPatch),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _neighbourMergePerformed(solverPatch.getNeighbourMergePerformed()),
  _isSkeletonJob(isSkeletonJob) {
  // copy the neighbour merge performed array
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::operator()() {
  _solver.fusedTimeStepBody(
      _solverPatch,_cellDescriptionsIndex,_element,
      false,false,_isSkeletonJob,false,_neighbourMergePerformed);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}

exahype::solvers::LimitingADERDGSolver::AdjustSolutionDuringMeshRefinementJob::AdjustSolutionDuringMeshRefinementJob(
  LimitingADERDGSolver& solver,
  SolverPatch&          solverPatch,
  const int             cellDescriptionsIndex,
  const bool            isInitialMeshRefinement):
  _solver(solver),
  _solverPatch(solverPatch),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _isInitialMeshRefinement(isInitialMeshRefinement)
{
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::AdjustSolutionDuringMeshRefinementJob::operator()() {
  _solver.adjustSolutionDuringMeshRefinementBody(_solverPatch,_cellDescriptionsIndex,_isInitialMeshRefinement);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
