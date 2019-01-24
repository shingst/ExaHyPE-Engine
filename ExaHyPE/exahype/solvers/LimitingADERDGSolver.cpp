/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
  N 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 *
 * \author Dominic E. Charrier, Tobias Weinzierl
 **/

#include <algorithm> // copy_n
#include <iomanip>
#include <chrono>

#include "LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"
#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "kernels/finitevolumes/commons/c/commons.h" // TODO measurements

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
    SolverPatch& solverPatch,
    Solver::CellInfo& cellInfo) const {
  _solver->synchroniseTimeStepping(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  _solver->startNewTimeStep();
  ensureLimiterTimeStepDataIsConsistent();

  logDebug("startNewTimeStep()","getMeshUpdateEvent()="<<Solver::toString(getMeshUpdateEvent())<<
           ",getNextMeshUpdateEvent()="<<Solver::toString(getNextMeshUpdateEvent()));
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    const bool isFirstTimeStepOfBatch,
    const bool isLastTimeStepOfBatch) {
  _solver->startNewTimeStepFused(isFirstTimeStepOfBatch,isLastTimeStepOfBatch);
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

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatch(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  assertion1(limiterElement!=NotFound,solverPatch.toString());
  LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return cellInfo._FiniteVolumesCellDescriptions[limiterElement];
}

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatch(
    const SolverPatch& solverPatch,CellInfo& cellInfo,const int limiterElement) const {
  assertion1(limiterElement!=NotFound,solverPatch.toString());
  LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return cellInfo._FiniteVolumesCellDescriptions[limiterElement];
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
    const int solverNumber,
    CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement!=NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    if (
        solverPatch.getType()==SolverPatch::Type::Cell &&
        solverPatch.getRefinementStatus()<ADERDGSolver::Erase
    ) {
      logError("determineMinAndMax(...)","solverPatch.getRefinementStatus()<-1 for cell="<<solverPatch.toString());
      std::abort();
    }
    // only done when doing initial refinement
    const bool newLimiterPatchAllocated =
        getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested &&
        ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    if ( newLimiterPatchAllocated ) {
      assertion1(tarch::la::equals(solverPatch.getCorrectorTimeStamp(),0.0),solverPatch.toString());
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
      adjustLimiterSolution(solverPatch,limiterPatch);
    }
    ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);

    if ( getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested ) {
      solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());
    }

    // must come last because of assertion in patch allocation
    _solver->finaliseStateUpdates(solverNumber,cellInfo);
  }
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(SolverPatch& solverPatch,CellInfo& cellInfo) {
  double admissibleTimeStepSize = _solver->startNewTimeStep(solverPatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
  return admissibleTimeStepSize;
}

double exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    SolverPatch& solverPatch,
    CellInfo&    cellInfo,
    const bool   isFirstTimeStepOfBatch,
    const bool   isLastTimeStepOfBatch)  {
  double admissibleTimeStepSize =_solver->startNewTimeStepFused(solverPatch,isFirstTimeStepOfBatch,isLastTimeStepOfBatch);
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
  return admissibleTimeStepSize;
}

double exahype::solvers::LimitingADERDGSolver::updateTimeStepSizes(
    const int solverNumber,CellInfo& cellInfo,const bool fused) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != Solver::NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    double admissibleTimeStepSize = std::numeric_limits<double>::max();
    admissibleTimeStepSize = _solver->updateTimeStepSizes(solverPatch,fused);
    ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
    return admissibleTimeStepSize;
  }
  else {
    return std::numeric_limits<double>::max();
  }
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    SolverPatch& solverPatch,CellInfo& cellInfo,const bool fused) const {
  synchroniseTimeStepping(solverPatch,cellInfo); // TODO(Dominic): Correct?
  if ( fused ) {
    _solver->rollbackToPreviousTimeStepFused(solverPatch);
  } else {
    _solver->rollbackToPreviousTimeStep(solverPatch);
  }
  ensureLimiterPatchTimeStepDataIsConsistent(solverPatch,cellInfo);
}

void exahype::solvers::LimitingADERDGSolver::adjustSolutionDuringMeshRefinementBody(
    SolverPatch& solverPatch,
    CellInfo& cellInfo,
    const bool isInitialMeshRefinement) {
  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    if (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::Prolongating) {
      _solver->prolongateVolumeData(solverPatch,isInitialMeshRefinement);
      solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
    }
    _solver->adjustSolution(solverPatch);

    determineSolverMinAndMax(solverPatch,false);
    if ( !evaluatePhysicalAdmissibilityCriterion(solverPatch,solverPatch.getCorrectorTimeStamp()) ) {
       solverPatch.setRefinementStatus(_solver->getMinRefinementStatusForTroubledCell());
       solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
    } else {
      _solver->markForRefinement(solverPatch);
    }
  }

  solverPatch.setHasCompletedLastStep(true);
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if ( limiterElement!=Solver::NotFound ) {
    LimiterPatch& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[limiterElement];
//    assertion2(solverPatch.getPreviousRefinementStatus() >=_solver->getMinRefinementStatusForSeparationCell() ||
//               solverPatch.getRefinementStatus() >=_solver->getMinRefinementStatusForSeparationCell(),solverPatch.toString(),limiterPatch.toString());
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

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStepOrRestrict(
    const int  solverNumber,
    CellInfo&  cellInfo,
    const bool isFirstTimeStepOfBatch,
    const bool isLastTimeStepOfBatch,
    const bool isAtRemoteBoundary) {
  const int element        = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[element];
    synchroniseTimeStepping(solverPatch,cellInfo);
    solverPatch.setHasCompletedLastStep(false);

    // Write the previous limiter status back onto the patch for all cell description types
    solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

    if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
      const bool isAMRSkeletonCell     = solverPatch.getHasVirtualChildren();
      const bool isSkeletonCell        = isAMRSkeletonCell || isAtRemoteBoundary;
      const bool mustBeDoneImmediately = isSkeletonCell && PredictionSweeps==1;

      if (
          (SpawnUpdateAsBackgroundJob || (SpawnPredictionAsBackgroundJob && !isLastTimeStepOfBatch)) &&
          !mustBeDoneImmediately
      ) {
        peano::datatraversal::TaskSet( new FusedTimeStepJob(
            *this,solverPatch,cellInfo,isFirstTimeStepOfBatch,isLastTimeStepOfBatch,isSkeletonCell) );
        return UpdateResult();
      } else {
        return fusedTimeStepBody(
            solverPatch, cellInfo,
            solverPatch.getNeighbourMergePerformed(),
            isFirstTimeStepOfBatch,isLastTimeStepOfBatch,
            isSkeletonCell,
            mustBeDoneImmediately);
      }
    }
    else { // other cell types
      UpdateResult result;
      if (
          solverPatch.getType()==SolverPatch::Type::Descendant &&
          solverPatch.getCommunicationStatus()>=ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication
      ) {
        _solver->restrictToTopMostParent(solverPatch);
      }
      ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
      _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
      result._meshUpdateEvent = _solver->evaluateRefinementCriteriaAfterSolutionUpdate(
          solverPatch,solverPatch.getNeighbourMergePerformed()); // must be done by all cell types TODO(Dominic): Clean up

      solverPatch.setHasCompletedLastStep(true);
      return result;
    }
  } else {
    return UpdateResult();
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStepBody(
    SolverPatch&                                               solverPatch,
    CellInfo&                                                  cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
    const bool                                                 isFirstTimeStepOfBatch,
    const bool                                                 isLastTimeStepOfBatch,
    const bool                                                 isSkeletonCell,
    const bool                                                 mustBeDoneImmediately) {
  UpdateResult result;
  updateSolution(solverPatch,cellInfo,neighbourMergePerformed,isFirstTimeStepOfBatch);
  const bool isTroubled = checkIfCellIsTroubledAndDetermineMinAndMax(solverPatch,cellInfo);
  revisitSolverPatchesInBuffer(solverPatch,cellInfo,isTroubled,neighbourMergePerformed,true);
  result._timeStepSize    = startNewTimeStep(solverPatch,cellInfo);
  result._meshUpdateEvent = determineRefinementStatusAfterSolutionUpdate(solverPatch,cellInfo,isTroubled,neighbourMergePerformed);

  if (
      solverPatch.getRefinementStatus()<_solver->getMinRefinementStatusForTroubledCell() &&
      SpawnPredictionAsBackgroundJob &&
      !mustBeDoneImmediately         &&
      isLastTimeStepOfBatch
  ) {
    const int element = cellInfo.indexOfADERDGCellDescription(solverPatch.getSolverNumber());
    peano::datatraversal::TaskSet( new ADERDGSolver::PredictionJob(
        *_solver.get(),solverPatch/*the reductions are delegated to _solver anyway*/,
        cellInfo._cellDescriptionsIndex,element,
        solverPatch.getCorrectorTimeStamp(),  // corrector time step data is correct; see docu
        solverPatch.getCorrectorTimeStepSize(),false/*is uncompressed*/,isSkeletonCell ));
  }
  else if ( solverPatch.getRefinementStatus()<_solver->getMinRefinementStatusForTroubledCell() ){
    _solver->performPredictionAndVolumeIntegralBody(
        solverPatch,
        solverPatch.getCorrectorTimeStamp(),  // corrector time step data is correct; see docu
        solverPatch.getCorrectorTimeStepSize(),
        false, isSkeletonCell );
  } else {
    solverPatch.setHasCompletedLastStep(true);
  }

  return result;
}

void exahype::solvers::LimitingADERDGSolver::performPredictionAndVolumeIntegral(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    SolverPatch& cellDescription = cellInfo._ADERDGCellDescriptions[element];

    if ( cellDescription.getRefinementStatus()<_solver->getMinRefinementStatusForTroubledCell() ) {
      _solver->performPredictionAndVolumeIntegral(solverNumber,cellInfo,isAtRemoteBoundary);
    }
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::updateBody(
    SolverPatch&                                               solverPatch,
    CellInfo&                                                  cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
    const bool                                                 isAtRemoteBoundary){
  if (CompressionAccuracy>0.0) { uncompress(solverPatch,cellInfo); }

  // the actual computations
  UpdateResult result;
  updateSolution(solverPatch,cellInfo,neighbourMergePerformed,true);
  const bool isTroubled = checkIfCellIsTroubledAndDetermineMinAndMax(solverPatch,cellInfo);
  revisitSolverPatchesInBuffer(solverPatch,cellInfo,isTroubled,neighbourMergePerformed,true);
  result._timeStepSize    = startNewTimeStep(solverPatch,cellInfo); // uses DG solution to compute time step size; might be result of FV->DG projection
  result._meshUpdateEvent = determineRefinementStatusAfterSolutionUpdate(solverPatch,cellInfo,isTroubled,neighbourMergePerformed);

  if (CompressionAccuracy>0.0) { compress(solverPatch,cellInfo,isAtRemoteBoundary); }

  solverPatch.setHasCompletedLastStep(true); // required as prediction checks the flag too. Field should be renamed "setHasCompletedLastOperation(...)".
  return result;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::updateOrRestrict(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool isAtRemoteBoundary){
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    synchroniseTimeStepping(solverPatch,cellInfo);
    solverPatch.setHasCompletedLastStep(false);

    if ( solverPatch.getType()==SolverPatch::Type::Cell && SpawnUpdateAsBackgroundJob ) {
      peano::datatraversal::TaskSet( new UpdateJob(
          *this, solverPatch,cellInfo,isAtRemoteBoundary ) );
      return UpdateResult();
    }
    else if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
      return updateBody(solverPatch,cellInfo,solverPatch.getNeighbourMergePerformed(),isAtRemoteBoundary);
    }
    else { // other cell types
      UpdateResult result;
      if (
          solverPatch.getType()==SolverPatch::Type::Descendant &&
          solverPatch.getCommunicationStatus()>=ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication
      ) {
        _solver->restrictToTopMostParent(solverPatch);
      }
      ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
      _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
      result._meshUpdateEvent = _solver->evaluateRefinementCriteriaAfterSolutionUpdate(
          solverPatch,solverPatch.getNeighbourMergePerformed()); // must be done by all cell types
      solverPatch.setHasCompletedLastStep(true);
      return result;
    }
  } else {
    return UpdateResult();
  }
}

void exahype::solvers::LimitingADERDGSolver::uncompress(
    SolverPatch& solverPatch,CellInfo& cellInfo) const {
  _solver->uncompress(solverPatch);
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo,limiterElement);
    _limiter->uncompress(limiterPatch);
  }
}

void exahype::solvers::LimitingADERDGSolver::compress(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    compress(solverPatch,cellInfo,isAtRemoteBoundary);
  }
}

void exahype::solvers::LimitingADERDGSolver::compress(
    SolverPatch& solverPatch,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) const {
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    _solver->compress(solverPatch,isAtRemoteBoundary);
    const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
    if ( limiterElement!=exahype::solvers::Solver::NotFound ) {
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo,limiterElement);
      _limiter->compress(limiterPatch,isAtRemoteBoundary);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustSolutionDuringMeshRefinement(
    const int solverNumber,
    CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    synchroniseTimeStepping(solverPatch,cellInfo);

    const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
    if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
      solverPatch.setHasCompletedLastStep(false);
      peano::datatraversal::TaskSet(
        new AdjustSolutionDuringMeshRefinementJob(*this,solverPatch,cellInfo,isInitialMeshRefinement) );
    } else {
      adjustSolutionDuringMeshRefinementBody(solverPatch,cellInfo,isInitialMeshRefinement);
    }
  }
}

int exahype::solvers::LimitingADERDGSolver::getMaxiumRefinementStatusInNeighbourhood(
    SolverPatch& solverPatch,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed) const {
  int mergedLimiterStatus = solverPatch.getRefinementStatus(); // allows increasing the status if we are "suddenly" are next to FV solver
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    if ( neighbourMergePerformed[i] ) {
      mergedLimiterStatus = std::max(mergedLimiterStatus,solverPatch.getFacewiseRefinementStatus(i)-1);
    }
  }
  return mergedLimiterStatus;
}

void exahype::solvers::LimitingADERDGSolver::updateSolution(
    SolverPatch&                                               solverPatch,
    CellInfo&                                                  cellInfo,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
    const bool                                                 backupPreviousSolution) {
  if (solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
    assertion(solverPatch.getRefinementStatus()>=ADERDGSolver::Erase);
    const int mergedLimiterStatus = getMaxiumRefinementStatusInNeighbourhood(solverPatch,neighbourMergePerformed);

    // update
    if ( mergedLimiterStatus >= _solver->_minRefinementStatusForTroubledCell-1 ) { // limiter update
      assertion1(solverPatch.getRefinementStatus()>0,solverPatch.toString());
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);

      _limiter->updateSolution(limiterPatch,neighbourMergePerformed,cellInfo._cellDescriptionsIndex,backupPreviousSolution);
      _solver->swapSolutionAndPreviousSolution(solverPatch);
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch); // TODO(Dominic): If we do static limiting, this is not necessary in troubled cells
    }
    else {
      _solver->updateSolution(solverPatch,neighbourMergePerformed,backupPreviousSolution);

      if (
          mergedLimiterStatus >=_solver->getMinRefinementStatusForSeparationCell() &&
          solverPatch.getRefinementStatus() < _solver->getMinRefinementStatusForBufferCell()
      ) {
        if ( solverPatch.getRefinementStatus()<_solver->getMinRefinementStatusForSeparationCell() ) {
          ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,mergedLimiterStatus);
        }
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
      else if ( mergedLimiterStatus < _solver->getMinRefinementStatusForSeparationCell() ) {
        ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);
      }
    }
  } else {
    _solver->updateSolution(solverPatch,neighbourMergePerformed,backupPreviousSolution);
  }
}

bool
exahype::solvers::LimitingADERDGSolver::checkIfCellIsTroubledAndDetermineMinAndMax(
    SolverPatch& solverPatch,
    CellInfo&    cellInfo) {
  bool isTroubled =
      !evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch) ||
      !evaluatePhysicalAdmissibilityCriterion(solverPatch,
          solverPatch.getCorrectorTimeStamp()+solverPatch.getCorrectorTimeStepSize()); // after min and max was found

  if ( // above call computes DG min and max on-the-fly. We use the FV min and max if the cell is troubled
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell
  ) {
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    determineLimiterMinAndMax(solverPatch,limiterPatch);
  }

  return isTroubled;
}

void exahype::solvers::LimitingADERDGSolver::revisitSolverPatchesInBuffer(
    SolverPatch&                                               solverPatch,
    CellInfo&                                                  cellInfo,
    const bool                                                 isTroubled,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed,
    const bool                                                 backupPreviousSolution) {
  if (solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
    const int mergedLimiterStatus = getMaxiumRefinementStatusInNeighbourhood(solverPatch,neighbourMergePerformed);

    const bool isADERDGCellInBuffer =
        mergedLimiterStatus                < _solver->getMinRefinementStatusForTroubledCell()-1 &&
        solverPatch.getRefinementStatus() >= _solver->getMinRefinementStatusForBufferCell();
    if ( isADERDGCellInBuffer && isTroubled ) { // limiter update
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
      _limiter->updateSolution(limiterPatch,neighbourMergePerformed,cellInfo._cellDescriptionsIndex,backupPreviousSolution);
      _solver->swapSolutionAndPreviousSolution(solverPatch);
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
    }
    else if ( isADERDGCellInBuffer ) { // just project onto limiter solution
      LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }
  }
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::LimitingADERDGSolver::determineRefinementStatusAfterSolutionUpdate(
    SolverPatch&                                               solverPatch,
    CellInfo&                                                  cellInfo,
    const bool                                                 isTroubled,
    const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed) {
  assertion1(solverPatch.getType()==SolverPatch::Type::Cell,solverPatch.toString());

  // pre-update mesh update events
  MeshUpdateEvent meshUpdateEvent = MeshUpdateEvent::None;
  if (
      isTroubled &&
      (solverPatch.getRefinementStatus()        < _solver->getMinRefinementStatusForBufferCell() ||
      solverPatch.getPreviousRefinementStatus() < _solver->getMinRefinementStatusForBufferCell())
  ) {
    meshUpdateEvent = MeshUpdateEvent::IrregularLimiterDomainChange;
  } else if (
      isTroubled &&
      solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
  ) {
    meshUpdateEvent = MeshUpdateEvent::RefinementRequested;
  }

  // update refinement status
  solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

  _solver->updateRefinementStatus(solverPatch,neighbourMergePerformed);
  if ( isTroubled ) {
    solverPatch.setRefinementStatus(_solver->_minRefinementStatusForTroubledCell);
  }
  else if ( // cool down
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getPreviousRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell
  ) {
    int newRefinementStatus = _solver->_minRefinementStatusForTroubledCell-2; // allows increasing the status if we are "suddenly" are next to FV solver
    for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) { // TODO TODO TODO
      if ( neighbourMergePerformed[i] ) { // TODO TODO TODO
        newRefinementStatus = std::max(newRefinementStatus,solverPatch.getFacewiseRefinementStatus(i)-1); // TODO TODO TODO
      }
    }
    solverPatch.setRefinementStatus(newRefinementStatus);
  }

  // post-update mesh update events (pure ADER-DG functionality
  if (
      solverPatch.getLevel()<getMaximumAdaptiveMeshLevel() ||
      solverPatch.getRefinementStatus() < _solver->getMinRefinementStatusForSeparationCell()
  ) {
    meshUpdateEvent = Solver::mergeMeshUpdateEvents(
        meshUpdateEvent,
        _solver->evaluateRefinementCriteriaAfterSolutionUpdate(solverPatch,neighbourMergePerformed));
  }

  return meshUpdateEvent;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = static_cast<double*>(solverPatch.getSolution());

  const int numberOfObservables = _solver->getDMPObservables();
  if ( numberOfObservables>0 ) {
    double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax());

    // 1. Check if the DMP is satisfied and search for the min and max
    // Write the new min and max to the storage reserved for face 0
    bool dmpIsSatisfied = discreteMaximumPrincipleAndMinAndMaxSearch(solution, observablesMin,observablesMax);

    // 1. Copy the result on the other faces as well
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

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(
    SolverPatch& solverPatch,const double timeStamp) {
  double* observablesMin = nullptr;
  double* observablesMax = nullptr;

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables > 0) {
    observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    observablesMax = static_cast<double*>(solverPatch.getSolutionMax());
  }

  const double* const solution = static_cast<double*>(solverPatch.getSolution());

  return _solver->isPhysicallyAdmissible(
      solution,
      observablesMin,observablesMax,
      solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      timeStamp);
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int solverNumber,Solver::CellInfo& cellInfo) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    if (
        solverPatch.getType()==SolverPatch::Type::Cell &&
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()
    ) {
      assertion1( solverPatch.getRefinementStatus()>=ADERDGSolver::Erase,solverPatch.getRefinementStatus() );

      if ( solverPatch.getRefinementStatus()<ADERDGSolver::Erase ) {
        logError("determineMinAndMax(...)","solverPatch.getRefinementStatus()<-1 for cell="<<solverPatch.toString());
        std::abort();
      } else if (solverPatch.getRefinementStatus()<_solver->getMinRefinementStatusForTroubledCell()) {
        determineSolverMinAndMax(solverPatch,true);
      } else {
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
    } else if (solverPatch.getType()==SolverPatch::Type::Cell) {
      determineSolverMinAndMax(solverPatch,true);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(
    SolverPatch& solverPatch, const bool validate) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionIndex()),
            solverPatch.toString());
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMinIndex()),
            solverPatch.toString());

    const double* const solution = static_cast<double*>(solverPatch.getSolution());

    double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax());

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
    double* limiterSolution = static_cast<double*>(limiterPatch.getSolution());

    double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin());
    double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax());

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
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  assertion(limiterElement!=Solver::NotFound);
  LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo,limiterElement);
  limiterPatch.setType(LimiterPatch::Type::Erased);
  _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch); // has its own lock
  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  cellInfo._FiniteVolumesCellDescriptions.erase(cellInfo._FiniteVolumesCellDescriptions.begin()+limiterElement);
  lock.free();
}

void exahype::solvers::LimitingADERDGSolver::ensureNoLimiterPatchIsAllocatedOnHelperCell(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (
      limiterElement!=exahype::solvers::Solver::NotFound &&
      solverPatch.getType()!=SolverPatch::Type::Cell
  ) {
    deallocateLimiterPatch(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
    const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (
      limiterElement!=Solver::NotFound      &&
      solverPatch.getRefinementStatus()         < _solver->getMinRefinementStatusForSeparationCell() &&
      solverPatch.getPreviousRefinementStatus() < _solver->getMinRefinementStatusForSeparationCell()
  ) {
    deallocateLimiterPatch(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustLimiterSolution(
    SolverPatch& solverPatch,
    LimiterPatch& limiterPatch) const {
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  _limiter->adjustSolution(limiterPatch);
}

void exahype::solvers::LimitingADERDGSolver::allocateLimiterPatch(const SolverPatch& solverPatch,CellInfo& cellInfo) const {
  #if defined(Asserts)
  const int previousLimiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  #endif
  assertion(previousLimiterElement==Solver::NotFound);

  exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
      solverPatch.getSolverNumber(),
      cellInfo,
      LimiterPatch::Type::Cell,
      LimiterPatch::RefinementEvent::None,
      solverPatch.getLevel(),
      solverPatch.getParentIndex(),
      solverPatch.getSize(),
      solverPatch.getOffset());

  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolutionIndex()),solverPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionIndex()),solverPatch.toString());

  LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
  _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);

  assertion1(DataHeap::getInstance().isValidIndex(limiterPatch.getPreviousSolutionIndex()),limiterPatch.toString());
  assertion1(DataHeap::getInstance().isValidIndex(limiterPatch.getSolutionIndex()),limiterPatch.toString());
  assertion1(limiterPatch.getSolution()!=nullptr,limiterPatch.toString());
  assertion1(limiterPatch.getPreviousSolution()!=nullptr,limiterPatch.toString());
}

bool exahype::solvers::LimitingADERDGSolver::ensureRequiredLimiterPatchIsAllocated(
        const SolverPatch& solverPatch,
        CellInfo&          cellInfo,
        const int          limiterStatus) const {
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  if (
      limiterElement        ==Solver::NotFound              &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType() ==SolverPatch::Type::Cell       &&
      limiterStatus         >= _solver->getMinRefinementStatusForSeparationCell()
  ) {
    assertion1(solverPatch.getPreviousRefinementStatus()<_solver->getMinRefinementStatusForSeparationCell(),
               solverPatch.toString());
    allocateLimiterPatch(solverPatch,cellInfo);
    return true;
  }
  return false;
}


void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = static_cast<double*>(solverPatch.getSolution());
  double*       limiterSolution = static_cast<double*>(limiterPatch.getSolution());

  projectOnFVLimiterSpace(solverSolution, limiterSolution);
}

void exahype::solvers::LimitingADERDGSolver::rollbackSolutionGlobally(
    const int  solverNumber,
    CellInfo&  cellInfo,
    const bool fusedTimeStepping) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement!=NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    rollbackToPreviousTimeStep(solverPatch,cellInfo,fusedTimeStepping);

    // 1. Ensure limiter patch is allocated (based on previous limiter status
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getPreviousRefinementStatus());

    // 2. Roll solution back to previous time step
    if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
      if ( solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() ) {
        assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

        _solver->swapSolutionAndPreviousSolution(solverPatch);   // roll back solver

        if ( solverPatch.getPreviousRefinementStatus() >= _solver->getMinRefinementStatusForSeparationCell() ) {
          LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
          _limiter->swapSolutionAndPreviousSolution(limiterPatch); // roll back limiter (must exist!)
        } else {
          const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
          if ( limiterElement!=Solver::NotFound ) {
            LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
            copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
            projectDGSolutionOnFVSpace(solverPatch,limiterPatch); // project DG solution on patch
          }
        }
      }
      else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
        _solver->swapSolutionAndPreviousSolution(solverPatch);
      }
    }

    // 3. Reset the current to the previous refinement status.
    solverPatch.setRefinementStatus(solverPatch.getPreviousRefinementStatus());
  }
}

void exahype::solvers::LimitingADERDGSolver::rollbackSolutionLocally(
    const int  solverNumber,
    CellInfo&  cellInfo,
    const bool fusedTimeStepping) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    rollbackToPreviousTimeStep(solverPatch,cellInfo,fusedTimeStepping);

    // 1. Ensure limiter patch is allocated (based on current limiter status)
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());

    // 2. Now roll back to the last valid solution
    if (
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2
    ) { // this is one of the important differences to the global recomputation where we rollback also cells with limiter status == 0
      assertion(solverPatch.getType()==SolverPatch::Type::Cell);
      assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);
      if ( solverPatch.getType()!=SolverPatch::Type::Cell ) {
        logError("rollbackSolutionLocally(..)","type is not cell but is=" << solverPatch.toString());
        std::abort();
      }

      if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minRefinementStatusForTroubledCell-1 ) {
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      }
      else { // We need to project limiter data for the previous stage
        LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
    }

    // 3. Only after the reinitialisation, it is safe to deallocate the limiter patch
    ensureNoLimiterPatchIsAllocatedOnHelperCell(solverPatch,cellInfo);
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);
  }
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = static_cast<double*>(limiterPatch.getSolution());
  double*       solverSolution  = static_cast<double*>(solverPatch.getSolution());

  projectOnDGSpace(limiterSolution, solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
    SolverPatch& solverPatch,CellInfo& cellInfo) {
  const bool fineGridCellInvolvedInLocalRecomputation =
      solverPatch.getType()==SolverPatch::Type::Cell        &&
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-2;

  if (
      fineGridCellInvolvedInLocalRecomputation &&
      solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForTroubledCell-1
  ) { // these guys are recomputing with the limiter
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);

    _limiter->updateSolution(limiterPatch,limiterPatch.getNeighbourMergePerformed(),cellInfo._cellDescriptionsIndex,true);
    projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
  }
  else if ( fineGridCellInvolvedInLocalRecomputation ) { // these guys are just swapping and projecting
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    _solver->swapSolutionAndPreviousSolution(solverPatch);
    _limiter->swapSolutionAndPreviousSolution(limiterPatch);
    projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
  }

  // limiterStatus==0 or cell is not on finest level
  #if defined(Asserts)
  const int limiterElement = cellInfo.indexOfFiniteVolumesCellDescription(solverPatch.getSolverNumber());
  assertion1(
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() || limiterElement==Solver::NotFound,
      solverPatch.toString());
  #endif
}

double exahype::solvers::LimitingADERDGSolver::recomputeSolutionLocally(
    const int solverNumber,Solver::CellInfo& cellInfo,
    const bool isAtRemoteBoundary,
    const bool fusedTimeStepping) {
  const int solverElement  = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement!=NotFound ) {
    // 1. Perform the local recomputation
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    recomputeSolution(solverPatch,cellInfo);

    // 1.1. Compute the predictor if the fused scheme is used
    double admissibleTimeStepSize = std::numeric_limits<double>::max();
    if ( fusedTimeStepping ) {
      admissibleTimeStepSize = startNewTimeStepFused(solverPatch,cellInfo,true,true); // TODO(Dominic): Check
      if (
          solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()                                     &&
          solverPatch.getType()==SolverPatch::Type::Cell                                            &&
          solverPatch.getPredictorTimeStepSize() > 0 // TODO(Dominic): Why is the dt > 0 necessary? Historic?
          &&
          solverPatch.getRefinementStatus() < _solver->getMinRefinementStatusForTroubledCell()
          &&                                                                                                  // is not troubled and
          (solverPatch.getRefinementStatus()        == _solver->getMinRefinementStatusForTroubledCell()-1 ||  // holds an active FV patch                                                                                                 // or
          solverPatch.getPreviousRefinementStatus() >= _solver->getMinRefinementStatusForTroubledCell())      // was previously troubled but is no more
      ) {
        _solver->performPredictionAndVolumeIntegral(
            solverNumber,cellInfo,
            solverPatch.getCorrectorTimeStamp(), // corrector time step data is correct; see docu
            solverPatch.getCorrectorTimeStepSize(),
            false/*already uncompressed*/,isAtRemoteBoundary);
      }
    } else {
      admissibleTimeStepSize = startNewTimeStep(solverPatch,cellInfo);
    }

    // 2. Set the previous limiter status to the current one
    solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

    return admissibleTimeStepSize;
  } else {
    return std::numeric_limits<double>::max();
  }
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursDataDuringLocalRecomputation(
    const int                                  solverNumber,
    Solver::CellInfo&                          cellInfo1,
    Solver::CellInfo&                          cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>&  pos1,
    const tarch::la::Vector<DIMENSIONS, int>&  pos2) {
  const int solverElement1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int solverElement2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);

  if ( solverElement1!=Solver::NotFound && solverElement2!=Solver::NotFound) {
    SolverPatch& solverPatch1 = cellInfo1._ADERDGCellDescriptions[solverElement1];
    SolverPatch& solverPatch2 = cellInfo2._ADERDGCellDescriptions[solverElement2];

    if (solverPatch1.getType()==SolverPatch::Type::Cell &&
        solverPatch2.getType()==SolverPatch::Type::Cell &&
        solverPatch1.getRefinementStatus()>=_solver->getMinRefinementStatusForTroubledCell()-2 &&
        solverPatch2.getRefinementStatus()>=_solver->getMinRefinementStatusForTroubledCell()-2) {
      assertion2(solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());
      assertion2(solverPatch2.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());

      waitUntilCompletedLastStep<SolverPatch>(solverPatch1,false,false); // must come before any other operation
      waitUntilCompletedLastStep<SolverPatch>(solverPatch2,false,false);

      _limiter->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataDuringLocalRecomputation(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  assertion2(tarch::la::countEqualEntries(posCell,posBoundary)==(DIMENSIONS-1),posCell.toString(),posBoundary.toString());
  Solver::BoundaryFaceInfo face(posCell,posBoundary);

  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != Solver::NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    if (solverPatch.getType()==SolverPatch::Type::Cell &&
        solverPatch.getRefinementStatus()>=_solver->getMinRefinementStatusForTroubledCell()-1) {
      assertion1(solverPatch.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch.toString());

      waitUntilCompletedLastStep<SolverPatch>(solverPatch,false,false); // must come before any other operation

      _limiter->mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighboursData(
    const int                                  solverNumber,
    Solver::CellInfo&                          cellInfo1,
    Solver::CellInfo&                          cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>&  pos1,
    const tarch::la::Vector<DIMENSIONS, int>&  pos2) {
  const int solverElement1 = cellInfo1.indexOfADERDGCellDescription(solverNumber);
  const int solverElement2 = cellInfo2.indexOfADERDGCellDescription(solverNumber);

  if ( solverElement1!=Solver::NotFound && solverElement2!=Solver::NotFound) {
    SolverPatch& solverPatch1 = cellInfo1._ADERDGCellDescriptions[solverElement1];
    SolverPatch& solverPatch2 = cellInfo2._ADERDGCellDescriptions[solverElement2];

    waitUntilCompletedLastStep<SolverPatch>(solverPatch1,false,false); // must come before any other operation
    waitUntilCompletedLastStep<SolverPatch>(solverPatch2,false,false);

    //
    // 1. Solve the Riemann problems/copy boundary layers
    //
    // We only limit on the finest mesh level
    if ( solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion2(solverPatch2.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());

      if (solverPatch1.getRefinementStatus()<_solver->getMinRefinementStatusForTroubledCell() &&
          solverPatch2.getRefinementStatus()<_solver->getMinRefinementStatusForTroubledCell()) {
        // assumes that face fluxes are not directly added to cell update or the update is cleared if FV update is performe
        _solver->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
      }

      if (solverPatch1.getRefinementStatus()>=_solver->getMinRefinementStatusForSeparationCell() &&
          solverPatch2.getRefinementStatus()>=_solver->getMinRefinementStatusForSeparationCell()) {
        _limiter->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
      }

      const int minStatus = std::min(solverPatch1.getRefinementStatus(),solverPatch2.getRefinementStatus());
      const int maxStatus = std::max(solverPatch1.getRefinementStatus(),solverPatch2.getRefinementStatus());
      if ( minStatus<=_solver->getMinRefinementStatusForSeparationCell() &&
           maxStatus>=_solver->getMinRefinementStatusForTroubledCell() ) {
        logError("mergeNeighboursBasedOnLimiterStatus(...)","Neighbours cannot communicate. " <<
            std::endl << "cell1=" << solverPatch1.toString() <<
            std::endl << ".cell2=" << solverPatch2.toString());
        std::terminate();
      }
    // On the other levels, we work with the ADER-DG solver only
    } else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      _solver->mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
    }

    // 2. Merge the min and max of both cell description's solver's
    // solution value. This is not done during recomputations.
    if ( _solver->getDMPObservables() > 0 ) {
      Solver::InterfaceInfo face(pos1,pos2);
      mergeSolutionMinMaxOnFace(solverPatch1,solverPatch2,face);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  Solver::InterfaceInfo& face) const {
  if ( ADERDGSolver::communicateWithNeighbour(solverPatch1,face._faceIndex1) ) {
    assertion1( ADERDGSolver::communicateWithNeighbour(solverPatch2,face._faceIndex2),solverPatch2.toString() );
    assertion( solverPatch1.getSolverNumber() == solverPatch2.getSolverNumber() );
    const int numberOfObservables = _solver->getDMPObservables();
    double* min1 = static_cast<double*>(solverPatch1.getSolutionMin()) + face._faceIndex1 * numberOfObservables;
    double* min2 = static_cast<double*>(solverPatch2.getSolutionMin()) + face._faceIndex2 * numberOfObservables;
    double* max1 = static_cast<double*>(solverPatch1.getSolutionMax()) + face._faceIndex1 * numberOfObservables;
    double* max2 = static_cast<double*>(solverPatch2.getSolutionMax()) + face._faceIndex2 * numberOfObservables;

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
    Solver::CellInfo&                         cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  assertion2(tarch::la::countEqualEntries(posCell,posBoundary)==(DIMENSIONS-1),posCell.toString(),posBoundary.toString());
  Solver::BoundaryFaceInfo face(posCell,posBoundary);

  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != Solver::NotFound ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    waitUntilCompletedLastStep<SolverPatch>(solverPatch,false,false); // must come before any other operation

    if (
        solverPatch.getType()==SolverPatch::Type::Cell &&
        (solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel() || // if then
        solverPatch.getRefinementStatus()<_solver->_minRefinementStatusForTroubledCell) // must be able to switch to DG on the fly
    ) {
      _solver->mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary);
    }

    if (
        solverPatch.getType()==SolverPatch::Type::Cell &&
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getRefinementStatus()>=_solver->_minRefinementStatusForSeparationCell
    ) { // must be able to switch to FV on the fly
      _limiter->mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary);
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

/////////////////////////////////////
// NEIGHBOUR - Local Recomputation
/////////////////////////////////////

void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourDuringLocalRecomputation(
        const int                                    toRank,
        const int                                    solverNumber,
        CellInfo&                                    cellInfo,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(src,dest);
  if (
      level==getMaximumAdaptiveMeshLevel() &&
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    waitUntilCompletedLastStep<SolverPatch>(solverPatch,true,true); // must come before any other operation

    logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " to rank="<<toRank<<",x="<<x<<",level="<<level);

    if ( solverPatch.getRefinementStatus() >= _solver->getMinRefinementStatusForTroubledCell()-2 ) {
      _limiter->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
    } else {
      _limiter->sendEmptyDataToNeighbour(toRank,x,level);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataDuringLocalRecomputation(
    const int                                    fromRank,
    const int                                    solverNumber,
    CellInfo&                                    cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(dest,src); // ! order of arguments
  if (
      level==getMaximumAdaptiveMeshLevel() &&
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " to rank="<<fromRank<<",x="<<x<<",level="<<level);

    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];

    waitUntilCompletedLastStep<SolverPatch>(solverPatch,true,true); // ! must come before any other operation

    if ( solverPatch.getRefinementStatus() >= _solver->getMinRefinementStatusForTroubledCell()-1 ) {
      _limiter->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
    } else {
      _limiter->dropNeighbourData(fromRank,x,level);
    }
  }
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
  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( element != NotFound ) {
    sendMinAndMaxToNeighbour(toRank,cellInfo._ADERDGCellDescriptions[element],src,dest,x,level);

    sendDataToNeighbourBasedOnLimiterStatus(
        toRank,solverNumber,cellInfo,src,dest,x,level);
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
  BoundaryFaceInfo face(src,dest);
  if ( numberOfObservables>0 && ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex) ) {
    Solver::BoundaryFaceInfo face(src,dest);
    if( ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex) ){
      assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMinIndex()));
      assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMaxIndex()));
      const double* observablesMin = static_cast<double*>(solverPatch.getSolutionMin()) + (face._faceIndex * numberOfObservables);
      const double* observablesMax = static_cast<double*>(solverPatch.getSolutionMax()) + (face._faceIndex * numberOfObservables);

      DataHeap::getInstance().sendData(
          observablesMin, numberOfObservables, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().sendData(
          observablesMax, numberOfObservables, toRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    solverNumber,
        CellInfo&                                    cellInfo,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(src,dest);
  if (
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
    waitUntilCompletedLastStep<SolverPatch>(solverPatch,true,true); // must come before any other operation

    // solver sends
    _solver->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);

    // limiter sends (receive order must be inverted)
    if ( level==getMaximumAdaptiveMeshLevel() ) {
      logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " to rank="<<toRank<<",x="<<x<<",level="<<level);

      if ( solverPatch.getRefinementStatus()>=_solver->getMinRefinementStatusForSeparationCell() ) {
        _limiter->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
      } else {
        _limiter->sendEmptyDataToNeighbour(toRank,x,level);
      }
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    CellInfo&                                    cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  logDebug("mergeWithNeighbourData(...)", "receive for solver " << _identifier <<" from rank="<<fromRank<<"x="<<x<<",level="<<level);

  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  if ( solverElement != NotFound ) {
    mergeWithNeighbourDataBasedOnLimiterStatus(
        fromRank,solverNumber,cellInfo,src,dest,x,level);

    mergeWithNeighbourMinAndMax(fromRank,cellInfo._ADERDGCellDescriptions[solverElement],src,dest,x,level);
  }

  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataBasedOnLimiterStatus(
    const int                                    fromRank,
    const int                                    solverNumber,
    CellInfo&                                    cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(dest,src);
  if (
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    if ( level == getMaximumAdaptiveMeshLevel() ) {
      SolverPatch& solverPatch = cellInfo._ADERDGCellDescriptions[solverElement];
      assertion1(solverPatch.getRefinementStatus()>=ADERDGSolver::Pending,solverPatch.toString());

      // limiter !!! Receive order must be inverted in neighbour comm.
      if (
          solverPatch.getRefinementStatus()                        >= _solver->getMinRefinementStatusForSeparationCell() &&
          solverPatch.getFacewiseRefinementStatus(face._faceIndex) >= _solver->getMinRefinementStatusForSeparationCell() // requires metadata exchange
      ) {
        assertion1(cellInfo.indexOfFiniteVolumesCellDescription(solverNumber)!=Solver::NotFound,solverPatch.toString());
        _limiter->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
      } else {
        _limiter->dropNeighbourData(fromRank,x,level);
      }

      // solver
      if ( solverPatch.getRefinementStatus() < _solver->getMinRefinementStatusForTroubledCell() ) {
        _solver->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
      }
      else {
        _solver->dropNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
      }
    } else {
      _solver->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMinAndMax(
    const int                                    fromRank,
    SolverPatch&                                 solverPatch,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  BoundaryFaceInfo face(dest,src);
  if (
      numberOfObservables>0 &&
      ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex)
  ) {
    // Inverted send-receive order: TODO(Dominic): Add to docu
    // Send order:    min,max
    // Receive order; max,min
    DataHeap::getInstance().receiveData(
        const_cast<double*>(_receivedMax.data()), numberOfObservables, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    DataHeap::getInstance().receiveData(
        const_cast<double*>(_receivedMin.data()), numberOfObservables, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

    mergeSolutionMinMaxOnFace(solverPatch,face,_receivedMin.data(),_receivedMax.data());
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  solverPatch,
  Solver::BoundaryFaceInfo& face,
  const double* const min,
  const double* const max) const {
  assertion1(ADERDGSolver::ADERDGSolver::communicateWithNeighbour(solverPatch,face._faceIndex) ,solverPatch.toString());

  double* solutionMin = static_cast<double*>(solverPatch.getSolutionMin());
  double* solutionMax = static_cast<double*>(solverPatch.getSolutionMax());

  const int numberOfObservables = _solver->getDMPObservables();
  for (int i=0; i<numberOfObservables; i++) {
    solutionMin[i+face._faceIndex*numberOfObservables]  = std::min( solutionMin[i+face._faceIndex*numberOfObservables], min[i] );
    solutionMax[i+face._faceIndex*numberOfObservables]  = std::max( solutionMax[i+face._faceIndex*numberOfObservables], max[i] );
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                    fromRank,
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  const int solverElement = cellInfo.indexOfADERDGCellDescription(solverNumber);
  BoundaryFaceInfo face(dest,src);
  if (
      solverElement != NotFound &&
      ADERDGSolver::communicateWithNeighbour(cellInfo._ADERDGCellDescriptions[solverElement],face._faceIndex)
  ) {
    // send order:   minAndMax,solver,limiter
    // receive order limiter,solver,minAndMax
    if ( level==getMaximumAdaptiveMeshLevel() ) {
      _limiter->dropNeighbourData(fromRank,x,level);
    }
    _solver->dropNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);

    const int numberOfObservables = _solver->getDMPObservables();
    if (numberOfObservables>0) {
      for(int receives=0; receives<2; ++receives)
        DataHeap::getInstance().receiveData(
            fromRank, x, level,
            peano::heap::MessageType::NeighbourCommunication);
    }
  }
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

#if defined(DistributedStealing)
void exahype::solvers::LimitingADERDGSolver::startStealingManager() {
  _solver->startStealingManager();
}

void exahype::solvers::LimitingADERDGSolver::stopStealingManager() {
  _solver->stopStealingManager();
}
#endif

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

exahype::solvers::Solver::CellProcessingTimes exahype::solvers::LimitingADERDGSolver::measureCellProcessingTimes(const int numberOfRuns) {
  // Setup
  const int cellDescriptionsIndex = ADERDGSolver::Heap::getInstance().createData(0,1);
  FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionsIndex,0,1);

  Solver::CellInfo cellInfo(cellDescriptionsIndex);
  _solver->addNewCellDescription(
      0,cellInfo,SolverPatch::Type::Cell,SolverPatch::RefinementEvent::None,
      getMaximumAdaptiveMeshLevel(), /* needs to be on the fine grid for the limiter cells */-1,
      getCoarsestMeshSize(),
      _domainOffset);

  SolverPatch& solverPatch   = cellInfo._ADERDGCellDescriptions[0];
  _solver->ensureNecessaryMemoryIsAllocated(solverPatch);

  adjustSolutionDuringMeshRefinementBody(solverPatch,cellInfo,true);
  solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
  updateTimeStepSizes(0,cellInfo,false);

  // ADER-DG specific setup ( all Riemanns have been performed, cell is surrounded by other Cell type cells )
  solverPatch.setNeighbourMergePerformed(true);
  solverPatch.setAugmentationStatus(0);
  solverPatch.setFacewiseAugmentationStatus(0);
  solverPatch.setCommunicationStatus(ADERDGSolver::CellCommunicationStatus);
  solverPatch.setFacewiseCommunicationStatus(ADERDGSolver::CellCommunicationStatus);

  // MEASUREMENTS
  CellProcessingTimes result;

  // measure ADERDG STP
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    int numberOfPicardIterations = std::numeric_limits<int>::max();
    for (int it=0; it<numberOfRuns; it++) {
      numberOfPicardIterations = _solver->performPredictionAndVolumeIntegralBody(solverPatch,solverPatch.getPredictorTimeStamp(),solverPatch.getPredictorTimeStepSize(),false,true);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._minTimePredictor = time_sec / numberOfRuns / numberOfPicardIterations;
    result._maxTimePredictor = result._minTimePredictor * _solver->getNodesPerCoordinateAxis(); // * (order+1)
  }

  // measure ADER-DG cells
  {
    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->_refineOrKeepOnFineGrid);
      solverPatch.setFacewiseRefinementStatus(_solver->_refineOrKeepOnFineGrid); // assumed  to be very cheap

      updateBody(solverPatch,cellInfo,solverPatch.getNeighbourMergePerformed(),true);

      _solver->swapSolutionAndPreviousSolution(solverPatch); // assumed  to be very cheap
      _solver->rollbackToPreviousTimeStep(solverPatch);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDGUpdate = time_sec / numberOfRuns;
  }

  // measure ADER-DG -> FV cells
  solverPatch.setRefinementStatus(_solver->getMinRefinementStatusForSeparationCell());
  solverPatch.setFacewiseRefinementStatus(_solver->getMinRefinementStatusForSeparationCell());
  {
    // FV specific setup ( all copies have been performed, impose periodic boundary conditions )
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    adjustLimiterSolution(solverPatch,limiterPatch);
    limiterPatch.setNeighbourMergePerformed(true);
    tarch::la::Vector<DIMENSIONS,int> center(1);
    dfor3(neighbour) // periodic BCs
      if ( tarch::la::countEqualEntries(center,neighbour)==DIMENSIONS-1 ) { // only consider faces
        double* FVSolution = static_cast<double*>(limiterPatch.getSolution());
        _limiter->ghostLayerFilling(FVSolution,FVSolution,neighbour-center);
      }
    enddforx

    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->_refineOrKeepOnFineGrid);
      solverPatch.setFacewiseRefinementStatus(_solver->getMinRefinementStatusForSeparationCell()); // assumed  to be very cheap

      updateBody(solverPatch,cellInfo,solverPatch.getNeighbourMergePerformed(),true);

      _solver->swapSolutionAndPreviousSolution(solverPatch); // assumed  to be very cheap
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      _solver->rollbackToPreviousTimeStep(solverPatch);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeADERDG2FVUpdate = time_sec / numberOfRuns;
  }

  // measure FV -> ADERDG cells
  solverPatch.setRefinementStatus(_solver->getMinRefinementStatusForTroubledCell()-1);
  solverPatch.setFacewiseRefinementStatus(_solver->getMinRefinementStatusForTroubledCell());
  {
    // FV specific setup ( all copies have been performed, impose periodic boundary conditions )
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    adjustLimiterSolution(solverPatch,limiterPatch);
    limiterPatch.setNeighbourMergePerformed(true);
    tarch::la::Vector<DIMENSIONS,int> center(1);
    dfor3(neighbour) // periodic BCs
      if ( tarch::la::countEqualEntries(center,neighbour)==DIMENSIONS-1 ) { // only consider faces
        double* FVSolution = static_cast<double*>(limiterPatch.getSolution());
        _limiter->ghostLayerFilling(FVSolution,FVSolution,neighbour-center);
      }
    enddforx

    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->getMinRefinementStatusForTroubledCell()-1);
      solverPatch.setFacewiseRefinementStatus(_solver->getMinRefinementStatusForTroubledCell()); // assumed  to be very cheap

      updateBody(solverPatch,cellInfo,solverPatch.getNeighbourMergePerformed(),true);

      _solver->swapSolutionAndPreviousSolution(solverPatch);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      _solver->rollbackToPreviousTimeStep(solverPatch); // assumed  to be very cheap
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeFV2ADERDGUpdate = time_sec / numberOfRuns;
  }

  // measure troubled / FV cells
  solverPatch.setRefinementStatus(_solver->getMinRefinementStatusForTroubledCell());
  solverPatch.setFacewiseRefinementStatus(_solver->getMinRefinementStatusForTroubledCell());
  {
    // FV specific setup ( all copies have been performed, impose periodic boundary conditions )
    ensureRequiredLimiterPatchIsAllocated(solverPatch,cellInfo,solverPatch.getRefinementStatus());
    LimiterPatch& limiterPatch = getLimiterPatch(solverPatch,cellInfo);
    adjustLimiterSolution(solverPatch,limiterPatch);
    limiterPatch.setNeighbourMergePerformed(true);
    tarch::la::Vector<DIMENSIONS,int> center(1);
    dfor3(neighbour) // periodic BCs
      if ( tarch::la::countEqualEntries(center,neighbour)==DIMENSIONS-1 ) { // only consider faces
        double* FVSolution = static_cast<double*>(limiterPatch.getSolution());
        _limiter->ghostLayerFilling(FVSolution,FVSolution,neighbour-center);
      }
    enddforx

    const std::chrono::high_resolution_clock::time_point timeStart = std::chrono::high_resolution_clock::now();
    for (int it=0; it<numberOfRuns; it++) {
      solverPatch.setRefinementStatus(_solver->getMinRefinementStatusForTroubledCell());
      solverPatch.setFacewiseRefinementStatus(_solver->getMinRefinementStatusForTroubledCell());

      updateBody(solverPatch,cellInfo,solverPatch.getNeighbourMergePerformed(),true);

      _solver->swapSolutionAndPreviousSolution(solverPatch);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
      _solver->rollbackToPreviousTimeStep(solverPatch);
    }
    const double time_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-timeStart).count() * 1e-9;
    result._timeFVUpdate = time_sec / numberOfRuns;
  }

  // Clean up
  solverPatch.setRefinementStatus(0);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(solverPatch,cellInfo);
  solverPatch.setType(SolverPatch::Type::Erased);
  _solver->ensureNoUnnecessaryMemoryIsAllocated(solverPatch);

  DataHeap::getInstance().deleteAllData();
  ADERDGSolver::Heap::getInstance().deleteAllData();
  FiniteVolumesSolver::Heap::getInstance().deleteAllData();

  return result;
}
