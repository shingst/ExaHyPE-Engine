/**
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

int exahype::solvers::LimitingADERDGSolver::getMaxMinimumLimiterStatusForTroubledCell() {
  int result = 0;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if ( solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG ) {
      result = std::max(
          result,
          static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->
          getMinimumRefinementStatusForTroubledCell());
    }
  }
  return result;
}

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
  assertion(_solver->getNumberOfParameters() == 0);
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

void exahype::solvers::LimitingADERDGSolver::updateNextAttainedStableState(const bool& attainedStableState)  {
  _solver->updateNextAttainedStableState(attainedStableState);
}
bool exahype::solvers::LimitingADERDGSolver::getNextAttainedStableState() const  {
  return _solver->getNextAttainedStableState();
}
bool exahype::solvers::LimitingADERDGSolver::getAttainedStableState() const  {
  return _solver->getAttainedStableState();
}
void exahype::solvers::LimitingADERDGSolver::setNextAttainedStableState()  {
  _solver->setNextAttainedStableState();
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
  _coarsestMeshSize  = coarsestMeshInfo.first;
  _coarsestMeshLevel = coarsestMeshInfo.second;

  _meshUpdateEvent=MeshUpdateEvent::RegularRefinementRequested;

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
      isMergingMetadata = getMeshUpdateEvent()!=MeshUpdateEvent::None;
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

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  _solver->startNewTimeStep();
  ensureLimiterTimeStepDataIsConsistent();

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<static_cast<int>(_meshUpdateEvent)<<
           ",nextLimiterDomainChange="<<static_cast<int>(_nextMeshUpdateEvent));
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  _solver->startNewTimeStepFused(isFirstIterationOfBatch,isLastIterationOfBatch);
  ensureLimiterTimeStepDataIsConsistent();

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<static_cast<int>(_meshUpdateEvent)<<
           ",nextLimiterDomainChange="<<static_cast<int>(_nextMeshUpdateEvent));
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

void exahype::solvers::LimitingADERDGSolver::updateRefinementStatusDuringRefinementStatusSpreading(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  synchroniseTimeStepping(cellDescriptionsIndex,solverElement);
  if ( solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForTroubledCell ) {
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
  }
  _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
  ensureNoLimiterPatchIsAllocatedOnHelperCell(cellDescriptionsIndex,solverElement);
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const bool initialGrid,
    const int solverNumber) {
  return
      _solver->progressMeshRefinementInEnterCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVerticesEnumerator,
          initialGrid,solverNumber);

}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  return
      _solver->progressMeshRefinementInLeaveCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,fineGridPositionOfCell,solverNumber);
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
  const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
  const int solverElement = _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  bool limiterRefinementDone = true;
  if ( solverElement!=exahype::solvers::Solver::NotFound ) {
    SolverPatch& solverPatch =
    _solver->getCellDescription(cellDescriptionsIndex,solverElement);
    limiterRefinementDone = 
         !(solverPatch.getType()==SolverPatch::Type::Descendant  &&
           solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
           (solverPatch.getRefinementStatus()>0                     ||
           solverPatch.getPreviousRefinementStatus()>0));
  } 
  if ( limiterRefinementDone )
    return 
        _solver->attainedStableState(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,solverNumber);
  else return false;
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

    // for global re-computation and mesh refinement
    const bool newLimiterPatchAllocated =
        ensureRequiredLimiterPatchIsAllocated(
            cellDescriptionsIndex,solverElement,solverPatch.getRefinementStatus());
    if (newLimiterPatchAllocated) {
      const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverNumber);
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      adjustLimiterSolution(solverPatch,limiterPatch);
    }
    ensureNoLimiterPatchIsAllocatedOnHelperCell(
        cellDescriptionsIndex,solverElement);
    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
        cellDescriptionsIndex,solverElement);
  }
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement)  {
  double admissibleTimeStepSize =
      _solver->startNewTimeStep(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

  return admissibleTimeStepSize;
}

double exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    const int cellDescriptionsIndex,
    const int solverElement,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch)  {
  double admissibleTimeStepSize =
      _solver->startNewTimeStepFused(cellDescriptionsIndex,solverElement,
                                     isFirstIterationOfBatch,isLastIterationOfBatch);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);

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
    const int cellDescriptionsIndex, const int solverElement) const {
  _solver->zeroTimeStepSizes(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
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
    const int cellDescriptionsIndex,
    const int solverElement,
    const bool isInitialMeshRefinement) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  zeroTimeStepSizes(cellDescriptionsIndex,solverElement);      // TODO(Dominic): Still necessary?
  synchroniseTimeStepping(cellDescriptionsIndex,solverElement);

  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    if (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::Prolongating) {
      _solver->prolongateVolumeData(solverPatch,isInitialMeshRefinement);
      solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
    }
    assertion2(solverPatch.getRefinementEvent()==SolverPatch::None || solverPatch.getRefinementEvent()==SolverPatch::RefiningRequested,solverPatch.toString(),cellDescriptionsIndex);

    _solver->adjustSolution(solverPatch);

    determineSolverMinAndMax(solverPatch);
    if ( !evaluatePhysicalAdmissibilityCriterion(solverPatch) ) {
      solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
      solverPatch.setRefinementStatus(_solver->getMinimumRefinementStatusForTroubledCell());
    } else {
      _solver->markForRefinement(solverPatch); // TODO This code probably overwrites the
      // refinement status during the iterations.
    }

    const int limiterElement =
        tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
    if (limiterElement!=Solver::NotFound) {
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
      _limiter->adjustSolution(limiterPatch);
    } // TODO(Dominic): Add to docu: We adjust the limiter patch here but we do not allocate it.
  }
}

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatchForSolverPatch(
    const int cellDescriptionsIndex, const SolverPatch& solverPatch) const {
  assertion(solverPatch.getRefinementStatus()>=0);
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  return getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex,limiterElement);
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

exahype::solvers::LimitingADERDGSolver::LimiterPatch& exahype::solvers::LimitingADERDGSolver::getLimiterPatchForSolverPatch(
    const SolverPatch& solverPatch, const int cellDescriptionsIndex, const int limiterElement) const {
  assertion(limiterElement!=Solver::NotFound);
  LimiterPatch& limiterPatch =
      FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
  // Ensure time stamps and step sizes are consistent
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  return limiterPatch;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStepBody(
    const int   cellDescriptionsIndex,
    const int   element,
    const bool  isFirstIterationOfBatch,
    const bool  isLastIterationOfBatch,
    const bool  isSkeletonJob,
    const bool  mustBeDoneImmediately,
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed) {
  // synchroniseTimeStepping(cellDescriptionsIndex,element); // assumes this was done in neighbour merge
  updateSolution(cellDescriptionsIndex,element,isFirstIterationOfBatch);
  UpdateResult result;
  result._meshUpdateEvent = updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
                                    cellDescriptionsIndex,element,neighbourMergePerformed);
  // This is important to memorise before calling startNewTimeStepFused
  // TODO(Dominic): Add to docu and/or make cleaner
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  const double memorisedPredictorTimeStamp    = solverPatch.getPredictorTimeStamp();
  const double memorisedPredictorTimeStepSize = solverPatch.getPredictorTimeStepSize();
  result._timeStepSize = startNewTimeStepFused(
      cellDescriptionsIndex,element,isFirstIterationOfBatch,isLastIterationOfBatch);

  if ( solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForTroubledCell() ) {   // TODO(Dominic): Add to docu. This will spawn or do a compression job right afterwards and must thus come last. This order is more natural anyway
    _solver->performPredictionAndVolumeIntegral(
        cellDescriptionsIndex, element,
        memorisedPredictorTimeStamp,memorisedPredictorTimeStepSize,
        false/*already uncompressed*/,isSkeletonJob);
  }
  return result;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::fusedTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isAtRemoteBoundary) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

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
          cellDescriptionsIndex,element,
          isFirstIterationOfBatch,isLastIterationOfBatch,
          isSkeletonCell,
          mustBeDoneImmediately,
          solverPatch.getNeighbourMergePerformed());
    } else {
      FusedTimeStepJob fusedTimeStepJob( *this, cellDescriptionsIndex, element,
          solverPatch.getNeighbourMergePerformed(), isSkeletonCell );
      Solver::submitPredictionJob(fusedTimeStepJob,isSkeletonCell);
      return UpdateResult();
    }
  } else {
    UpdateResult result;
    result._meshUpdateEvent =
        exahype::solvers::Solver::mergeMeshUpdateEvents(
            result._meshUpdateEvent,
            _solver->evaluateRefinementCriteriaAfterSolutionUpdate(
                solverPatch,solverPatch.getNeighbourMergePerformed()));
    return result;
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::update(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary){
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  
  // Write the previous limiter status back onto the patch for all cell description types
  solverPatch.setPreviousRefinementStatus(solverPatch.getRefinementStatus());

  UpdateResult result;
  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
      // uncompress
      if (CompressionAccuracy>0.0) {
        _solver->uncompress(solverPatch);
        const int limiterElement =
            tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
        if (limiterElement!=exahype::solvers::Solver::NotFound) {
          LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
          _limiter->uncompress(limiterPatch);
        }
      }

      // the actual computations
      updateSolution(cellDescriptionsIndex,element,true);
      result._timeStepSize    = startNewTimeStep(cellDescriptionsIndex,element);
      result._meshUpdateEvent = updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
                                       cellDescriptionsIndex,element,solverPatch.getNeighbourMergePerformed());  // !!! limiter status must be updated before refinement criterion is evaluated
      // compress again
      if (CompressionAccuracy>0.0) {
        compress(cellDescriptionsIndex,element,isAtRemoteBoundary);
      }
  } else if ( solverPatch.getType()==SolverPatch::Type::Descendant ) {
    _solver->updateRefinementStatus(solverPatch,solverPatch.getNeighbourMergePerformed());
    result._meshUpdateEvent =
        _solver->evaluateRefinementCriteriaAfterSolutionUpdate(solverPatch,solverPatch.getNeighbourMergePerformed());
  }
  return result;
}

void exahype::solvers::LimitingADERDGSolver::compress(
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    _solver->compress(solverPatch,isAtRemoteBoundary);
    const int limiterElement =
        tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      _limiter->compress(limiterPatch,isAtRemoteBoundary);
    }
  }
}


void exahype::solvers::LimitingADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    const bool backupPreviousSolution) {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  
  // 1. Erase old cells; now it's safe (TODO(Dominic): Add to docu)
  ensureNoLimiterPatchIsAllocatedOnHelperCell(cellDescriptionsIndex,element);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,element);

  // 2. Update the solution in the cells
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(solverPatch.getRefinementStatus()>=-1);

      if (solverPatch.getRefinementStatus()       < _solver->_minimumRefinementStatusForPassiveFVPatch) {
        _solver->updateSolution(solverPatch,backupPreviousSolution);
      }
      else if ( solverPatch.getRefinementStatus() < _solver->_minimumRefinementStatusForActiveFVPatch ) {
        _solver->updateSolution(solverPatch,backupPreviousSolution);

        LimiterPatch& limiterPatch =
            getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
      else { // solverPatch.getRefinementStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        assertion1(limiterElement!=Solver::NotFound,solverPatch.toString());
        _limiter->updateSolution(cellDescriptionsIndex,limiterElement,backupPreviousSolution);

        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        projectFVSolutionOnDGSpace(solverPatch,limiterPatch); // TODO(Dominic): Required for healing
      }
    } else {
      _solver->updateSolution(solverPatch,backupPreviousSolution);
    }
  }
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::LimitingADERDGSolver::updateRefinementStatusAndMinAndMaxAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int solverElement,
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed
) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  MeshUpdateEvent meshUpdateEvent = MeshUpdateEvent::None;
  if ( solverPatch.getType()==SolverPatch::Type::Cell ) {
    const bool solutionIsValid =
        evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch) &&
        evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found
    if (
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch()
    ) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
    } // else: Keep the previously computed min and max values

    meshUpdateEvent =
        determineRefinementStatusAfterSolutionUpdate(solverPatch,!solutionIsValid,neighbourMergePerformed);

    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);

    bool limiterPatchAllocated =
        ensureRequiredLimiterPatchIsAllocated(
            cellDescriptionsIndex,solverPatch.getSolverNumber(),solverPatch.getRefinementStatus());
    if ( limiterPatchAllocated ) {
      assertion(solverPatch.getLevel()==getMaximumAdaptiveMeshLevel());
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);

      copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }
  } else {
    _solver->updateRefinementStatus(solverPatch,neighbourMergePerformed);
    ensureNoLimiterPatchIsAllocatedOnHelperCell(cellDescriptionsIndex,solverElement);
  }
  return meshUpdateEvent;
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::LimitingADERDGSolver::determineRefinementStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,
    const bool isTroubled,
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed) const {
  assertion1(solverPatch.getType()==SolverPatch::Type::Cell,solverPatch.toString());

  MeshUpdateEvent meshUpdateEvent = MeshUpdateEvent::None;
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
      meshUpdateEvent = MeshUpdateEvent::IrregularRefinementRequested;
    }
  } else {
    if ( solverPatch.getPreviousRefinementStatus() >= _solver->getMinimumRefinementStatusForTroubledCell() ) {
      solverPatch.setRefinementStatus(_solver->getMinimumRefinementStatusForTroubledCell());
      solverPatch.setIterationsToCureTroubledCell(
          solverPatch.getIterationsToCureTroubledCell()-1);
      if (solverPatch.getIterationsToCureTroubledCell()==0) {
        solverPatch.setRefinementStatus(_solver->getMinimumRefinementStatusForTroubledCell()-1);
        solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1); // TODO(Dominic): Probably not necessary
      }
    } else {
      _solver->updateRefinementStatus(solverPatch,neighbourMergePerformed);
      if ( solverPatch.getRefinementStatus() < _solver->_minimumRefinementStatusForPassiveFVPatch ) {
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
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();

  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

    // 1. Check if the DMP is satisfied and search for the min and max
    // Write the new min and max to the storage reserved for face 0
    bool dmpIsSatisfied = discreteMaximumPrincipleAndMinAndMaxSearch(solution, observablesMin,observablesMax);

    // TODO(Dominic):
//    // 2. Copy the result on the other faces as well
//    for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
//      std::copy_n(
//          observablesMin,numberOfObservables, // past-the-end element
//          observablesMin+i*numberOfObservables);
//      std::copy_n(
//          observablesMax,numberOfObservables, // past-the-end element
//          observablesMax+i*numberOfObservables);
//    }

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
    observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();
  }

  const double* const solution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

  return _solver->isPhysicallyAdmissible(
      solution,
      observablesMin,observablesMax,numberOfObservables,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      solverPatch.getCorrectorTimeStamp(),solverPatch.getCorrectorTimeStepSize());
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int solverElement) {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      assertion1( solverPatch.getRefinementStatus()>=-1,solverPatch.getRefinementStatus() );
      if (solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch()) {
        determineSolverMinAndMax(solverPatch);
      } else { // solverPatch.getRefinementStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
    } else {
      determineSolverMinAndMax(solverPatch);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(SolverPatch& solverPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),
            solverPatch.toString());
    assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()),
            solverPatch.toString());

    const double* const solution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

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

    for (int i=0; i<DIMENSIONS_TIMES_TWO*numberOfObservables; ++i) {
      assertion2(*(observablesMin+i)<std::numeric_limits<double>::max(),i,solverPatch.toString());
      assertion2(*(observablesMax+i)>-std::numeric_limits<double>::max(),i,solverPatch.toString());
    } // Dead code elimination will get rid of this loop
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    double* limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    double* observablesMin = DataHeap::getInstance().getData(
        solverPatch.getSolutionMin()).data();
    double* observablesMax = DataHeap::getInstance().getData(
        solverPatch.getSolutionMax()).data();

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
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=Solver::NotFound);
  const int limiterElement =
        tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
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
    const int cellDescriptionsIndex,
    const int solverElement) const {
  assertion(solverElement!=Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (
      limiterElement!=exahype::solvers::Solver::NotFound &&
      solverPatch.getType()!=SolverPatch::Type::Cell
  ) {
    deallocateLimiterPatch(cellDescriptionsIndex,solverElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  assertion(solverElement!=Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (
      limiterElement!=Solver::NotFound      &&
      solverPatch.getRefinementStatus()         < _solver->_minimumRefinementStatusForPassiveFVPatch &&
      solverPatch.getPreviousRefinementStatus() < _solver->_minimumRefinementStatusForPassiveFVPatch
  ) {
    deallocateLimiterPatch(cellDescriptionsIndex,solverElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::adjustLimiterSolution(
    SolverPatch& solverPatch,
    LimiterPatch& limiterPatch) const {
  copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
  _limiter->adjustSolution(limiterPatch);
}

int exahype::solvers::LimitingADERDGSolver::allocateLimiterPatch(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  #if defined(Asserts)
  const int previouslimiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  #endif
  assertion(previouslimiterElement==Solver::NotFound);

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

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  assertion(limiterElement!=Solver::NotFound);
  LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(
      solverPatch,cellDescriptionsIndex,limiterElement);
  _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);

  return limiterElement;
}

bool exahype::solvers::LimitingADERDGSolver::ensureRequiredLimiterPatchIsAllocated(
        const int cellDescriptionsIndex,
        const int solverElement,
        const int limiterStatus) const {
  assertion(solverElement!=Solver::NotFound);
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      limiterElement==Solver::NotFound                      &&
      solverPatch.getType()==SolverPatch::Type::Cell        &&
      limiterStatus>=_solver->_minimumRefinementStatusForPassiveFVPatch
  ) {
    assertion1(solverPatch.getPreviousRefinementStatus()<_solver->_minimumRefinementStatusForPassiveFVPatch,
               solverPatch.toString());
    allocateLimiterPatch(cellDescriptionsIndex,solverElement);
    return true;
  }
  return false;
}


void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
  double*       limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();

  projectOnFVLimiterSpace(solverSolution, limiterSolution);
}

// TODO(Dominic): Check that we have rolled back in time as well
void exahype::solvers::LimitingADERDGSolver::rollbackSolutionGlobally(
    const int cellDescriptionsIndex, const int solverElement,
    const bool fusedTimeStepping) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  if (fusedTimeStepping) { // TODO merge with synchronisation
    rollbackToPreviousTimeStepFused(cellDescriptionsIndex,solverElement);
  } else {
    rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
  }

  // 1. Ensure limiter patch is allocated (based on previous limiter status
  ensureRequiredLimiterPatchIsAllocated(
      cellDescriptionsIndex,solverElement,
      solverPatch.getPreviousRefinementStatus());

  // 2. Rollback solution to previous time step
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if ( solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

      _solver->swapSolutionAndPreviousSolution(solverPatch);   // roll back solver

      if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch ) {
        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch); // roll back limiter (must exist!)
      } else {
        const int limiterElement =
            tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
        if (limiterElement!=Solver::NotFound) {
          LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
          copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
          projectDGSolutionOnFVSpace(solverPatch,limiterPatch); // project DG solution on patch
        }
      }
    }
    else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      _solver->swapSolutionAndPreviousSolution(solverPatch);
    }
  }
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
      cellDescriptionsIndex,solverElement,solverPatch.getRefinementStatus());

  // 2. Now roll back to the last valid solution
  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForPassiveFVPatch
  ) { // this is one of the important differences to the global recomputation where we rollback also cells with limiter status == 0
    assertion(solverPatch.getType()==SolverPatch::Type::Cell);
    assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

    if ( solverPatch.getPreviousRefinementStatus() >= _solver->_minimumRefinementStatusForPassiveFVPatch ) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      _limiter->swapSolutionAndPreviousSolution(limiterPatch);
    }
    else { // We need to project limiter data for the previous stage
      _solver->swapSolutionAndPreviousSolution(solverPatch);
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }

    // 2.1 Reset the iterationsToCure on all troubled cells to maximum value if cell is troubled
    if (solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForTroubledCell()) {
      solverPatch.setIterationsToCureTroubledCell(1+_iterationsToCureTroubledCell);
    }
  }

  // 3. Only after the reinitialisation, it is safe to deallocate the limiter patch
  ensureNoLimiterPatchIsAllocatedOnHelperCell(cellDescriptionsIndex,solverElement);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
  double*       solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

  projectOnDGSpace(limiterSolution, solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolutionLocally(
        const int cellDescriptionsIndex, const int solverElement) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForPassiveFVPatch
  ) {
    if (
        solverPatch.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch()
    ) { // these guys are recomputing with the limiter
      const int limiterElement =
          tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
      assertion(limiterElement!=Solver::NotFound);
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(
          solverPatch,cellDescriptionsIndex,limiterElement);

      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,true);
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    }
    else { // these guys are just swapping and projecting
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);

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
  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  assertion1(
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() || limiterElement==Solver::NotFound,
      solverPatch.toString());
  #endif
}

void exahype::solvers::LimitingADERDGSolver::recomputePredictorLocally(
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getRefinementStatus() < _solver->getMinimumRefinementStatusForTroubledCell()  &&
      solverPatch.getPredictorTimeStepSize() > 0
  ) {
    if (
        solverPatch.getRefinementStatus() >= _solver->getMinimumRefinementStatusForActiveFVPatch()
        || // TODO(Dominic): Reassess what is this doing
        (solverPatch.getRefinementStatus() < _solver->getMinimumRefinementStatusForActiveFVPatch() &&
         solverPatch.getPreviousRefinementStatus()>=_solver->getMinimumRefinementStatusForTroubledCell())
    ) {
      _solver->performPredictionAndVolumeIntegral(
          cellDescriptionsIndex,
          element,
          solverPatch.getPredictorTimeStamp(),
          solverPatch.getPredictorTimeStepSize(),
          false/*already uncompressed*/,isAtRemoteBoundary);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::prolongateFaceData(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->prolongateFaceData(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::restriction(
        const int cellDescriptionsIndex,
        const int element) {
  _solver->restriction(cellDescriptionsIndex,element);
}


///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursMetadata(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {
  _solver->mergeNeighboursMetadata(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1. Solve the riemann problems
  mergeNeighboursBasedOnLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,
      pos1,pos2,
      false /* isRecomputation */);

  // 2. Merge the min and max of both cell description's solver's
  // solution value.
  mergeSolutionMinMaxOnFace(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighboursBasedOnLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    const bool                                isRecomputation) const {
  synchroniseTimeStepping(cellDescriptionsIndex1,element1);
  synchroniseTimeStepping(cellDescriptionsIndex2,element2);

  SolverPatch& solverPatch1 = ADERDGSolver::getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = ADERDGSolver::getCellDescription(cellDescriptionsIndex2,element2);
  const int limiterElement1 = tryGetLimiterElement(cellDescriptionsIndex1,solverPatch1.getSolverNumber());
  const int limiterElement2 = tryGetLimiterElement(cellDescriptionsIndex2,solverPatch2.getSolverNumber());

  // We only limit on the finest mesh level
  if (solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel()) {
    assertion2(solverPatch1.getLevel()==getMaximumAdaptiveMeshLevel(),solverPatch1.toString(),solverPatch2.toString());
    // 1. Merge solver solution or limiter solution values in
    // non-overlapping parts of solver and limiter domain:
    if (solverPatch1.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch() &&
        solverPatch2.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch()) {
      if (!isRecomputation) {
        _solver->mergeNeighbours(
            cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
         // Which one is left and right is checked internally again.
      }
    }
    else if (solverPatch1.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch() &&
             solverPatch2.getRefinementStatus()>=_solver->getMinimumRefinementStatusForActiveFVPatch()) {
      _limiter->mergeNeighbours(
          cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2);
          // Which one is left and right is checked internally again.
    }
    // 2. Merge limiter solution values in overlapping part
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
      _limiter->mergeNeighbours(
          cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2);
          // Which one is left and right is checked internally again.
      if (!isRecomputation) {
        _solver->mergeNeighbours(
            cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
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
    if (!isRecomputation) {
      _solver->mergeNeighbours(
          cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
          // Which one is left and right is checked internally again.
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  if ( _solver->getDMPObservables() > 0 ) {
    SolverPatch& solverPatch1 = ADERDGSolver::getCellDescription(cellDescriptionsIndex1,element1);
    SolverPatch& solverPatch2 = ADERDGSolver::getCellDescription(cellDescriptionsIndex2,element2);

    // 1.2. Merge min/max of both solver patches
    const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
    const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
    const int orientation2 = 1-orientation1;

    const int faceIndex1 = 2*direction+orientation1;
    const int faceIndex2 = 2*direction+orientation2;


    if (
        (solverPatch1.getCommunicationStatus()==ADERDGSolver::CellCommunicationStatus &&
        solverPatch1.getFacewiseCommunicationStatus(faceIndex1) >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
        solverPatch1.getFacewiseAugmentationStatus(faceIndex1)  <  ADERDGSolver::MaximumAugmentationStatus) // excludes Ancestors
        ||
        (solverPatch2.getCommunicationStatus()==ADERDGSolver::CellCommunicationStatus &&
        solverPatch2.getFacewiseCommunicationStatus(faceIndex2) >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
        solverPatch2.getFacewiseAugmentationStatus(faceIndex2)  <  ADERDGSolver::MaximumAugmentationStatus) // excludes Ancestors
    ) {
      mergeSolutionMinMaxOnFace(
          solverPatch1,solverPatch2,faceIndex1,faceIndex2);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  const int faceIndex1,
  const int faceIndex2
) const {
  assertion( solverPatch1.getSolverNumber() == solverPatch2.getSolverNumber() );
  const int numberOfObservables = _solver->getDMPObservables();
  double* min1 = DataHeap::getInstance().getData( solverPatch1.getSolutionMin()  ).data() + faceIndex1 * numberOfObservables;
  double* min2 = DataHeap::getInstance().getData( solverPatch2.getSolutionMin()  ).data() + faceIndex2 * numberOfObservables;
  double* max1 = DataHeap::getInstance().getData( solverPatch1.getSolutionMax()  ).data() + faceIndex1 * numberOfObservables;
  double* max2 = DataHeap::getInstance().getData( solverPatch2.getSolutionMax()  ).data() + faceIndex2 * numberOfObservables;

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

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  mergeWithBoundaryDataBasedOnLimiterStatus(
      cellDescriptionsIndex,element,
      solverPatch.getRefinementStatus(),
      posCell,posBoundary,
      false /* isRecomputation */);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const int                                 limiterStatus, //TODO(Dominic): Still necessary?
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      const bool                                isRecomputation) {
  synchroniseTimeStepping(cellDescriptionsIndex,element);

  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

      if (solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForActiveFVPatch()) {
        assertion(limiterStatus<=_solver->_minimumRefinementStatusForPassiveFVPatch || limiterElement!=Solver::NotFound);
        if (!isRecomputation) {
          _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary);
        }
      }
      else {
        assertion(limiterElement!=Solver::NotFound);
        _limiter->mergeWithBoundaryData(cellDescriptionsIndex,limiterElement,posCell,posBoundary);
      }
    }
    else { // solverPatch.getLevel()!=getMaximumAdaptiveMeshLevel()
      if (!isRecomputation) {
        _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary);
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

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMetadata(
      const exahype::MetadataHeap::HeapEntries& metadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int cellDescriptionsIndex,
      const int element) const {
  _solver->mergeWithNeighbourMetadata(metadata,src,dest,cellDescriptionsIndex,element);
}

///////////////////////////////////
// NEIGHBOUR - Time marching
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  sendMinAndMaxToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);

  sendDataToNeighbourBasedOnLimiterStatus(
      toRank,cellDescriptionsIndex,element,src,dest,
      false,/* isRecomputation */
      x,level);

  // send order:   minAndMax,solver,limiter
  // receive order limiter,solver,minAndMax
}

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int numberOfObservables = _solver->getDMPObservables();
  if (numberOfObservables>0) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
    if (ADERDGSolver::holdsFaceData(solverPatch)) {
      const int direction   = tarch::la::equalsReturnIndex(src, dest);
      const int orientation = (1 + dest(direction) - src(direction))/2;
      const int faceIndex   = 2*direction+orientation;

      assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()));
      assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMax()));
      const double* observablesMin = DataHeap::getInstance().getData(
          solverPatch.getSolutionMin()).data() +
          (faceIndex * numberOfObservables);
      const double* observablesMax = DataHeap::getInstance().getData(
          solverPatch.getSolutionMax()).data() +
          (faceIndex * numberOfObservables);

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
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const {
  if ( level==getMaximumAdaptiveMeshLevel() ) {
    SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

    logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " to rank " <<
                 toRank << " at vertex x=" << x << ", level=" << level <<
                 ", source=" << src << ", destination=" << dest <<", limiterStatus="<<solverPatch.getRefinementStatus());

    // solver sends
    if (
        solverPatch.getRefinementStatus()<_solver->getMinimumRefinementStatusForTroubledCell() &&
        isRecomputation==false
    ) {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
    } else {
      _solver->sendEmptyDataToNeighbour(toRank,x,level);
    }

    // limiter sends (receive order must be inverted)
    if (solverPatch.getRefinementStatus()>=_solver->_minimumRefinementStatusForPassiveFVPatch) {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion1(limiterElement!=Solver::NotFound,solverPatch.toString());
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
    } else {
      _limiter->sendEmptyDataToNeighbour(toRank,x,level);
    }

  } else {

    if ( isRecomputation==false ) {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
    } else {
      _solver->sendEmptyDataToNeighbour(toRank,x,level);
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
    assertion1(solverPatch.getRefinementStatus()>=0,solverPatch.toString());

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


    const int direction   = tarch::la::equalsReturnIndex(src, dest);
    const int orientation = (1 + src(direction) - dest(direction))/2;
    const int faceIndex   = 2*direction+orientation;
    if(
        (solverPatch.getCommunicationStatus()                 == ADERDGSolver::CellCommunicationStatus &&
        solverPatch.getFacewiseCommunicationStatus(faceIndex) >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
        solverPatch.getFacewiseAugmentationStatus(faceIndex)  <  ADERDGSolver::MaximumAugmentationStatus)
        ||
        (solverPatch.getFacewiseCommunicationStatus(faceIndex) == ADERDGSolver::CellCommunicationStatus &&
        solverPatch.getCommunicationStatus()                   >= ADERDGSolver::MinimumCommunicationStatusForNeighbourCommunication &&
        solverPatch.getAugmentationStatus()                    <  ADERDGSolver::MaximumAugmentationStatus)
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
          solverPatch,faceIndex,_receivedMin.data(),_receivedMax.data());
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
  const int           faceIndex,
  const double* const min, const double* const max) const {
  assertion1(ADERDGSolver::holdsFaceData(solverPatch),solverPatch.toString());

  double* solutionMin = DataHeap::getInstance().getData( solverPatch.getSolutionMin()  ).data();
  double* solutionMax = DataHeap::getInstance().getData( solverPatch.getSolutionMax()  ).data();

  const int numberOfObservables = _solver->getDMPObservables();
  for (int i=0; i<numberOfObservables; i++) {
    solutionMin[i+faceIndex*numberOfObservables]  = std::min( solutionMin[i+faceIndex*numberOfObservables], min[i] );
    solutionMax[i+faceIndex*numberOfObservables]  = std::max( solutionMax[i+faceIndex*numberOfObservables], max[i] );
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
    const bool initialGrid,
    const int solverNumber) {
  _solver->progressMeshRefinementInPrepareSendToWorker(
      workerRank, fineGridCell, fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell, coarseGridVerticesEnumerator,
      initialGrid, solverNumber);
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

void exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement,
    const bool initialGrid) {
  _solver->progressMeshRefinementInMergeWithWorker(
      localCellDescriptionsIndex,
      receivedCellDescriptionsIndex,receivedElement,
      initialGrid);
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
    const int                                    level) {
  return _solver->progressMeshRefinementInMergeWithMaster(
      worker,localElement,localElement,coarseGridCellDescriptionsIndex,x,level);
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
  const int             cellDescriptionsIndex,
  const int             element,
  const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed,
  const bool            isSkeletonJob):
  _solver(solver),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _neighbourMergePerformed(neighbourMergePerformed),
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
      _cellDescriptionsIndex,_element,
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

exahype::solvers::LimitingADERDGSolver::AdjustLimiterSolutionJob::AdjustLimiterSolutionJob(
  LimitingADERDGSolver& solver,
  SolverPatch&          solverPatch,
  LimiterPatch&         limiterPatch) :
  _solver(solver),
  _solverPatch(solverPatch),
  _limiterPatch(limiterPatch) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::AdjustLimiterSolutionJob::operator()() {
  _solver.adjustLimiterSolution(_solverPatch,_limiterPatch);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
