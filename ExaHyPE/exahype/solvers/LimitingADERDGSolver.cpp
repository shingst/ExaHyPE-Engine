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
          getMinimumLimiterStatusForTroubledCell());
    }
  }
  return result;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLimiterStatusSpreading() {
  bool result = false;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    result |=
        solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        &&
        (static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        !=exahype::solvers::LimiterDomainChange::Regular ||
        solver->getMeshUpdateRequest());
  }
  return result;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalOrGlobalRecomputation() {
  bool result = false;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    result |=
        solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        !=exahype::solvers::LimiterDomainChange::Regular;
  }
  return result;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation(){
  bool result = false;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    result |=
        solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        ==exahype::solvers::LimiterDomainChange::Irregular;
  }
  return result;
}

bool exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation(){
  bool result = false;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    result |=
        solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
        ==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate;
  }
  return result;
}

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
    std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
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
          _iterationsToCureTroubledCell(iterationsToCureTroubledCell),
          _limiterDomainChange(LimiterDomainChange::Regular),
          _nextLimiterDomainChange(LimiterDomainChange::Regular)
{
  assertion(_solver->getNumberOfParameters() == 0);
  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());

  // TODO(WORKAROUND)
  const int numberOfObservables = _solver->getDMPObservables();
  _invalidObservables.resize(numberOfObservables);
  std::fill_n(_invalidObservables.data(),_invalidObservables.size(),-1);
}

void exahype::solvers::LimitingADERDGSolver::updateNextMeshUpdateRequest(const bool& meshUpdateRequest)  {
  _solver->updateNextMeshUpdateRequest(meshUpdateRequest);
}
bool exahype::solvers::LimitingADERDGSolver::getNextMeshUpdateRequest() const  {
  return _solver->getNextMeshUpdateRequest();
}
bool exahype::solvers::LimitingADERDGSolver::getMeshUpdateRequest() const  {
  return _solver->getMeshUpdateRequest();
}
void exahype::solvers::LimitingADERDGSolver::setNextMeshUpdateRequest()  {
  _solver->setNextMeshUpdateRequest();
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
  _coarsestMeshLevel =
      exahype::solvers::Solver::computeMeshLevel(_maximumMeshSize,boundingBoxSize[0]);

  _limiterDomainChange=LimiterDomainChange::Regular;

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
    case exahype::State::AlgorithmSection::LimiterStatusSpreading:
      isMergingMetadata = getLimiterDomainChange()!=LimiterDomainChange::Regular;
      break;
    case exahype::State::AlgorithmSection::MeshRefinement:
      isMergingMetadata = getMeshUpdateRequest();
      assertion( getLimiterDomainChange()!=LimiterDomainChange::IrregularRequiringMeshUpdate || getMeshUpdateRequest());
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

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<static_cast<int>(_limiterDomainChange)<<
           ",nextLimiterDomainChange="<<static_cast<int>(_nextLimiterDomainChange));
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStepFused(
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  _solver->startNewTimeStepFused(isFirstIterationOfBatch,isLastIterationOfBatch);
  ensureLimiterTimeStepDataIsConsistent();

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<static_cast<int>(_limiterDomainChange)<<
           ",nextLimiterDomainChange="<<static_cast<int>(_nextLimiterDomainChange));
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

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::getNextLimiterDomainChange() const {
  return _nextLimiterDomainChange;
}

void exahype::solvers::LimitingADERDGSolver::setNextLimiterDomainChange() {
  _limiterDomainChange     = _nextLimiterDomainChange;
  _nextLimiterDomainChange = LimiterDomainChange::Regular;
}
void exahype::solvers::LimitingADERDGSolver::updateNextLimiterDomainChange(
    exahype::solvers::LimiterDomainChange limiterDomainChange) {
  _nextLimiterDomainChange =
      std::max( _nextLimiterDomainChange, limiterDomainChange );
}
exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::getLimiterDomainChange() const {
  return _limiterDomainChange;
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
void exahype::solvers::LimitingADERDGSolver::markForRefinementBasedOnLimiterStatus(
    SolverPatch& solverPatch) const {
  const int minimumLimiterStatusForRefinement =
      computeMinimumLimiterStatusForRefinement(solverPatch.getLevel());

  if (
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getLimiterStatus() >= minimumLimiterStatusForRefinement &&
      solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
  ) {
    solverPatch.setRefinementRequest(SolverPatch::RefinementRequest::Refine);
  } else if (
      solverPatch.getType()==SolverPatch::Type::Cell &&
      ( // getLimiterDomainChange()!=LimiterDomainChange::Regular || TODO(Dominic): Do we need this?
      solverPatch.getLimiterStatus() >= minimumLimiterStatusForRefinement ||
      solverPatch.getPreviousLimiterStatus() >= minimumLimiterStatusForRefinement)
  ) {
    solverPatch.setRefinementRequest(SolverPatch::RefinementRequest::Keep);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateLimiterStatusDuringLimiterStatusSpreading(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  synchroniseTimeStepping(cellDescriptionsIndex,solverElement);

  if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell()) {
    solverPatch.setLimiterStatus(
        ADERDGSolver::determineLimiterStatus(solverPatch,solverPatch.getNeighbourMergePerformed()));
  } else {
    solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
  }
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
  const int cellDescriptionsIndex = fineGridCell.getCellDescriptionsIndex();
  const int solverElement = _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

    updateLimiterStatusDuringLimiterStatusSpreading(cellDescriptionsIndex,solverElement);

    if (solverPatch.getRefinementRequest()!=SolverPatch::RefinementRequest::Pending) {
      markForRefinementBasedOnLimiterStatus(solverPatch);
    }
  }

  return
      _solver->progressMeshRefinementInEnterCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVerticesEnumerator,
          initialGrid,solverNumber);

}

void exahype::solvers::LimitingADERDGSolver::vetoErasingChildrenRequestBasedOnLimiterStatus(
    SolverPatch& fineGridSolverPatch) const {
  const int minLimiterStatusForRefinement =
      computeMinimumLimiterStatusForRefinement(fineGridSolverPatch.getLevel());
  if (
      (fineGridSolverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::ErasingChildrenRequested ||
      fineGridSolverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::ChangeChildrenToVirtualChildrenRequested)
      &&
      (getLimiterDomainChange()!=LimiterDomainChange::Regular ||
      fineGridSolverPatch.getLimiterStatus()>=minLimiterStatusForRefinement ||
      fineGridSolverPatch.getPreviousLimiterStatus()>=minLimiterStatusForRefinement)
      // TODO(Dominic): Add to docu: This is necessary for not erasing cells in global recomputation and it further adds some laziness in erasing.
      // TODO(Dominic): Add to docu: We always veto based on the previous limiter status too
      // in case we need to go back in time
  ) {
    fineGridSolverPatch.setRefinementEvent(SolverPatch::RefinementEvent::None);
  }
}

bool exahype::solvers::LimitingADERDGSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  const int parentSolverElement =
      _solver->tryGetElement(coarseGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
      solverElement!=Solver::NotFound &&
      parentSolverElement!=Solver::NotFound
  ) {
    SolverPatch& solverPatch = _solver->getCellDescription(
        fineGridCell.getCellDescriptionsIndex(),solverElement);
    _solver->restrictToNextParent(solverPatch,parentSolverElement);
  }

  return
      _solver->progressMeshRefinementInLeaveCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,fineGridPositionOfCell,solverNumber);
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::LimitingADERDGSolver::eraseOrRefineAdjacentVertices(
      const int& cellDescriptionsIndex,
      const int& solverNumber,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) const {
  return _solver->eraseOrRefineAdjacentVertices(
             cellDescriptionsIndex,solverNumber,cellSize);
}

bool exahype::solvers::LimitingADERDGSolver::attainedStableState(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int solverNumber) const {
  return _solver->attainedStableState(
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
    // only for global re-computation
    if ( getLimiterDomainChange()==exahype::solvers::LimiterDomainChange::IrregularRequiringMeshUpdate ) {
      solverPatch.setLimiterStatus(solverPatch.getPreviousLimiterStatus());
    }

    // for global re-computation and mesh refinement
    const bool newLimiterPatchAllocated =
        ensureRequiredLimiterPatchIsAllocated(
            cellDescriptionsIndex,solverElement,solverPatch.getLimiterStatus());
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
int exahype::solvers::LimitingADERDGSolver::computeMinimumLimiterStatusForRefinement(
    int level) const {
  const int levelDelta = getMaximumAdaptiveMeshLevel() - level;

  assertion(levelDelta>=0);
  const int status = _solver->getMinimumLimiterStatusForActiveFVPatch()-1 + (levelDelta-1);
  return (_solver->getMinimumLimiterStatusForTroubledCell()-status>=2) ?
          status : _solver->getMinimumLimiterStatusForTroubledCell()-2;
}

bool exahype::solvers::LimitingADERDGSolver::evaluateLimiterStatusRefinementCriterion(
    const SolverPatch& solverPatch) const {
  return
      solverPatch.getType()==SolverPatch::Type::Cell
      &&
      solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
      &&
      solverPatch.getLimiterStatus() >=
      computeMinimumLimiterStatusForRefinement(solverPatch.getLevel());
}

bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleRefinementCriterion(
    const int cellDescriptionsIndex,
    const int solverElement) const {
//  if (solverElement!=Solver::NotFound) {
//    SolverPatch& solverPatch =
//        ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
//    const int numberOfObservables = _solver->getDMPObservables();
//    if (solverPatch.getType()==SolverPatch::Type::Cell
//        &&
//        numberOfObservables > 0
//    ) {
//      double* observablesMin = DataHeap::getInstance().getData(
//          solverPatch.getSolutionMin()).data();
//      double* observablesMax = DataHeap::getInstance().getData(
//          solverPatch.getSolutionMax()).data();
//
//      bool discreteMaximumPrincipleSatisfied = true;
//      if (
//          solverPatch.getLevel() < getMaximumAdaptiveMeshLevel()
//          &&
//          (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
//         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
//         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
//      ) {
//        const double* localMinPerObservable    = observablesMin+0;
//        const double* localMaxPerObservable    = observablesMin+0;
//        const double* boundaryMinPerObservable = observablesMin+numberOfObservables;
//        const double* boundaryMaxPerObservable = observablesMax+numberOfObservables;
//
//        const double differenceScaling   = _DMPDifferenceScaling;
//        const double relaxationParameter = _DMPMaximumRelaxationParameter * 0.1; // TODO(Dominic): Should be a parameter in the spec file (optional)
//
//        for(int v = 0; v < numberOfObservables; v++) {
//          double boundaryMin = boundaryMinPerObservable[v];
//          double boundaryMax = boundaryMaxPerObservable[v];
//          double scaledDifference = (boundaryMax - boundaryMin) * differenceScaling;
//
//          assertion5(tarch::la::greaterEquals(scaledDifference,0.0),scaledDifference,boundaryMin,boundaryMax,localMinPerObservable[v],localMaxPerObservable[v]);
//          scaledDifference = std::max( scaledDifference, relaxationParameter );
//
//          if((localMinPerObservable[v] < (boundaryMin - scaledDifference)) ||
//              (localMaxPerObservable[v] > (boundaryMax + scaledDifference))) {
//            discreteMaximumPrincipleSatisfied=false;
//          }
//        }
//
//        return false; // TODO(Dominic): Reenable
////        return !discreteMaximumPrincipleSatisfied;
//      }
//
//      // copy local min and max onto boundary
//      //  TODO(Dominic):
//      // 2. Copy the result on the other faces as well
//      for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
//        std::copy_n(
//            observablesMin,numberOfObservables, // past-the-end element
//            observablesMin+i*numberOfObservables);
//        std::copy_n(
//            observablesMax,numberOfObservables, // past-the-end element
//            observablesMax+i*numberOfObservables);
//      }
//    }
//  }
  return false;
}

bool
exahype::solvers::LimitingADERDGSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,const int element) {
  // TODO(Dominic): Reenable if it makes sense
  // refinementRequested |=
  //    evaluateDiscreteMaximumPrincipleRefinementCriterion(cellDescriptionsIndex,element);
  SolverPatch& solverPatch =
      ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  return
      evaluateLimiterStatusRefinementCriterion(solverPatch)
      ||
      (solverPatch.getLimiterStatus()<
      computeMinimumLimiterStatusForRefinement(solverPatch.getLevel())
      &&
      _solver->evaluateRefinementCriterionAfterSolutionUpdate(
          cellDescriptionsIndex,element));
}


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
  _solver->rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
  ensureLimiterPatchTimeStepDataIsConsistent(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStepFused(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  _solver->rollbackToPreviousTimeStepFused(cellDescriptionsIndex,solverElement);
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
    assertion2(solverPatch.getRefinementEvent()==SolverPatch::None,solverPatch.toString(),cellDescriptionsIndex);

    _solver->adjustSolution(solverPatch);

    determineSolverMinAndMax(solverPatch);
    if ( !evaluatePhysicalAdmissibilityCriterion(solverPatch) ) {
      solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
      solverPatch.setRefinementRequest(SolverPatch::RefinementRequest::Keep);
      if ( solverPatch.getLevel()<getMaximumAdaptiveMeshLevel() ) {
        solverPatch.setRefinementRequest(SolverPatch::RefinementRequest::Refine);
      }
    } else {
      _solver->markForRefinement(solverPatch);
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
  assertion(solverPatch.getLimiterStatus()>=0);
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  return getLimiterPatchForSolverPatch(solverPatch,cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::ensureLimiterPatchTimeStepDataIsConsistent(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=Solver::NotFound) {
    assertion1(solverPatch.getPreviousLimiterStatus()>0 || solverPatch.getLimiterStatus()>0,solverPatch.toString());
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
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed,
    const bool  vetoSpawnPredictionJob) {
  // synchroniseTimeStepping(cellDescriptionsIndex,element); // assumes this was done in neighbour merge
  updateSolution(cellDescriptionsIndex,element,isFirstIterationOfBatch);
  UpdateResult result;
  result._limiterDomainChange = updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
                                  cellDescriptionsIndex,element,neighbourMergePerformed);
  // This is important to memorise before calling startNewTimeStepFused
  // TODO(Dominic): Add to docu and/or make cleaner
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  const double memorisedPredictorTimeStamp    = solverPatch.getPredictorTimeStamp();
  const double memorisedPredictorTimeStepSize = solverPatch.getPredictorTimeStepSize();
  result._timeStepSize = startNewTimeStepFused(
      cellDescriptionsIndex,element,isFirstIterationOfBatch,isLastIterationOfBatch);
  result._refinementRequested = evaluateRefinementCriterionAfterSolutionUpdate(cellDescriptionsIndex,element);

  if ( solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell() ) {   // TODO(Dominic): Add to docu. This will spawn or do a compression job right afterwards and must thus come last. This order is more natural anyway
    _solver->performPredictionAndVolumeIntegral( solverPatch,
        memorisedPredictorTimeStamp,memorisedPredictorTimeStepSize,
        false/*already uncompressed*/,vetoSpawnPredictionJob);
  } else { // just perform a restriction of the limiter status to the next parent
    const int parentElement = tryGetElement(
        solverPatch.getParentIndex(),solverPatch.getSolverNumber());
    if (parentElement!=exahype::solvers::Solver::NotFound) {
      _solver->restrictToNextParent(solverPatch,parentElement); // TODO(Dominic): This job should be not stolen
    }
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

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    bool vetoSpawnBackgroundJobs =
        !SpawnPredictionAsBackgroundJob ||
        #if !defined(Parallel) || !defined(SharedMemoryParallelisation)
        isAtRemoteBoundary ||
        #endif
        ADERDGSolver::isInvolvedInProlongationOrRestriction(solverPatch);

    if (
        isFirstIterationOfBatch ||
        isLastIterationOfBatch  ||
        vetoSpawnBackgroundJobs
    ) {
      return fusedTimeStepBody(
          cellDescriptionsIndex,element,
          isFirstIterationOfBatch,isLastIterationOfBatch,
          solverPatch.getNeighbourMergePerformed(),vetoSpawnBackgroundJobs);
    } else {
      int& jobCounter = (isAtRemoteBoundary) ? NumberOfSkeletonJobs: NumberOfEnclaveJobs;
      FusedTimeStepJob fusedTimeStepJob( *this, cellDescriptionsIndex, element,
          solverPatch.getNeighbourMergePerformed(), jobCounter );
      peano::datatraversal::TaskSet spawnedSet( fusedTimeStepJob, peano::datatraversal::TaskSet::TaskType::Background );
      return UpdateResult();
    }
  } else {
    return UpdateResult();
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::LimitingADERDGSolver::update(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary){
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
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
      UpdateResult result;
      updateSolution(cellDescriptionsIndex,element,true);
      result._timeStepSize         = startNewTimeStep(cellDescriptionsIndex,element);
      result._limiterDomainChange  = updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
          cellDescriptionsIndex,element,solverPatch.getNeighbourMergePerformed());  // !!! limiter status must be updated before refinement criterion is evaluated
      result._refinementRequested |= evaluateRefinementCriterionAfterSolutionUpdate(
          cellDescriptionsIndex,element);

      // compress again
      if (CompressionAccuracy>0.0) {
        compress(cellDescriptionsIndex,element,isAtRemoteBoundary);
      }
      return result;
  } else {
    UpdateResult result;
    return result;
  }
}

void exahype::solvers::LimitingADERDGSolver::compress(
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    bool vetoSpawnBackgroundJobs =
      #if !defined(Parallel) || !defined(SharedMemoryParallelisation)
      isAtRemoteBoundary ||
      #endif
      ADERDGSolver::isInvolvedInProlongationOrRestriction(solverPatch);

    _solver->compress(solverPatch,vetoSpawnBackgroundJobs,isAtRemoteBoundary);
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

  // 0. Erase old cells; now it's safe (TODO(Dominic): Add to docu)
  ensureNoLimiterPatchIsAllocatedOnHelperCell(cellDescriptionsIndex,element);
  ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,element);

  // 1. Write back the limiter status to the previous limiter status field
  solverPatch.setPreviousLimiterStatus(solverPatch.getLimiterStatus());

  // 2. Update the solution in the cells
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel()) {
      const int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(solverPatch.getLimiterStatus()>=0);

      if (solverPatch.getLimiterStatus()==0) {
        _solver->updateSolution(solverPatch,backupPreviousSolution);
      }
      else if ( solverPatch.getLimiterStatus()<_solver->_minimumLimiterStatusForActiveFVPatch ) {
        _solver->updateSolution(solverPatch,backupPreviousSolution);

        LimiterPatch& limiterPatch =
            getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _limiter->swapSolutionAndPreviousSolution(limiterPatch);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      }
      else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
        assertion1(limiterElement!=Solver::NotFound,solverPatch.toString());
        _limiter->updateSolution(cellDescriptionsIndex,limiterElement,backupPreviousSolution);

        LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
        _solver->swapSolutionAndPreviousSolution(solverPatch);
        projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
      }
    } else {
      _solver->updateSolution(solverPatch,backupPreviousSolution);
    }
  }
}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int solverElement,
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed
) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  LimiterDomainChange limiterDomainChange = LimiterDomainChange::Regular;
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const bool solutionIsValid =
        evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch) &&
        evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found
    if (
        solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
        solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()
    ) {
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
    } // else: Keep the previously computed min and max values

    limiterDomainChange =
        determineLimiterStatusAfterSolutionUpdate(solverPatch,!solutionIsValid,neighbourMergePerformed);

    ensureNoUnrequiredLimiterPatchIsAllocatedOnComputeCell(cellDescriptionsIndex,solverElement);

    bool limiterPatchAllocated =
        ensureRequiredLimiterPatchIsAllocated(
            cellDescriptionsIndex,solverPatch.getSolverNumber(),solverPatch.getLimiterStatus());
    if ( limiterPatchAllocated ) {
      assertion(solverPatch.getLevel()==getMaximumAdaptiveMeshLevel());
      LimiterPatch& limiterPatch = getLimiterPatchForSolverPatch(cellDescriptionsIndex,solverPatch);

      copyTimeStepDataFromSolverPatch(solverPatch,limiterPatch);
      projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
    }
  } else {
    solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch,neighbourMergePerformed));
    ensureNoLimiterPatchIsAllocatedOnHelperCell(cellDescriptionsIndex,solverElement);
  }
  return limiterDomainChange;
}

exahype::solvers::LimiterDomainChange
exahype::solvers::LimitingADERDGSolver::determineLimiterStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,
    const bool isTroubled,
    const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed) const {
  assertion1(solverPatch.getType()==SolverPatch::Type::Cell,solverPatch.toString());

  LimiterDomainChange limiterDomainChange = LimiterDomainChange::Regular;
  if (isTroubled) {
    solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1);
    if (solverPatch.getLimiterStatus()>_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      // TODO(Dominic): Add to docu: width=2: 5->5 and 4->5 are okay. width=1: 3->3 is okay
      limiterDomainChange = LimiterDomainChange::Regular;
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    } else {
      limiterDomainChange = LimiterDomainChange::Irregular;
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
    }
    if (solverPatch.getLevel()<getMaximumAdaptiveMeshLevel()) {
      limiterDomainChange = LimiterDomainChange::IrregularRequiringMeshUpdate;
    }
  } else {
    if (solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
      solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell());
      solverPatch.setIterationsToCureTroubledCell(
          solverPatch.getIterationsToCureTroubledCell()-1);
      if (solverPatch.getIterationsToCureTroubledCell()==0) {
        solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell()-1);
        solverPatch.setIterationsToCureTroubledCell(_iterationsToCureTroubledCell+1); // TODO(Dominic): Probably not necessary
      }
    } else {
      solverPatch.setLimiterStatus(
          ADERDGSolver::determineLimiterStatus(solverPatch,neighbourMergePerformed));
    }
  }

  assertion2(
      !isTroubled ||
      solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell(),
      isTroubled,
      solverPatch.getLimiterStatus());
  return limiterDomainChange;
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
      assertion(solverPatch.getLimiterStatus()>=0);
      if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        determineSolverMinAndMax(solverPatch);
      } else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch
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
      solverPatch.getLimiterStatus()==0     &&
      solverPatch.getPreviousLimiterStatus()==0
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
      limiterStatus>0
  ) {
    assertion1(solverPatch.getPreviousLimiterStatus()==0,solverPatch.toString());
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

//void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpaceAfterSolutionUpdate(
//    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
//  const double* solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
//  double*       limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
//
//  const double* solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
//    double*       limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
//
//  projectOnFVLimiterSpace(solverSolution, limiterSolution);
//}

// TODO(Dominic): Check that we have rolled back in time as well
void exahype::solvers::LimitingADERDGSolver::rollbackSolutionGlobally(
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  // 1. Ensure limiter patch is allocated (based on previous limiter status
  ensureRequiredLimiterPatchIsAllocated(
      cellDescriptionsIndex,solverElement,
      solverPatch.getPreviousLimiterStatus());

  // 2. Rollback solution to previous time step
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    if ( solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() ) {
      assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

      _solver->swapSolutionAndPreviousSolution(solverPatch);   // roll back solver

      if ( solverPatch.getPreviousLimiterStatus() > 0 ) {
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
    const int cellDescriptionsIndex, const int solverElement) const {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,solverElement);

  // 1. Ensure limiter patch is allocated (based on current limiter status)
  ensureRequiredLimiterPatchIsAllocated(
      cellDescriptionsIndex,solverElement,solverPatch.getLimiterStatus());

  // 2. Now roll back to the last valid solution
  if (
      solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getLimiterStatus()>0
  ) { // this is one of the important differences to the global recomputation where we rollback also cells with limiter status == 0
    assertion(solverPatch.getType()==SolverPatch::Type::Cell);
    assertion(solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None);

    if ( solverPatch.getPreviousLimiterStatus() > 0 ) {
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
    if (solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
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
      solverPatch.getLimiterStatus()>0
  ) {
    if (
        solverPatch.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()
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

      if ( solverPatch.getPreviousLimiterStatus() > 0 ) { // these did just do a swap
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

  if (solverPatch.getLevel()==getMaximumAdaptiveMeshLevel() &&
      solverPatch.getType()==SolverPatch::Type::Cell &&
      solverPatch.getLimiterStatus() < _solver->getMinimumLimiterStatusForTroubledCell()  &&
      solverPatch.getPredictorTimeStepSize() > 0) {
    assertion(solverPatch.getLimiterStatus()>=0);
    if (
        solverPatch.getLimiterStatus() >= _solver->getMinimumLimiterStatusForActiveFVPatch()
        ||
        (solverPatch.getLimiterStatus() < _solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell())
    ) {
      _solver->performPredictionAndVolumeIntegral(
          solverPatch,
          solverPatch.getPredictorTimeStamp(),
          solverPatch.getPredictorTimeStepSize(),
          false/*already uncompressed*/,isAtRemoteBoundary);
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::prolongateAndPrepareRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->prolongateAndPrepareRestriction(cellDescriptionsIndex,element);

  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);
  if (
      solverPatch.getType()==SolverPatch::Type::Ancestor &&
      solverPatch.getPreviousLimiterStatus()>=_solver->getMinimumLimiterStatusForTroubledCell()) {
    solverPatch.setLimiterStatus(_solver->getMinimumLimiterStatusForTroubledCell()-1);
    solverPatch.setFacewiseLimiterStatus(0);
  }
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
    if (solverPatch1.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch() &&
        solverPatch2.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      if (!isRecomputation) {
        _solver->mergeNeighbours(
            cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
         // Which one is left and right is checked internally again.
      }
    }
    else if (solverPatch1.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch() &&
             solverPatch2.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch()) {
      _limiter->mergeNeighbours(
          cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2);
          // Which one is left and right is checked internally again.
    }
    // 2. Merge limiter solution values in overlapping part
    // of solver and limiter domain:
    else if (
        (solverPatch1.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch2.getLimiterStatus() <_solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch2.getLimiterStatus() > 0)
        ||
        (solverPatch1.getLimiterStatus() > 0                                                 &&
         solverPatch1.getLimiterStatus() <_solver->getMinimumLimiterStatusForActiveFVPatch() &&
         solverPatch2.getLimiterStatus()>=_solver->getMinimumLimiterStatusForActiveFVPatch())
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

  SolverPatch& solverPatch1 = ADERDGSolver::getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = ADERDGSolver::getCellDescription(cellDescriptionsIndex2,element2);

  // 1.2. Merge min/max of both solver patches
  const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;

  const int faceIndex1 = 2*direction+orientation1;
  const int faceIndex2 = 2*direction+orientation2;

  mergeSolutionMinMaxOnFace(
      solverPatch1,solverPatch2,faceIndex1,faceIndex2);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& solverPatch1,
  SolverPatch& solverPatch2,
  const int faceIndex1,
  const int faceIndex2
) const {
  if (
      _solver->getDMPObservables() > 0
      &&
      (solverPatch1.getType()==SolverPatch::Cell ||
      solverPatch2.getType()==SolverPatch::Cell)
  ) {
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
  } // else do nothing
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  SolverPatch& solverPatch = ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  mergeWithBoundaryDataBasedOnLimiterStatus(
      cellDescriptionsIndex,element,
      solverPatch.getLimiterStatus(),
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

      if (solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch()) {
        assertion(limiterStatus==0 || limiterElement!=Solver::NotFound);
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
            _invalidObservables, toRank, x, level,
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
                 ", source=" << src << ", destination=" << dest <<", limiterStatus="<<solverPatch.getLimiterStatus());

    // solver sends
    if (
        solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForTroubledCell() &&
        isRecomputation==false
    ) {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
    } else {
      _solver->sendEmptyDataToNeighbour(toRank,x,level);
    }

    // limiter sends (receive order must be inverted)
    if (solverPatch.getLimiterStatus()>0) {
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
          _invalidObservables, toRank, x, level,
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
        ", source=" << src << ", destination=" << dest << ",limiterStatus=" << solverPatch.getLimiterStatus());
    assertion1(solverPatch.getLimiterStatus()>=0,solverPatch.toString());

    if ( solverPatch.getLimiterStatus()<_solver->getMinimumLimiterStatusForActiveFVPatch() ) {
      _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      if (!isRecomputation) {
        _solver->mergeWithNeighbourData(
            fromRank,cellDescriptionsIndex,element,
            src,dest,x,level);
      } else {
        _solver->dropNeighbourData(fromRank,src,dest,x,level);
      }
    } else { // solverPatch.getLimiterStatus()>=ADERDGSolver::MinimumLimiterStatusForActiveFVPatch) {
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

    if (solverPatch.getType()==SolverPatch::Type::Cell) {
      const int direction   = tarch::la::equalsReturnIndex(src, dest);
      const int orientation = (1 + src(direction) - dest(direction))/2;
      const int faceIndex   = 2*direction+orientation;

      const int receivedMaxIndex = DataHeap::getInstance().createData(numberOfObservables, numberOfObservables);
      const int receivedMinIndex = DataHeap::getInstance().createData(numberOfObservables, numberOfObservables);
      DataHeap::HeapEntries& receivedMax = DataHeap::getInstance().getData(receivedMaxIndex);
      DataHeap::HeapEntries& receivedMin = DataHeap::getInstance().getData(receivedMinIndex);
      assertionEquals(DataHeap::getInstance().getData(receivedMaxIndex).size(),static_cast<size_t>(numberOfObservables));
      assertionEquals(DataHeap::getInstance().getData(receivedMinIndex).size(),static_cast<size_t>(numberOfObservables));

      // Inverted send-receive order: TODO(Dominic): Add to docu
      // Send order:    min,max
      // Receive order; max,min
      DataHeap::getInstance().receiveData(
          receivedMax.data(), numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
      DataHeap::getInstance().receiveData(
          receivedMin.data(), numberOfObservables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);

      mergeSolutionMinMaxOnFace(
          solverPatch,faceIndex,receivedMin.data(),receivedMax.data());

      DataHeap::getInstance().deleteData(receivedMinIndex,true);
      DataHeap::getInstance().deleteData(receivedMaxIndex,true);
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
    const int localCellDescriptionsIndex,    const int localElement,
    const int receivedCellDescriptionsIndex, const int receivedElement,
    const bool initialGrid) {
  _solver->progressMeshRefinementInMergeWithWorker(
      localCellDescriptionsIndex,localElement,
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

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerMetadata(
      const MetadataHeap::HeapEntries& receivedMetadata,
      const int                        cellDescriptionsIndex,
      const int                        element) {
  _solver->mergeWithWorkerMetadata(receivedMetadata,cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries messageForMaster =
      _solver->compileMessageForMaster(4);

  // Send additional data to master
  messageForMaster.push_back(
      exahype::solvers::convertToDouble(_limiterDomainChange));

  assertion1(messageForMaster.size()==4,messageForMaster.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending time step data: " <<
        "data[0]=" << messageForMaster[0] <<
        ",data[1]=" << messageForMaster[1] <<
        ",data[2]=" << messageForMaster[2] <<
        ",data[3]=" << messageForMaster[3]);
  }
  DataHeap::getInstance().sendData(
      messageForMaster.data(), messageForMaster.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries messageFromWorker(4); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      messageFromWorker.data(),messageFromWorker.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  // pass message to the ADER-DG solver
  _solver->mergeWithWorkerData(messageFromWorker);

  // merge own flags
  const int firstEntry=3;
  LimiterDomainChange workerLimiterDomainChange =
      exahype::solvers::convertToLimiterDomainChange(messageFromWorker[firstEntry]);
  updateNextLimiterDomainChange(workerLimiterDomainChange); // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Received data from worker:" <<
             " messageFromWorker[0]=" << messageFromWorker[0] <<
             " messageFromWorker[1]=" << messageFromWorker[1] <<
             " messageFromWorker[2]=" << messageFromWorker[2] <<
             " messageFromWorker[3]=" << messageFromWorker[3]);
    logDebug("mergeWithWorkerData(...)","nextLimiterDomainChange=" << static_cast<int>(_nextLimiterDomainChange));
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendDataToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);
  // limiter is only active on the finest mesh level
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToMaster(masterRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const exahype::MetadataHeap::HeapEntries&    workerMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(
      workerRank,workerMetadata,cellDescriptionsIndex,element,x,level);
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->dropWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries messageForWorker = _solver->compileMessageForWorker(8);

  // append additional data
  messageForWorker.push_back(
      exahype::solvers::convertToDouble(_limiterDomainChange));

  assertion1(messageForWorker.size()==8,messageForWorker.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","sending data to worker: " <<
             "data[0]=" << messageForWorker[0] <<
             ",data[1]=" << messageForWorker[1] <<
             ",data[2]=" << messageForWorker[2] <<
             ",data[3]=" << messageForWorker[3] <<
             ",data[4]=" << messageForWorker[4] <<
             ",data[5]=" << messageForWorker[5] <<
             ",data[6]=" << messageForWorker[6] <<
             ",data[7]=" << messageForWorker[7]);
  }
  DataHeap::getInstance().sendData(
      messageForWorker.data(), messageForWorker.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries messageFromMaster(8); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      messageFromMaster.data(),messageFromMaster.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  // pass message to the ADER-DG solver
  _solver->mergeWithMasterData(messageFromMaster);

  // merge own data
  const int firstEntry=7;
  _limiterDomainChange = exahype::solvers::convertToLimiterDomainChange(messageFromMaster[firstEntry]);

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","received data from master:" <<
             " messageFromMaster[0]=" << messageFromMaster[0] <<
             " messageFromMaster[1]=" << messageFromMaster[1] <<
             " messageFromMaster[2]=" << messageFromMaster[2] <<
             " messageFromMaster[3]=" << messageFromMaster[3] <<
             " messageFromMaster[4]=" << messageFromMaster[4] <<
             " messageFromMaster[5]=" << messageFromMaster[5] <<
             " messageFromMaster[6]=" << messageFromMaster[6] <<
             " messageFromMaster[7]=" << messageFromMaster[7]);
    logDebug("mergeWithWorkerData(...)","_limiterDomainChange=" << static_cast<int>(_limiterDomainChange));
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorker(
      workerRank,cellDescriptionsIndex,element,x,level);
  // limiter is only active on the finest mesh level
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToWorker(workerRank,x,level);
  // limiter is only active on the finest mesh level
}

void exahype::solvers::LimitingADERDGSolver::receiveDataFromMaster(
      const int                                    masterRank,
      std::deque<int>&                             receivedHeapDataIndices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const {
  _solver->receiveDataFromMaster(masterRank,receivedHeapDataIndices,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const exahype::MetadataHeap::HeapEntries&     masterMetadata,
    std::deque<int>&                              receivedDataHeapIndices,
    const int                                     cellDescriptionsIndex,
    const int                                     element) const {
  _solver->mergeWithMasterData(
      masterMetadata,receivedDataHeapIndices,cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::dropMasterData(
    std::deque<int>& receivedDataHeapIndices) const {
  _solver->dropMasterData(receivedDataHeapIndices);
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
  int&                  jobCounter):
  _solver(solver),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _neighbourMergePerformed(neighbourMergePerformed),
  _jobCounter(jobCounter) {
  // copy the neighbour merge performed array
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    _jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::LimitingADERDGSolver::FusedTimeStepJob::operator()() {
  _solver.fusedTimeStepBody(_cellDescriptionsIndex,_element,false,false,_neighbourMergePerformed,true);
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    _jobCounter--;
    assertion( _jobCounter>=0 );
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
