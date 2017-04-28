/*
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitingADERDGSolver.h"

#include "kernels/limiter/generic/Limiter.h"

#include "exahype/VertexOperations.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

bool exahype::solvers::LimitingADERDGSolver::irregularChangeOfLimiterDomainOfOneSolver() {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::LimitingADERDGSolver::computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) {
  return _solver->computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);
}

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
    std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
    const double DMPRelaxationParameter,
    const double DMPDifferenceScaling)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNodesPerCoordinateAxis(), solver->getMaximumMeshSize(),
          solver->getMaximumAdaptiveMeshDepth(),
          solver->getTimeStepping()),
          _solver(std::move(solver)),
          _limiter(std::move(limiter)),
          _limiterDomainHasChanged(false),
          _nextLimiterDomainHasChanged(false),
          _DMPMaximumRelaxationParameter(DMPRelaxationParameter),
          _DMPDifferenceScaling(DMPDifferenceScaling)
{
  assertion(_solver->getNumberOfParameters() == 0);

  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());
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

void exahype::solvers::LimitingADERDGSolver::initSolver(const double timeStamp, const tarch::la::Vector<DIMENSIONS,double>& boundingBox) {
  _coarsestMeshLevel = exahype::solvers::Solver::computeMeshLevel(_maximumMeshSize,boundingBox[0]);

  _solver->initSolver(timeStamp, boundingBox);
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->synchroniseTimeStepping(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  _solver->startNewTimeStep();

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min

  setNextGridUpdateRequested();

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<_limiterDomainHasChanged<<",_nextLimiterDomainHasChanged="<<_nextLimiterDomainHasChanged);

  _limiterDomainHasChanged     = _nextLimiterDomainHasChanged;
  _nextLimiterDomainHasChanged = false;
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes() {
  _solver->zeroTimeStepSizes();
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseTimeStepData() {
  _solver->reinitialiseTimeStepData();
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback() {
  _solver->reconstructStandardTimeSteppingDataAfterRollback();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMinCellSize(double minCellSize) {
  _solver->updateNextMinCellSize(minCellSize);
  _nextMinCellSize = _solver->getNextMinCellSize();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMaxCellSize(double maxCellSize) {
  _solver->updateNextMaxCellSize(maxCellSize);
  _nextMaxCellSize = _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMinCellSize() const {
  return _solver->getNextMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMaxCellSize() const {
  return _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinCellSize() const {
  return _solver->getMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMaxCellSize() const {
  return _solver->getMaxCellSize();
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeWithLimiterStatusOfNeighbours(
    SolverPatch& solverPatch,
    const tarch::la::Vector<THREE_POWER_D, int>& neighbourCellDescriptionsIndices) const {
  tarch::la::Vector<DIMENSIONS,int> center(1);
  dfor3(v)
    if (SolverHeap::getInstance().isValidIndex(neighbourCellDescriptionsIndices[vScalar]) &&
        tarch::la::countEqualEntries(v,center)==DIMENSIONS-1) {
      const int neighbourElement =
          _solver->tryGetElement(neighbourCellDescriptionsIndices[vScalar],solverPatch.getSolverNumber());
      if (neighbourElement!=exahype::solvers::Solver::NotFound) {
        SolverPatch& neighbourSolverPatch =
            _solver->getCellDescription(neighbourCellDescriptionsIndices[vScalar],neighbourElement);
        SolverPatch::LimiterStatus neighbourLimiterStatus = neighbourSolverPatch.getLimiterStatus();

        const int direction   = tarch::la::equalsReturnIndex(v,center);
        const int orientation = (1 + v(direction) - center(direction))/2;
        mergeWithLimiterStatus(solverPatch,2*direction+orientation,neighbourLimiterStatus);
      }
    }
  enddforx
}

void exahype::solvers::LimitingADERDGSolver::mergeWithLimiterStatusOfNeighbours(
    const int cellDescriptionsIndex,
    const int solverElement,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&
  indicesAdjacentToFineGridVertices =
      VertexOperations::readCellDescriptionsIndex(
          fineGridVerticesEnumerator,fineGridVertices);

  if (multiscalelinkedcell::adjacencyInformationIsConsistent(
      indicesAdjacentToFineGridVertices)) {
    const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionsIndices =
        multiscalelinkedcell::getIndicesAroundCell(
            indicesAdjacentToFineGridVertices);

    mergeWithLimiterStatusOfNeighbours(solverPatch,neighbourCellDescriptionsIndices);
    solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  }
}

bool exahype::solvers::LimitingADERDGSolver::markForRefinementBasedOnLimiterStatus(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        _solver->getCellDescription(fineGridCell.getCellDescriptionsIndex(),solverElement);
    if (
        solverPatch.getType()==SolverPatch::Type::Cell
        &&
        solverPatch.getLevel() < _coarsestMeshLevel+_maximumAdaptiveMeshDepth
        &&
        (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
         solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
    ) {
      switch (solverPatch.getLimiterStatus()) {
        case SolverPatch::LimiterStatus::Troubled:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
          solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::RefiningRequested);
          return true;
        default:
          return false;
      }
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::markForRefinement(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber)  {
  bool refineFineGridCell = _solver->markForRefinement(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,
      initialGrid,
      solverNumber);

  refineFineGridCell |= markForRefinementBasedOnLimiterStatus(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
      fineGridPositionOfCell,
      initialGrid,
      solverNumber);

  const int solverElement = _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (solverElement!=exahype::solvers::Solver::NotFound) {
    mergeWithLimiterStatusOfNeighbours(
        fineGridCell.getCellDescriptionsIndex(),solverElement,
        fineGridVertices,fineGridVerticesEnumerator);
    allocateOrDeallocateLimiterPatch(
        fineGridCell.getCellDescriptionsIndex(),solverElement);
  }

  return refineFineGridCell;
}

bool exahype::solvers::LimitingADERDGSolver::updateStateInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber)  {
  bool refineFineGridCell =
      _solver->updateStateInEnterCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
          fineGridPositionOfCell,initialGrid,solverNumber);

  return refineFineGridCell;
}

bool exahype::solvers::LimitingADERDGSolver::updateStateInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
//  const int solverElement = _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
// TODO(Dominic): Check that this is done not too often
//  if (solverElement!=exahype::solvers::Solver::NotFound) {
//    mergeLimiterStatusWithAncestors(fineGridCell.getCellDescriptionsIndex(),solverElement);
//  }

  bool eraseFineGridCell =
      _solver->updateStateInLeaveCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridCell,coarseGridVertices,coarseGridVerticesEnumerator,
          fineGridPositionOfCell,solverNumber);

  return eraseFineGridCell;
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

  const int solverElement =
      _solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (solverElement!=exahype::solvers::Solver::NotFound) {
    SolverPatch& solverPatch =
        _solver->getCellDescription(fineGridCell.getCellDescriptionsIndex(),solverElement);

    solverPatch.setPreviousLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  }
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
bool exahype::solvers::LimitingADERDGSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element) {
  return _solver->evaluateRefinementCriterionAfterSolutionUpdate(cellDescriptionsIndex,element);
}


double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double*   tempEigenvalues)  {
  double admissibleTimeStepSize =
      _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::zeroTimeStepSizes(const int cellDescriptionsIndex, const int solverElement) {
  _solver->zeroTimeStepSizes(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback(
      const int cellDescriptionsIndex,
      const int element) const {
  _solver->reconstructStandardTimeSteppingDataAfterRollback(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int solverElement,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  _solver->setInitialConditions(
      cellDescriptionsIndex,solverElement,
      fineGridVertices,fineGridVerticesEnumerator);

  const int limiterElement =
      tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->setInitialConditions(
        cellDescriptionsIndex,limiterElement,
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

/**
 * This method assumes the ADERDG solver's cell-local limiter status has
 * already been determined.
 */
void exahype::solvers::LimitingADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    double** tempStateSizedVectors,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator)  {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  const int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

  // 0. Write back the limiter status to the previous limiter status field
  solverPatch.setPreviousLimiterStatus(solverPatch.getLimiterStatus());

  //   1. Update the solution in the cells
  SolverPatch::LimiterStatus limiterStatus =
      solverPatch.getPreviousLimiterStatus();

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    break;
  case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
  case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);

    assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
    LimiterPatch& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    double* solverSolution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();
    _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);
    double* limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();
    kernels::limiter::generic::c::projectOnFVLimiterSpace(
        solverSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        limiterSolution);
    } break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
  case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
    if (solverPatch.getType()==SolverPatch::Type::Cell) {
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());

      LimiterPatch& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

      limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());

      _limiter->updateSolution(
          cellDescriptionsIndex,limiterElement,
          tempStateSizedVectors,
          tempUnknowns,
          fineGridVertices,fineGridVerticesEnumerator);

      double* solverSolution = DataHeap::getInstance().getData(
          solverPatch.getSolution()).data();
      double* limiterSolution = DataHeap::getInstance().getData(
          limiterPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution);
    }
    } break;
  }

  // 2. Write current limiter status to previous limiter status field
  solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));

  // 3. Allocate or deallocate a limiter patch according to the updated limiter status
  allocateOrDeallocateLimiterPatch(cellDescriptionsIndex,solverPatch.getSolverNumber());
}

//void printSolutionMinOrMax(const double* minOrMax,const int numberOfVariables,const char* identifier) {
//  std::cout << identifier << "=" ;
//
//  for (int i = 0; i < DIMENSIONS_TIMES_TWO*numberOfVariables; ++i) {
//    std::cout << minOrMax[i] << ",";
//  }
//  std::cout << std::endl;
//}
//
//void printNormalFluxes(const double* flux,const int numberOfFaceUnknowns,const char* identifier) {
//  std::cout << identifier << "=" ;
//
//  for (int d=0; d<DIMENSIONS_TIMES_TWO;d++) {
//    for (int i = 0; i < numberOfFaceUnknowns; ++i) {
//      std::cout << flux[i+d*numberOfFaceUnknowns] << ",";
//    }
//    std::cout << "|";
//  }
//  std::cout << std::endl;
//}

bool exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    bool solutionIsValid = evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch) &&
        evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found

    switch (ADERDGSolver::determineLimiterStatus(solverPatch)) {
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
      const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion2(limiterElement!=exahype::solvers::Solver::NotFound,limiterElement,cellDescriptionsIndex);
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
    }
    break;
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      // Already computed.
      break;
    }

    return determineLimiterStatusAfterSolutionUpdate(solverPatch,!solutionIsValid);
  }

  return false;
}

bool exahype::solvers::LimitingADERDGSolver::updateLimiterStatusAndMinAndMaxAfterSetInitialConditions(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    determineSolverMinAndMax(solverPatch);

    bool limiterDomainHasChanged =
        determineLimiterStatusAfterSolutionUpdate(
            solverPatch,
            !evaluatePhysicalAdmissibilityCriterion(solverPatch)); // only evaluate PAD here

    return limiterDomainHasChanged;
  }

  return false;
}

bool exahype::solvers::LimitingADERDGSolver::determineLimiterStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,const bool isTroubled) const {
  bool irregularLimiterDomainChange=false;

  SolverPatch::LimiterStatus limiterStatus =
      exahype::solvers::ADERDGSolver::determineLimiterStatus(solverPatch);
  solverPatch.setLimiterStatus(limiterStatus);

  if (isTroubled) {
    switch (limiterStatus) {
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      irregularLimiterDomainChange=true;
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::Troubled);
      ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      // do nothing
      break;
    }
  } else {
    switch (solverPatch.getPreviousLimiterStatus()) {
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      // do nothing
      break;
    case SolverPatch::LimiterStatus::Troubled:
      solverPatch.setLimiterStatus(SolverPatch::LimiterStatus::NeighbourOfTroubled1);
      ADERDGSolver::writeLimiterStatusOnBoundary(solverPatch);
      break;
    }
  }

  return irregularLimiterDomainChange;
}


bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // 1. Check if the DMP is satisfied and search for the min and max
  // Write the new min and max to the storage reserved for face 0
  bool dmpIsSatisfied = kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch(
        solution,_solver->getNumberOfVariables(),_solver->getNodesPerCoordinateAxis(),
        _DMPMaximumRelaxationParameter, _DMPDifferenceScaling,
        solutionMin,solutionMax);

  // 2. Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy_n(
        solutionMin,_numberOfVariables,
        solutionMin+i*_numberOfVariables);
    std::copy_n(
        solutionMax,_numberOfVariables,
        solutionMax+i*_numberOfVariables);
  }

  return dmpIsSatisfied;
}

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch) {
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  return _solver->isPhysicallyAdmissible(
      solutionMin,solutionMax,
      solverPatch.getOffset()+0.5*solverPatch.getSize(),solverPatch.getSize(),
      solverPatch.getCorrectorTimeStamp(),solverPatch.getCorrectorTimeStepSize());
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

    switch(solverPatch.getPreviousLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      determineSolverMinAndMax(solverPatch);
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
      break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();

  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // Write the result to the face with index "0"
  kernels::limiter::generic::c::findCellLocalMinAndMax(
      solution,_numberOfVariables,_nodesPerCoordinateAxis,solutionMin,solutionMax);

  // Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy_n(
        solutionMin,_numberOfVariables,
        solutionMin+i*_numberOfVariables);
    std::copy_n(
        solutionMax,_numberOfVariables,
        solutionMax+i*_numberOfVariables);
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  double* limiterSolution = DataHeap::getInstance().getData(
      limiterPatch.getSolution()).data();

  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // Write the result to the face with index "0"
  kernels::limiter::generic::c::findCellLocalLimiterMinAndMax(
      limiterSolution,_numberOfVariables,_nodesPerCoordinateAxis,_limiter->getGhostLayerWidth(),solutionMin,solutionMax);

  // Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy_n(
        solutionMin,_numberOfVariables, // past-the-end element
        solutionMin+i*_numberOfVariables);
    std::copy_n(
        solutionMax,_numberOfVariables, // past-the-end element
        solutionMax+i*_numberOfVariables);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithLimiterStatus(
    SolverPatch& solverPatch,
    const int direction,
    const SolverPatch::LimiterStatus& otherLimiterStatus) const {
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();

  const int limiterStatusAsInt      = static_cast<int>(limiterStatus);
  const int otherLimiterStatusAsInt = static_cast<int>(otherLimiterStatus);

  const int newLimiterStatusAsInt =
          ( limiterStatus==SolverPatch::Troubled )
          ?
          limiterStatusAsInt
          :
          std::max (
                0,
                std::max( limiterStatusAsInt, otherLimiterStatusAsInt ) - 1
          );
  solverPatch.setFacewiseLimiterStatus( direction,static_cast<SolverPatch::LimiterStatus>(newLimiterStatusAsInt) );

  // On coarser grid level, we use only four limiter statuses (T,NT1,NT2,Ok).
  if (solverPatch.getLevel() < _coarsestMeshLevel+_maximumAdaptiveMeshDepth
      &&
      newLimiterStatusAsInt <
      static_cast<int>(SolverPatch::LimiterStatus::NeighbourOfTroubled2)
  ) {
    solverPatch.setFacewiseLimiterStatus(direction,SolverPatch::LimiterStatus::Ok);
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourLimiterStatus(
    SolverPatch& solverPatch,
    const int faceIndex,
    const SolverPatch::LimiterStatus& neighbourLimiterStatus,
    const SolverPatch::LimiterStatus& neighbourOfNeighbourLimiterStatus) const {
  SolverPatch::LimiterStatus limiterStatus = ADERDGSolver::determineLimiterStatus(solverPatch);

  mergeWithLimiterStatus(solverPatch,faceIndex,neighbourLimiterStatus);

  if (limiterStatus==SolverPatch::LimiterStatus::Ok) {
    if (neighbourLimiterStatus==SolverPatch::LimiterStatus::Ok
        && neighbourOfNeighbourLimiterStatus==SolverPatch::LimiterStatus::Troubled) {
      solverPatch.setFacewiseLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourOfTroubled1);
    }
  }
}

bool exahype::solvers::LimitingADERDGSolver::allocateOrDeallocateLimiterPatch(
        const int cellDescriptionsIndex,
        const int solverElement) const {
  assertion(solverElement!=exahype::solvers::Solver::NotFound);
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  const SolverPatch::LimiterStatus limiterStatus = ADERDGSolver::determineLimiterStatus(solverPatch);
  int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,solverElement);

  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    switch(solverPatch.getType()) {
    case SolverPatch::Type::EmptyAncestor:
    case SolverPatch::Type::Ancestor:
    case SolverPatch::Type::Descendant:
    case SolverPatch::Type::EmptyDescendant:
    case SolverPatch::Type::Erased: {
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      limiterPatch.setType(LimiterPatch::Type::Erased);
      _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch);

      tarch::multicore::Lock lock(_heapSemaphore);
      LimiterHeap::getInstance().getData(cellDescriptionsIndex).erase(
          LimiterHeap::getInstance().getData(cellDescriptionsIndex).begin()+limiterElement);
      lock.free();
    } break;
    case SolverPatch::Type::Cell:
      switch(limiterStatus) {
      case SolverPatch::LimiterStatus::Ok: {
        LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        limiterPatch.setType(LimiterPatch::Type::Erased);
        _limiter->ensureNoUnnecessaryMemoryIsAllocated(limiterPatch);

        tarch::multicore::Lock lock(_heapSemaphore);
        LimiterHeap::getInstance().getData(cellDescriptionsIndex).erase(
            LimiterHeap::getInstance().getData(cellDescriptionsIndex).begin()+limiterElement);
        lock.free();
      } break;
      default:
        // do nothing
        break;
      }
      break;
    }

  } else {
    if (solverPatch.getType()==SolverPatch::Type::Cell) {
      switch (limiterStatus) {
      case SolverPatch::LimiterStatus::Troubled:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
        tarch::multicore::Lock lock(_heapSemaphore);
        exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
            cellDescriptionsIndex,
            solverPatch.getSolverNumber(),
            LimiterPatch::Type::Cell,
            LimiterPatch::RefinementEvent::None,
            solverPatch.getLevel(),
            solverPatch.getParentIndex(),
            solverPatch.getSize(),
            solverPatch.getOffset());
        lock.free();

        assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolution()),solverPatch.toString());
        assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),solverPatch.toString());

        limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        _limiter->ensureNecessaryMemoryIsAllocated(limiterPatch);

        // Copy more geometry information from the solver patch
        limiterPatch.setIsInside(solverPatch.getIsInside());
      } return true;
      case SolverPatch::LimiterStatus::Ok:
        // do nothing
        break;
      }
    }
  }

  return false;
}

void exahype::solvers::LimitingADERDGSolver::projectDGSolutionOnFVSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
  double*       limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();

  // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
  kernels::limiter::generic::c::projectOnFVLimiterSpace(
      solverSolution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),
      _limiter->getGhostLayerWidth(),
      limiterSolution);
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers(
    const int cellDescriptionsIndex,
    const int solverElement,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);

  // 0. Update the limiter status
  solverPatch.setLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
  // 1. Ensure memory is allocated
  allocateOrDeallocateLimiterPatch(cellDescriptionsIndex,solverElement);

  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    int limiterElement = _limiter->tryGetElement(
        fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());

    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      // TODO(Dominic): Add to docu: We need to rollback the solution for this case since we need to supply the // NeighbourOfTroubled1 neighbours with solution values from the old time step.
      switch (solverPatch.getPreviousLimiterStatus()) {
      case SolverPatch::LimiterStatus::Troubled:  // TODO(Dominic): Add to docu: Here we work with a valid old FVM solution
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _limiter->rollbackSolution(
            fineGridCell.getCellDescriptionsIndex(),limiterElement,
            fineGridVertices,fineGridVerticesEnumerator);
      } break;
      case SolverPatch::LimiterStatus::Ok: // TODO(Dominic): Add to docu: Here we work with a valid old ADER-DG solution
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4: { // TODO(Dominic): Possibly dangerous if nans appear
        _solver->rollbackSolution(
            fineGridCell.getCellDescriptionsIndex(),solverElement,
            fineGridVertices,fineGridVerticesEnumerator);

        if (solverPatch.getPreviousLimiterStatus()==SolverPatch::LimiterStatus::Ok) {
          const int limiterElement = _limiter->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());
          assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
          assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getPreviousSolution()),solverPatch.toString());
          assertion1(DataHeap::getInstance().isValidIndex(solverPatch.getSolution()),solverPatch.toString());

          LimiterPatch& limiterPatch = _limiter->getCellDescription(fineGridCell.getCellDescriptionsIndex(),limiterElement);

          projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
        } else {
          _limiter->rollbackSolution(cellDescriptionsIndex,limiterElement,fineGridVertices,fineGridVerticesEnumerator);
        }
      } break;
      }
    } break;
    case SolverPatch::LimiterStatus::Ok: {
      // do nothing
    } break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::projectFVSolutionOnDGSpace(
    SolverPatch& solverPatch,LimiterPatch& limiterPatch) const {
  const double* limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
  double*       solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

  // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
  kernels::limiter::generic::c::projectOnDGSpace(
      limiterSolution,_solver->getNumberOfVariables(),
      _solver->getNodesPerCoordinateAxis(),
      _limiter->getGhostLayerWidth(),
      solverSolution);
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
        const int cellDescriptionsIndex, const int element,
        exahype::solvers::SolutionUpdateTemporaryVariables& solutionUpdateTemporaryVariables,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  if (solverPatch.getType()==SolverPatch::Type::Cell) {
    const int limiterElement =
        _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2: {
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      // Copy time step data from the solver patch
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      // 1. Evolve solution to desired  time step again
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
          solutionUpdateTemporaryVariables._tempStateSizedVectors[solverPatch.getSolverNumber()],
          solutionUpdateTemporaryVariables._tempUnknowns[solverPatch.getSolverNumber()],
          fineGridVertices,fineGridVerticesEnumerator);
      // 2. Project FV solution on ADER-DG space
      projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
    } break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: { // TODO(Dominic): Add to docu: Here, were just went back one time step to supply the NT neighbours with old limiter unknowns.
      switch (solverPatch.getPreviousLimiterStatus()) {
      case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
      case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      case SolverPatch::LimiterStatus::Ok: {
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _solver->swapSolutionAndPreviousSolution(cellDescriptionsIndex,element);
        LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        projectDGSolutionOnFVSpace(solverPatch,limiterPatch);
      } break;
      case SolverPatch::LimiterStatus::Troubled: // TODO(Dominic): Update docu
      case SolverPatch::LimiterStatus::NeighbourOfTroubled1:   // TODO(Dominic): Update docu
      case SolverPatch::LimiterStatus::NeighbourOfTroubled2: { // TODO(Dominic): Update docu
        assertion(limiterElement!=exahype::solvers::Solver::NotFound);
        _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);
        LimiterPatch& limiterPatch    = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        projectFVSolutionOnDGSpace(solverPatch,limiterPatch);
      } break;
      }
    } break;
    case SolverPatch::LimiterStatus::Ok: {
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
      // do nothing
    }  break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::recomputePredictor(
    const int cellDescriptionsIndex,
    const int element,
    exahype::solvers::PredictionTemporaryVariables& predictionTemporaryVariables,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getPredictorTimeStepSize() > 0) {
    switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      _solver->performPredictionAndVolumeIntegral(
          solverPatch,
          predictionTemporaryVariables._tempSpaceTimeUnknowns    [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempSpaceTimeFluxUnknowns[solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempUnknowns             [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempFluxUnknowns         [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempStateSizedVectors    [solverPatch.getSolverNumber()],
          predictionTemporaryVariables._tempPointForceSources    [solverPatch.getSolverNumber()]);
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      switch (solverPatch.getPreviousLimiterStatus()) {
      case SolverPatch::LimiterStatus::Troubled:
        _solver->performPredictionAndVolumeIntegral(
            solverPatch,
            predictionTemporaryVariables._tempSpaceTimeUnknowns    [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempSpaceTimeFluxUnknowns[solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempUnknowns             [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempFluxUnknowns         [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempStateSizedVectors    [solverPatch.getSolverNumber()],
            predictionTemporaryVariables._tempPointForceSources    [solverPatch.getSolverNumber()]);
        break;
      default:
        break;
      }
      break;
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::prepareNextNeighbourMerging(
    const int cellDescriptionsIndex,const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  _solver->prepareNextNeighbourMerging(
      cellDescriptionsIndex,element,
      fineGridVertices,fineGridVerticesEnumerator);

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->prepareNextNeighbourMerging(
        cellDescriptionsIndex,limiterElement,
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateLimiterStatus(
    const int cellDescriptionsIndex,const int element) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = ADERDGSolver::determineLimiterStatus(solverPatch);

  for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
    solverPatch.setFacewiseLimiterStatus(i,limiterStatus);
  }
}

void exahype::solvers::LimitingADERDGSolver::updatePreviousLimiterStatus(
    const int cellDescriptionsIndex,const int element) const  {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  solverPatch.setPreviousLimiterStatus(ADERDGSolver::determineLimiterStatus(solverPatch));
}

void exahype::solvers::LimitingADERDGSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->preProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  _limiter->preProcess(cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->postProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  _limiter->postProcess(cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  switch(ADERDGSolver::determineLimiterStatus(solverPatch)) {
    case SolverPatch::LimiterStatus::Ok:
      _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      // do nothing;
      break;
    case SolverPatch::LimiterStatus::Troubled:
      if (solverPatch.getType()==SolverPatch::Type::Ancestor) {
        _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);
      }
      break;
  }
}

/**
 * Should be called in leaveCell(...)
 */
void exahype::solvers::LimitingADERDGSolver::mergeLimiterStatusWithAncestors(
    const int cellDescriptionsIndex, const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  if (solverPatch.getType()==SolverPatch::Type::Cell
      || solverPatch.getType()==SolverPatch::Type::Ancestor) {
    int parentCellDescriptionsIndex = solverPatch.getParentIndex();
    int parentElement =
        _solver->tryGetElement(parentCellDescriptionsIndex,solverPatch.getSolverNumber());

    while (SolverHeap::getInstance().isValidIndex(parentCellDescriptionsIndex)) {
      SolverPatch& parentPatch = _solver->getCellDescription(parentCellDescriptionsIndex,parentElement);

      tarch::la::Vector<DIMENSIONS,int> subcellIndex =
          exahype::amr::computeSubcellIndex(solverPatch.getOffset(),solverPatch.getSize(),parentPatch.getOffset());

      const int levelDelta = solverPatch.getLevel() - parentPatch.getLevel();
      for (int d = 0; d < DIMENSIONS; d++) {
        // only merge if the cell's face is part of one of the parent cell's faces
        if (subcellIndex[d]==0 ||
            subcellIndex[d]==tarch::la::aPowI(levelDelta,3)-1) {
          const int faceIndex = 2*d + ((subcellIndex[d]==0) ? 0 : 1); // Do not remove brackets.
          tarch::multicore::Lock lock(_heapSemaphore); // TODO(Dominic): Use another semaphore?
          mergeWithLimiterStatus(parentPatch,faceIndex,solverPatch.getFacewiseLimiterStatus(faceIndex));
        }
      }

      // next parent
      if (parentPatch.getType()==SolverPatch::Type::EmptyAncestor) {
        parentCellDescriptionsIndex = solverPatch.getParentIndex();
        parentElement =
            _solver->tryGetElement(parentCellDescriptionsIndex,solverPatch.getSolverNumber());
      } else if (parentPatch.getType()==SolverPatch::Type::Ancestor) {
        parentCellDescriptionsIndex = -1; // this guy has to continue as soon as he is touched in leaveCell(...)
        parentElement = -1;
      }
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::restrictData(
    const int cellDescriptionsIndex,
    const int element,
    const int parentCellDescriptionsIndex,
    const int parentElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {

  _solver->restrictData(cellDescriptionsIndex,element,parentCellDescriptionsIndex,parentElement,subcellIndex);
}


///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& pLeft,
  SolverPatch& pRight,
  const int faceIndexLeft,
  const int faceIndexRight
) const {
  if (pLeft.getType()==SolverPatch::Cell ||
      pRight.getType()==SolverPatch::Cell) {
    assertion( pLeft.getSolverNumber() == pRight.getSolverNumber() );
    const int numberOfVariables = getNumberOfVariables();
    double* minLeft  = DataHeap::getInstance().getData( pLeft.getSolutionMin()  ).data()  + faceIndexLeft  * numberOfVariables;
    double* minRight = DataHeap::getInstance().getData( pRight.getSolutionMin()  ).data() + faceIndexRight * numberOfVariables;
    double* maxLeft  = DataHeap::getInstance().getData( pLeft.getSolutionMax()  ).data()  + faceIndexLeft  * numberOfVariables;
    double* maxRight = DataHeap::getInstance().getData( pRight.getSolutionMax()  ).data() + faceIndexRight * numberOfVariables;

    for (int i=0; i<numberOfVariables; i++) {
      const double min = std::min(
          *(minLeft+i),
          *(minRight+i)
      );
      const double max = std::max(
          *(maxLeft+i),
          *(maxRight+i)
      );

      *(minLeft+i)  = min;
      *(minRight+i) = min;

      *(maxLeft+i)  = max;
      *(maxRight+i) = max;
    }
  } // else do nothing
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1. Communicate the data between the solvers that is necessary for the solves
  mergeNeighboursBasedOnLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,
      pos1,pos2,
      false, /* isRecomputation */
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);

  // 2. Merge the min and max of both cell description's solver's
  // solution value.
  // 2.1 Determine "left" and "right" solverPatch
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(cellDescriptionsIndexRight,elementRight);

  // 2.2. Merge min/max of both solver patches
  mergeSolutionMinMaxOnFace(solverPatchLeft,solverPatchRight,faceIndexLeft,faceIndexRight);
}

// TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
// They depend on the isRecomputation value
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursBasedOnLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    const bool                                isRecomputation,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {

  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);
  const int limiterElement1 = tryGetLimiterElement(cellDescriptionsIndex1,solverPatch1.getSolverNumber());
  const int limiterElement2 = tryGetLimiterElement(cellDescriptionsIndex2,solverPatch2.getSolverNumber());

  // 0.a Merge the limiter status of the neighours
  const int direction    = tarch::la::equalsReturnIndex(pos1,pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;
  SolverPatch::LimiterStatus limiterStatus1 = solverPatch1.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus2 = solverPatch2.getLimiterStatus();
  mergeWithLimiterStatus(solverPatch1,2*direction+orientation1,limiterStatus2);
  mergeWithLimiterStatus(solverPatch2,2*direction+orientation2,limiterStatus1);

  // 1. Merge solver solution or limiter solution values in
  // non-overlapping parts of solver and limiter domain:
  switch (limiterStatus1) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::Ok:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
          if (!isRecomputation) {
            _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::Troubled:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
          if (solverPatch1.getType()==SolverPatch::Type::Cell &&
              solverPatch2.getType()==SolverPatch::Type::Cell) {
            _limiter->mergeNeighbours(
                cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
  }
  // 2. Merge limiter solution values in overlapping part
  // of solver and limiter domain:
  switch (limiterStatus1) {
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
          assertion2(limiterElement1!=exahype::solvers::Solver::NotFound,solverPatch1.toString(),solverPatch2.toString());
          assertion2(limiterElement2!=exahype::solvers::Solver::NotFound,solverPatch2.toString(),solverPatch1.toString());

          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          if (!isRecomputation) {
            _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
        case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
          assertion2(limiterElement1!=exahype::solvers::Solver::NotFound,solverPatch1.toString(),solverPatch2.toString());
          assertion2(limiterElement2!=exahype::solvers::Solver::NotFound,solverPatch2.toString(),solverPatch1.toString());

          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          if (!isRecomputation) {
            _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMax(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1.1. Determine "left" and "right" solverPatch
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(cellDescriptionsIndexRight,elementRight);

  // 1.2. Merge min/max of both solver patches
  mergeSolutionMinMaxOnFace(solverPatchLeft,solverPatchRight,faceIndexLeft,faceIndexRight);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  mergeWithBoundaryDataBasedOnLimiterStatus(
      cellDescriptionsIndex,element,
      solverPatch.getLimiterStatus(),
      posCell,posBoundary,
      false, /* isRecomputation */
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const SolverPatch::LimiterStatus&         limiterStatus,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      const bool                                isRecomputation,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4:
    case SolverPatch::LimiterStatus::Ok:
      assertion(limiterStatus==SolverPatch::LimiterStatus::Ok || limiterElement!=exahype::solvers::Solver::NotFound);
      if (!isRecomputation) {
        _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
                                       tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      }
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->mergeWithBoundaryData(cellDescriptionsIndex,limiterElement,posCell,posBoundary,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      break;
    default:
      break;
  }
}

#ifdef Parallel
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerNeighbourCommunication    = 1;
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 0;
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerMasterWorkerCommunication = 0;

///////////////////////////////////
// NEIGHBOUR - Mesh refinement
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMetadata(
      const exahype::MetadataHeap::HeapEntries& neighbourMetadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int cellDescriptionsIndex,
      const int element) {
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithNeighbourMetadata(neighbourMetadata,src,dest,cellDescriptionsIndex,limiterElement);
  } // There is no drop method for metadata necessary
  _solver->mergeWithNeighbourMetadata(neighbourMetadata,src,dest,cellDescriptionsIndex,element);

  // Refine according to the merged limiter status
  const SolverPatch::LimiterStatus neighbourLimiterStatus =
     static_cast<SolverPatch::LimiterStatus>(neighbourMetadata[exahype::MetadataLimiterStatus].getU());
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  // TODO(Dominic): Update
  if (
      neighbourLimiterStatus==SolverPatch::LimiterStatus::Troubled
      &&
      solverPatch.getLevel()<_coarsestMeshLevel+_maximumAdaptiveMeshDepth
      &&
      solverPatch.getType()==SolverPatch::Type::Cell
      &&
      (solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::None ||
      solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::DeaugmentingChildrenRequested ||
      solverPatch.getRefinementEvent()==SolverPatch::RefinementEvent::AugmentingRequested)
  ) {
    solverPatch.setRefinementEvent(SolverPatch::RefinementEvent::RefiningRequested);
  }

  // TODO(Dominic)
  if (tarch::la::equals(src,dest)==DIMENSIONS-1) {
    const int direction = peano::utils::dLinearised(src-dest+1,3);
    mergeWithLimiterStatus(solverPatch,direction,neighbourLimiterStatus);

    // old code keep for reference in case we go back to face index
    // const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
    // assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
    // const int faceIndex = 2 * normalOfExchangedFace +
    //        (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
    //                                                                            // |src|dest| : 1; |dest|src| : 0
  }
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
}

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
                                                                          // |src|dest| : 1; |dest|src| : 0

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()));
  assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMax()));

  // We append all the max values to the min values.
  // And then append the limiter status as double
  std::vector<double> minAndMaxToSend(2*_numberOfVariables);
  for (int i=0; i<_numberOfVariables; i++) {
    minAndMaxToSend[i]                    = DataHeap::getInstance().getData( solverPatch.getSolutionMin() )[faceIndex*_numberOfVariables+i];
    minAndMaxToSend[i+_numberOfVariables] = DataHeap::getInstance().getData( solverPatch.getSolutionMax() )[faceIndex*_numberOfVariables+i];
  }

  DataHeap::getInstance().sendData(
      minAndMaxToSend, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
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
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();

  logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " from rank " <<
               toRank << " at vertex x=" << x << ", level=" << level <<
               ", src=" << src << ", dest=" << dest <<", limiterStatus="<<limiterStatus);

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok: {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      } break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
    case SolverPatch::LimiterStatus::Troubled: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:{
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

  _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  mergeWithNeighbourDataBasedOnLimiterStatus(
      fromRank,neighbourMetadata,cellDescriptionsIndex,element,src,dest,
      false,/*isRecomputation*/
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);

  mergeWithNeighbourMinAndMax(fromRank,cellDescriptionsIndex,element,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataBasedOnLimiterStatus(
    const int                                    fromRank,
    const exahype::MetadataHeap::HeapEntries&    neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const bool                                   isRecomputation,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();

  logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
          fromRank << " at vertex x=" << x << ", level=" << level <<
          ", src=" << src << ", dest=" << dest << ",limiterStatus=" << SolverPatch::toString(limiterStatus));

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled3:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled4: {
      _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      if (!isRecomputation) {
        _solver->mergeWithNeighbourData(
                 fromRank,neighbourMetadata,cellDescriptionsIndex,element,
                 src,dest,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);
      } else {
        _solver->dropNeighbourData(fromRank,src,dest,x,level);
      }
      }break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled1:
    case SolverPatch::LimiterStatus::NeighbourOfTroubled2:{
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->mergeWithNeighbourData(
          fromRank,neighbourMetadata,cellDescriptionsIndex,limiterElement,
          src,dest,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);
      _solver->dropNeighbourData(fromRank,src,dest,x,level);
      } break;
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
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  const int receivedMinMaxIndex = DataHeap::getInstance().createData(2*_numberOfVariables, 2*_numberOfVariables);
  assertion(DataHeap::getInstance().getData(receivedMinMaxIndex).size()==static_cast<unsigned int>(2*_numberOfVariables));
  double* receivedMinAndMax = DataHeap::getInstance().getData(receivedMinMaxIndex).data();

  DataHeap::getInstance().receiveData(receivedMinAndMax, 2*_numberOfVariables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
  mergeSolutionMinMaxOnFace(solverPatch,faceIndex,receivedMinAndMax,receivedMinAndMax+_numberOfVariables);

  DataHeap::getInstance().deleteData(receivedMinMaxIndex,true);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  SolverPatch,
  const int           faceIndex,
  const double* const min, const double* const max) const {
  if (SolverPatch.getType() == SolverPatch::Cell ||
      SolverPatch.getType() == SolverPatch::Ancestor ||
      SolverPatch.getType() == SolverPatch::Descendant
      ) {
    double* solutionMin = DataHeap::getInstance().getData( SolverPatch.getSolutionMin()  ).data();
    double* solutionMax = DataHeap::getInstance().getData( SolverPatch.getSolutionMax()  ).data();

    for (int i=0; i<_numberOfVariables; i++) {
      solutionMin[i+faceIndex*_numberOfVariables]  = std::min( solutionMin[i+faceIndex*_numberOfVariables], min[i] );
      solutionMax[i+faceIndex*_numberOfVariables]  = std::max( solutionMax[i+faceIndex*_numberOfVariables], max[i] );
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
      DataHeap::getInstance().receiveData(
          fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

///////////////////////////////////////
// NEIGHBOUR - Limiter status spreading
///////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendLimiterStatusToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
                                                                          // |src|dest| : 1; |dest|src| : 0
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const double limiterStatusAsDouble =
      static_cast<double>(solverPatch.getFacewiseLimiterStatus(faceIndex));

  DataHeap::getInstance().sendData(
      &limiterStatusAsDouble,1,toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataInsteadOfLimiterStatusToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  std::vector<double> emptyMessage(0);
  DataHeap::getInstance().sendData(
      emptyMessage, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourLimiterStatus(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  // Merge the limiter status.
  const int receivedLimiterStatusIndex = DataHeap::getInstance().createData(0, 1);
  assertion(DataHeap::getInstance().getData(receivedLimiterStatusIndex).empty());
  DataHeap::getInstance().receiveData(receivedLimiterStatusIndex,fromRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
  const double* const receivedLimiterStatus =
      DataHeap::getInstance().getData(receivedLimiterStatusIndex).data();
  assertion1(std::lround(*receivedLimiterStatus) <= 3,*receivedLimiterStatus);
  assertion1(std::lround(*receivedLimiterStatus) >= 0,*receivedLimiterStatus);
  SolverPatch::LimiterStatus neighbourLimiterStatus =
      static_cast<SolverPatch::LimiterStatus>(std::lround(*receivedLimiterStatus));
  mergeWithLimiterStatus(solverPatch,faceIndex,neighbourLimiterStatus);

  // Clean up.
  DataHeap::getInstance().deleteData(receivedLimiterStatusIndex,true); // TODO(Dominic): We are dealing here with an array of size 1 here.
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourLimiterStatus(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::getInstance().receiveData(fromRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
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
  _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
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

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

///////////////////////////////////////
// FORK OR JOIN
///////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorkerOrMasterDueToForkOrJoin(
      toRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorkerOrMasterDueToForkOrJoin(
          toRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
          toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
  _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
      fromRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerOrMasterDataDueToForkOrJoin(
        fromRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
        fromRank,x,level);
  } // !!! Receive order must be the same in master<->worker comm.
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
  _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->sendDataToMaster(masterRank,x,level);
//  _limiter->sendDataToMaster(masterRank,x,level);

  // Send the information to master if limiter status has changed or not
  std::vector<double> dataToSend(0,1);
  dataToSend.push_back(_limiterDomainHasChanged ? 1.0 : -1.0); // TODO(Dominic): ugly
  assertion1(dataToSend.size()==1,dataToSend.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending data to master:" <<
             " data[0]=" << dataToSend[0]);
  }
  DataHeap::getInstance().sendData(
      dataToSend.data(), dataToSend.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(workerRank,x,level); // !!! Receive order must be the same in master<->worker comm.
//  _limiter->mergeWithWorkerData(workerRank,x,level); // TODO(Dominic): Revision

  // Receive the information if limiter status has changed
  std::vector<double> receivedData(1); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      receivedData.data(),receivedData.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion(tarch::la::equals(receivedData[0],1.0) ||
            tarch::la::equals(receivedData[0],-1.0)); // TODO(Dominic): ugly

  bool workerLimiterDomainHasChanged = tarch::la::equals(receivedData[0],1.0) ? true : false;
  updateNextLimiterDomainHasChanged(workerLimiterDomainHasChanged); // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Received data from worker:" <<
            " data[0]=" << receivedData[0]);
    logDebug("mergeWithWorkerData(...)","_nextLimiterDomainHasChanged=" << _nextLimiterDomainHasChanged);
  }
}

bool exahype::solvers::LimitingADERDGSolver::hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) {
  #if defined(Asserts) || defined(Debug)
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(_limiter->hasToSendDataToMaster(cellDescriptionsIndex,limiterElement));
  }
  #endif

  return _solver->hasToSendDataToMaster(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToMaster(
        masterRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToMaster(masterRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToMaster(masterRank,x,level);
  _limiter->sendEmptyDataToMaster(masterRank,x,level);
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

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerData( // !!! Receive order must be the same in master<->worker comm.
        workerRank,workerMetadata,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerData(workerRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropWorkerData(workerRank,x,level);// !!! Receive order must be the same in master<->worker comm.
  _limiter->dropWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->sendDataToWorker(workerRank,x,level);
//  _limiter->sendDataToWorker(workerRank,x,level); // TODO(Dominic): Revision

  // TODO(Dominic): Add information that limiter status has been
  // changed for this solver.
  // Send the information to master if limiter status has changed or not
  std::vector<double> dataToSend(0,1);
  dataToSend.push_back(_limiterDomainHasChanged ? 1.0 : -1.0);
  assertion1(dataToSend.size()==1,dataToSend.size());
  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","Sending data to worker:" <<
            " data[0]=" << dataToSend[0]);
  }
  DataHeap::getInstance().sendData(
      dataToSend.data(), dataToSend.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithMasterData(masterRank,x,level); // !!! Receive order must be the same in master<->worker comm.
//  _limiter->mergeWithMasterData(masterRank,x,level); // TODO(Dominic): Revision

  // Receive the information if limiter status has changed
  std::vector<double> receivedData(1); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      receivedData.data(),receivedData.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion(tarch::la::equals(receivedData[0],1.0) ||
            tarch::la::equals(receivedData[0],-1.0)); // TODO(Dominic): ugly

  bool masterLimiterDomainHasChanged = tarch::la::equals(receivedData[0],1.0) ? true : false;
  _limiterDomainHasChanged = masterLimiterDomainHasChanged; // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithMasterData(...)","Received data from worker:" <<
            " data[0]=" << receivedData[0]);
    logDebug("mergeWithMasterData(...)","_limiterDomainHasChanged=" << _limiterDomainHasChanged);
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

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorker(
        workerRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorker(workerRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToWorker(workerRank,x,level);
  _limiter->sendEmptyDataToWorker(workerRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                     masterRank,
    const exahype::MetadataHeap::HeapEntries&     masterMetadata,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->mergeWithMasterData(
      masterRank,masterMetadata,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithMasterData( // !!! Receive order must be the same in master<->worker comm.
        masterRank,masterMetadata,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropMasterData(masterRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropMasterData(masterRank,x,level);
  _limiter->dropMasterData(masterRank,x,level);
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
