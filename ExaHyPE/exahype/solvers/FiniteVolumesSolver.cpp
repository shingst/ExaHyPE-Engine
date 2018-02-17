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
 * \author Dominic E. Charrier, Tobias Weinzierl, Jean-Matthieu Gallard, Fabian Güra
 **/

#include "exahype/solvers/FiniteVolumesSolver.h"

#include <iomanip>
#include <string>
#include <limits>
#include <algorithm>

#include "peano/utils/Loop.h"

#include "tarch/multicore/Lock.h"
#include "peano/datatraversal/TaskSet.h"

#include "exahype/Cell.h"
#include "exahype/Vertex.h"
#include "exahype/VertexOperations.h"

#include "kernels/KernelUtils.h"

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "peano/heap/CompressedFloatingPointNumbers.h"

#include "tarch/multicore/Jobs.h"


namespace {
constexpr const char* tags[]{"solutionUpdate", "stableTimeStepSize"};
}  // namespace

tarch::logging::Log exahype::solvers::FiniteVolumesSolver::_log( "exahype::solvers::FiniteVolumesSolver");

void exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(const int cellDescriptionsIndex) {
  assertion(Heap::getInstance().isValidIndex(cellDescriptionsIndex));

  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    assertion(p.getType()==CellDescription::Type::Cell);

    auto *solver = exahype::solvers::RegisteredSolvers[p.getSolverNumber()];

    FiniteVolumesSolver* fvSolver = nullptr;
    if (solver->getType()==Solver::Type::FiniteVolumes) {
      fvSolver = static_cast<FiniteVolumesSolver*>(solver);
    }
    else if (solver->getType()==Solver::Type::LimitingADERDG) {
      fvSolver =
          static_cast<LimitingADERDGSolver*>(solver)->getLimiter().get();
    }
    assertion(fvSolver!=nullptr);
    p.setType(CellDescription::Type::Erased);
    fvSolver->ensureNoUnnecessaryMemoryIsAllocated(p);
  }

  Heap::getInstance().getData(cellDescriptionsIndex).clear();
}

exahype::solvers::FiniteVolumesSolver::FiniteVolumesSolver(
    const std::string& identifier, int numberOfVariables,
    int numberOfParameters, int nodesPerCoordinateAxis, int ghostLayerWidth,
    double maximumMeshSize, int maximumAdaptiveMeshDepth,
    exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, exahype::solvers::Solver::Type::FiniteVolumes,
             numberOfVariables, numberOfParameters, nodesPerCoordinateAxis,
             maximumMeshSize, maximumAdaptiveMeshDepth,
             timeStepping, std::move(profiler)),
            _previousMinTimeStamp( std::numeric_limits<double>::max() ),
            _previousMinTimeStepSize( std::numeric_limits<double>::max() ),
            _minTimeStamp( std::numeric_limits<double>::max() ),
            _minTimeStepSize( std::numeric_limits<double>::max() ),
            _minNextTimeStepSize( std::numeric_limits<double>::max() ),
            _ghostLayerWidth( ghostLayerWidth ) {
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }

  CompressedDataHeap::getInstance().setName("compressed-data");
}

int exahype::solvers::FiniteVolumesSolver::getDataPerPatch() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis, DIMENSIONS + 0);
}

int exahype::solvers::FiniteVolumesSolver::getGhostLayerWidth() const {
  return _ghostLayerWidth;
}

int exahype::solvers::FiniteVolumesSolver::getGhostDataPerPatch() const {
  return (_numberOfVariables+_numberOfParameters) * power(_nodesPerCoordinateAxis+2*_ghostLayerWidth, DIMENSIONS + 0) - getDataPerPatch();
}

int exahype::solvers::FiniteVolumesSolver::getDataPerPatchFace() const {
  return _ghostLayerWidth*(_numberOfVariables+_numberOfParameters)*power(_nodesPerCoordinateAxis, DIMENSIONS - 1);
}

int exahype::solvers::FiniteVolumesSolver::getDataPerPatchBoundary() const {
  return DIMENSIONS_TIMES_TWO *getDataPerPatchFace();
}

double exahype::solvers::FiniteVolumesSolver::getPreviousMinTimeStepSize() const {
  return _previousMinTimeStepSize;
}

double exahype::solvers::FiniteVolumesSolver::getMinTimeStamp() const {
  return _minTimeStamp;
}

double exahype::solvers::FiniteVolumesSolver::getMinTimeStepSize() const {
  return _minTimeStepSize;
}

void exahype::solvers::FiniteVolumesSolver::updateMinNextTimeStepSize(
    double value) {
  _minNextTimeStepSize = std::min(_minNextTimeStepSize, value);
}

void exahype::solvers::FiniteVolumesSolver::initSolver(
    const double timeStamp,
    const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
    const tarch::la::Vector<DIMENSIONS,double>& domainSize,
    const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize
) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  _coarsestMeshLevel =
      exahype::solvers::Solver::computeMeshLevel(_maximumMeshSize,boundingBoxSize[0]);


  _previousMinTimeStepSize = 0.0;
  _previousMinTimeStamp = timeStamp;

  _minTimeStepSize = 0.0;
  _minTimeStamp = timeStamp;

  _meshUpdateRequest = true; // for the initial mesh refinement
}

bool exahype::solvers::FiniteVolumesSolver::isPerformingPrediction(
    const exahype::State::AlgorithmSection& section) const {
  return false;
}

bool exahype::solvers::FiniteVolumesSolver::isMergingMetadata(
    const exahype::State::AlgorithmSection& section) const {
  return false;
}

void exahype::solvers::FiniteVolumesSolver::synchroniseTimeStepping(
    CellDescription& cellDescription) const {
  switch (_timeStepping) {
    case TimeStepping::Global:
      cellDescription.setPreviousTimeStepSize(_previousMinTimeStamp);
      cellDescription.setPreviousTimeStepSize(_previousMinTimeStepSize);
      cellDescription.setTimeStamp(_minTimeStamp);
      cellDescription.setTimeStepSize(_minTimeStepSize);
      break;
    case TimeStepping::GlobalFixed:
      cellDescription.setPreviousTimeStepSize(_previousMinTimeStamp);
      cellDescription.setPreviousTimeStepSize(_previousMinTimeStepSize);
      cellDescription.setTimeStamp(_minTimeStamp);
      cellDescription.setTimeStepSize(_minTimeStepSize);
      break;
  }
}

void exahype::solvers::FiniteVolumesSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  synchroniseTimeStepping(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::startNewTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _previousMinTimeStepSize  = _minTimeStepSize;
      _minTimeStamp            += _minTimeStepSize;
      _minTimeStepSize          = _minNextTimeStepSize;
      _minNextTimeStepSize      = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _previousMinTimeStepSize  = _minTimeStepSize;
      _minTimeStamp            += _minTimeStepSize;
      _minTimeStepSize          = _minNextTimeStepSize;
      break;
  }

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min
}

void exahype::solvers::FiniteVolumesSolver::startNewTimeStepFused(
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  // n-1
   if ( isFirstIterationOfBatch ) {
     _previousMinTimeStepSize  = _minTimeStepSize;
     _previousMinTimeStamp     = _minTimeStamp;
   }
   // n
   _minTimeStamp            += _minTimeStepSize;
   if ( isLastIterationOfBatch ) {
     switch (_timeStepping) {
       case TimeStepping::Global:
         _minTimeStepSize        = _minNextTimeStepSize;
         _minNextTimeStepSize    = std::numeric_limits<double>::max();
         break;
       case TimeStepping::GlobalFixed:
         _minTimeStepSize        = _minNextTimeStepSize;
         break;
     }

     _minCellSize     = _nextMinCellSize;
     _maxCellSize     = _nextMaxCellSize;
     _nextMinCellSize = std::numeric_limits<double>::max();
     _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min
   }
}

void exahype::solvers::FiniteVolumesSolver::updateTimeStepSizesFused() {
  updateTimeStepSizes();
}

void exahype::solvers::FiniteVolumesSolver::updateTimeStepSizes() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minTimeStepSize          = _minNextTimeStepSize;
      _minNextTimeStepSize      = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minTimeStepSize          = _minNextTimeStepSize;
      break;
  }

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min
}

/**
 * Zero predictor and corrector time step size.
 */
void exahype::solvers::FiniteVolumesSolver::zeroTimeStepSizes() {
  _minTimeStepSize = 0;
  assertionEquals(_minNextTimeStepSize,std::numeric_limits<double>::max());
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStep() {
  switch (_timeStepping) {
    case TimeStepping::Global:
      _minNextTimeStepSize     = std::numeric_limits<double>::max();

      _minTimeStamp    = _previousMinTimeStamp;
      _minTimeStepSize = _previousMinTimeStepSize;

      _previousMinTimeStepSize = std::numeric_limits<double>::max();
      break;
    case TimeStepping::GlobalFixed:
      _minTimeStamp     = _previousMinTimeStamp;
      _minTimeStepSize  = _previousMinTimeStepSize;
      break;
  }
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStepFused() {
  rollbackToPreviousTimeStep();
}

double exahype::solvers::FiniteVolumesSolver::getMinNextTimeStepSize() const {
  return _minNextTimeStepSize;
}

int exahype::solvers::FiniteVolumesSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if (Heap::getInstance().isValidIndex(cellDescriptionsIndex)) {
    int element=0;
    for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
      if (p.getSolverNumber()==solverNumber) {
        return element;
      }
      ++element;
    }
  }
  return NotFound;
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::FiniteVolumesSolver::computeSubcellPositionOfCellOrAncestor(
        const int cellDescriptionsIndex,
        const int element) const {
  SubcellPosition empty;
  return empty;
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::FiniteVolumesSolver::markForRefinement(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  // do nothing
  return false;
}

exahype::solvers::Solver::UpdateStateInEnterCellResult exahype::solvers::FiniteVolumesSolver::updateStateInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const bool initialGrid,
    const int solverNumber) {
  // Fine grid cell based uniform mesh refinement.
  UpdateStateInEnterCellResult result;

  int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (fineGridCellElement==exahype::solvers::Solver::NotFound &&
      tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),getMaximumMeshSize()) &&
      tarch::la::allGreater(coarseGridVerticesEnumerator.getCellSize(),getMaximumMeshSize())) {
    addNewCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
               multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
               solverNumber);

    result._newComputeCellAllocated = true;
    // Fine grid cell based adaptive mesh refinement operations are not implemented.
  } else if (fineGridCellElement!=exahype::solvers::Solver::NotFound) {
    // do nothing
  }

  return result;
}

void exahype::solvers::FiniteVolumesSolver::addNewCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  fineGridCell.addNewCellDescription(
              solverNumber,
              CellDescription::Cell,
              CellDescription::None,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCellDescriptionsIndex,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getVertexPosition());
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

  CellDescription& fineGridCellDescription =
      getCellDescription(fineGridCell.getCellDescriptionsIndex(),fineGridCellElement);
  ensureNecessaryMemoryIsAllocated(fineGridCellDescription);
  fineGridCellDescription.setIsInside(
      exahype::Cell::determineInsideAndOutsideFaces(
          fineGridVertices,fineGridVerticesEnumerator));
}

void exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
    const int cellDescriptionsIndex,
    const int solverNumber,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
    const exahype::records::FiniteVolumesCellDescription::RefinementEvent refinementEvent,
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
    const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion2(static_cast<unsigned int>(solverNumber) < solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  CellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setLevel(level);
  newCellDescription.setRefinementEvent(refinementEvent);
  // newCellDescription.setHelperCellNeedsToStoreFaceData(false); // TODO(Dominic): Add to FV cell descr.

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Default data field indices
  newCellDescription.setSolution(-1);
  newCellDescription.setSolutionAverages(-1);
  newCellDescription.setSolutionCompressed(-1);
  newCellDescription.setPreviousSolution(-1);
  newCellDescription.setPreviousSolutionAverages(-1);
  newCellDescription.setPreviousSolutionCompressed(-1);
  newCellDescription.setExtrapolatedSolution(-1);
  newCellDescription.setExtrapolatedSolutionAverages(-1);
  newCellDescription.setExtrapolatedSolutionCompressed(-1);

  newCellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  Heap::getInstance().getData(cellDescriptionsIndex).push_back(newCellDescription);
  lock.free();
}


void exahype::solvers::FiniteVolumesSolver::ensureNoUnnecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  if (DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
    switch (cellDescription.getType()) {
      case CellDescription::Erased: {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolution()));

          if (cellDescription.getSolution()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getSolution());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getSolution()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressed());
          }
          DataHeap::getInstance().deleteData(cellDescription.getSolutionAverages());

          if (cellDescription.getPreviousSolution()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getPreviousSolution());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getPreviousSolution()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionCompressed());
          }
          DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAverages());

          if (cellDescription.getExtrapolatedSolution()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolution());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getExtrapolatedSolution()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionCompressed());
          }
          DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionAverages());

          cellDescription.setSolution(-1);
          cellDescription.setPreviousSolution(-1);
          cellDescription.setExtrapolatedSolution(-1);

          cellDescription.setSolutionCompressed(-1);
          cellDescription.setPreviousSolutionCompressed(-1);
          cellDescription.setExtrapolatedSolutionCompressed(-1);

          cellDescription.setSolutionAverages(-1);
          cellDescription.setPreviousSolutionAverages(-1);
          cellDescription.setExtrapolatedSolutionAverages(-1);
        lock.free();
      } break;
      case CellDescription::Cell:
        // do nothing
        break;
      default:
        assertionMsg(false,"No other cell description types are supported at the moment!");
        break;
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::ensureNecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  switch (cellDescription.getType()) {
    case CellDescription::Cell:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
        // Allocate volume data
        const int patchSize         = getDataPerPatch()+getGhostDataPerPatch();

        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setSolution(        DataHeap::getInstance().createData( patchSize, patchSize ));
          cellDescription.setPreviousSolution(DataHeap::getInstance().createData( patchSize, patchSize ));

          cellDescription.setSolutionCompressed(-1);
          cellDescription.setPreviousSolutionCompressed(-1);

          cellDescription.setSolutionAverages(
              DataHeap::getInstance().createData( getNumberOfVariables()+getNumberOfParameters(), getNumberOfVariables()+getNumberOfParameters() ) );
          cellDescription.setPreviousSolutionAverages(
              DataHeap::getInstance().createData( getNumberOfVariables()+getNumberOfParameters(), getNumberOfVariables()+getNumberOfParameters() ) );

          // Zero out the solution and previous solution arrays. For our MUSCL-Hancock implementation which
          // does not take the corner neighbours into account e.g., it is important that the values in
          // the corner cells of the first ghost layer are set to zero.
          std::fill_n( DataHeap::getInstance().getData(cellDescription.getSolution()).data(),         patchSize, 0.0 );
          std::fill_n( DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data(), patchSize, 0.0 );

          // Allocate boundary data
          const int patchBoundarySize = getDataPerPatchBoundary();
          cellDescription.setExtrapolatedSolution(DataHeap::getInstance().createData( patchBoundarySize, patchBoundarySize ));
          std::fill_n( DataHeap::getInstance().getData(cellDescription.getExtrapolatedSolution()).data(), patchBoundarySize, 0.0 );

          cellDescription.setExtrapolatedSolutionCompressed(-1);
          cellDescription.setExtrapolatedSolutionAverages( DataHeap::getInstance().createData(
            (getNumberOfVariables()+getNumberOfParameters()) * 2 * DIMENSIONS,
            (getNumberOfVariables()+getNumberOfParameters()) * 2 * DIMENSIONS ) );
        lock.free();
      }
      break;
    case CellDescription::Erased:
      // do nothing
    default:
      assertionMsg(false,"No other cell description types are supported at the moment!");
      break;
  }
}

bool exahype::solvers::FiniteVolumesSolver::attainedStableState(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    const int solverNumber) const {
  return true;
}

bool exahype::solvers::FiniteVolumesSolver::updateStateInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber) {
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

  return fineGridCellElement!=exahype::solvers::Solver::NotFound;
}

void exahype::solvers::FiniteVolumesSolver::finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) {
  // do nothing
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
bool exahype::solvers::FiniteVolumesSolver::evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element) {
  // TODO(Dominic): Not implemented.
  return false;
}

double exahype::solvers::FiniteVolumesSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& p = getCellDescription(cellDescriptionsIndex,element);

  if (p.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
    //         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo
    double* solution = exahype::DataHeap::getInstance().getData(p.getSolution()).data();

    double admissibleTimeStepSize = stableTimeStepSize(solution, p.getSize());
    assertion(!std::isnan(admissibleTimeStepSize));

    // n-1
    p.setPreviousTimeStamp(p.getTimeStamp());
    p.setPreviousTimeStepSize(p.getTimeStepSize());

    // n
    p.setTimeStamp(p.getTimeStamp()+p.getTimeStepSize());
    p.setTimeStepSize(admissibleTimeStepSize);

    return admissibleTimeStepSize;
  } else {
    return std::numeric_limits<double>::max();
  }
}

double exahype::solvers::FiniteVolumesSolver::startNewTimeStepFused(
    const int cellDescriptionsIndex,
    const int element,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if (cellDescription.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
    //         assertion1(cellDescription.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,cellDescription.toString()); // todo
    double* solution = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    double admissibleTimeStepSize = stableTimeStepSize(solution, cellDescription.getSize());
    assertion(!std::isnan(admissibleTimeStepSize));

    // n-1
    if (isFirstIterationOfBatch) {
      cellDescription.setPreviousTimeStamp(cellDescription.getTimeStamp());
      cellDescription.setPreviousTimeStepSize(cellDescription.getTimeStepSize());
    }
    // n
    cellDescription.setTimeStamp(cellDescription.getTimeStamp()+cellDescription.getTimeStepSize());
    if (isLastIterationOfBatch) {
      cellDescription.setTimeStepSize(admissibleTimeStepSize);
    }

    return admissibleTimeStepSize;
  } else {
    return std::numeric_limits<double>::max();
  }
}

double exahype::solvers::FiniteVolumesSolver::updateTimeStepSizesFused(
          const int cellDescriptionsIndex,
          const int element) {
  return updateTimeStepSizes(cellDescriptionsIndex,element);
}

double exahype::solvers::FiniteVolumesSolver::updateTimeStepSizes(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& p = getCellDescription(cellDescriptionsIndex,element);

  if (p.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
    //         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo
    double* solution = exahype::DataHeap::getInstance().getData(p.getSolution()).data();

    double admissibleTimeStepSize = stableTimeStepSize(solution, p.getSize());

    assertion(!std::isnan(admissibleTimeStepSize));
    p.setTimeStepSize(admissibleTimeStepSize);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

void exahype::solvers::FiniteVolumesSolver::zeroTimeStepSizes(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell) {
    cellDescription.setTimeStepSize(0.0);
  }
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  cellDescription.setTimeStamp(cellDescription.getPreviousTimeStamp());
  cellDescription.setTimeStepSize(cellDescription.getPreviousTimeStepSize());

  cellDescription.setPreviousTimeStepSize(std::numeric_limits<double>::max());
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStepFused(
    const int cellDescriptionsIndex,
    const int element) const {
  rollbackToPreviousTimeStep(cellDescriptionsIndex,element);
}

void exahype::solvers::FiniteVolumesSolver::adjustSolution(
    const int cellDescriptionsIndex,
    const int element) {
  // reset helper variables
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell
//      && cellDescription.getRefinementEvent()==CellDescription::None
      ) {
    double* solution = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    adjustSolution(
          solution,
          cellDescription.getOffset()+0.5*cellDescription.getSize(),
          cellDescription.getSize(),
          cellDescription.getTimeStamp(),
          cellDescription.getTimeStepSize());

    for (int i=0; i<getDataPerPatch()+getGhostDataPerPatch(); i++) {
      assertion3(std::isfinite(solution[i]),cellDescription.toString(),"setInitialConditions(...)",i);
    } // Dead code elimination will get rid of this loop if Asserts/Debug flags are not set.
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::FiniteVolumesSolver::fusedTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool vetoSpawnPredictorAsBackgroundThread) {
  updateSolution(cellDescriptionsIndex,element,isFirstIterationOfBatch);

  UpdateResult result;
  result._timeStepSize = startNewTimeStepFused(cellDescriptionsIndex,element,
                                               isFirstIterationOfBatch,isLastIterationOfBatch);
  return result;
}

void exahype::solvers::FiniteVolumesSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    const bool backupPreviousSolution) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  double* newSolution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  double* solution    = DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data();
  if (backupPreviousSolution) {
    std::copy(newSolution,newSolution+getDataPerPatch()+getGhostDataPerPatch(),solution); // Copy (current solution) in old solution field.
  }
  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution");

  //    std::cout << "[pre] solution:" << std::endl;
  //    printFiniteVolumesSolution(cellDescription); // TODO(Dominic): remove
//  if (cellDescriptionsIndex==2516) {
//    std::cout << "[pre] solution:" << std::endl;
//    printFiniteVolumesSolution(cellDescription); // TODO(Dominic): remove
//
//    ADERDGSolver::CellDescription& aderPatch =
//        ADERDGSolver::getCellDescription(cellDescriptionsIndex,cellDescription.getSolverNumber());
//    std::cout << "aderPatch="<<aderPatch.toString() << std::endl;
//  }

  assertion1(cellDescription.getTimeStamp()<std::numeric_limits<double>::max(),cellDescription.toString());
  assertion1(cellDescription.getTimeStepSize()<std::numeric_limits<double>::max(),cellDescription.toString());
  double admissibleTimeStepSize=0;
  if (cellDescription.getTimeStepSize()>0) {
    solutionUpdate(
        newSolution,solution,
        cellDescription.getSize(),cellDescription.getTimeStepSize(),admissibleTimeStepSize);
  }

  // cellDescription.getTimeStepSize() = 0 is an initial condition
  assertion2( tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || !std::isnan(admissibleTimeStepSize), cellDescription.toString(), cellDescriptionsIndex );
  assertion2( tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || !std::isinf(admissibleTimeStepSize), cellDescription.toString(), cellDescriptionsIndex );
  assertion3( tarch::la::equals(cellDescription.getTimeStepSize(),0.0) || admissibleTimeStepSize<std::numeric_limits<double>::max(),
              admissibleTimeStepSize, cellDescription.toString(), cellDescriptionsIndex );

  if ( !tarch::la::equals(cellDescription.getTimeStepSize(), 0.0) && tarch::la::smaller(admissibleTimeStepSize,cellDescription.getTimeStepSize()) ) { //TODO JMG 1.001 factor to prevent same dt computation to throw logerror
    logWarning("updateSolution(...)","Finite volumes solver time step size harmed CFL condition. dt="<<
               cellDescription.getTimeStepSize()<<", dt_adm=" << admissibleTimeStepSize << ". cell=" <<cellDescription.toString());
  }

  adjustSolution(
      newSolution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp()+cellDescription.getTimeStepSize(),
      cellDescription.getTimeStepSize());

//  if (cellDescriptionsIndex==2516) {
//    std::cout << "[post] solution:" << std::endl;
//    printFiniteVolumesSolution(cellDescription); // TODO(Dominic): remove
//  }
  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution");
}

void exahype::solvers::FiniteVolumesSolver::swapSolutionAndPreviousSolution(
    CellDescription& cellDescription) const {
  // Simply swap the heap indices
  const int previousSolution = cellDescription.getPreviousSolution();
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolution(previousSolution);
}


void exahype::solvers::FiniteVolumesSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell)
    return;

  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::restrictToNextParent(
      const int fineGridCellDescriptionsIndex,
      const int fineGridElement,
      const int coarseGridCellDescriptionsIndex,
      const int coarseGridElement) const {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::restrictToTopMostParent(
    const int cellDescriptionsIndex,
    const int element,
    const int parentCellDescriptionsIndex,
    const int parentElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  assertionMsg(false,"Please implement!");
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::mergeNeighboursMetadata(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) const {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  CellDescription& cellDescription1 = getCellDescription(cellDescriptionsIndex1,element1);
  CellDescription& cellDescription2 = getCellDescription(cellDescriptionsIndex2,element2);

//  if (cellDescriptionsIndex1==2516 ||
//      cellDescriptionsIndex2==2516) {
//    std::cout << "cell1: "<< cellDescriptionsIndex1 << "," << cellDescription1.toString() << std::endl;
//    std::cout << "cell2: "<< cellDescriptionsIndex2 << "," << cellDescription2.toString() << std::endl;
//  }

  if (cellDescription1.getType()==CellDescription::Cell ||
      cellDescription2.getType()==CellDescription::Cell) {
    synchroniseTimeStepping(cellDescription1);
    synchroniseTimeStepping(cellDescription2);

    assertion1(cellDescription1.getTimeStamp()<std::numeric_limits<double>::max(),cellDescription1.toString());
    assertion1(cellDescription1.getTimeStepSize()<std::numeric_limits<double>::max(),cellDescription1.toString());
    assertion1(cellDescription2.getTimeStamp()<std::numeric_limits<double>::max(),cellDescription2.toString());
    assertion1(cellDescription2.getTimeStepSize()<std::numeric_limits<double>::max(),cellDescription2.toString());

    peano::datatraversal::TaskSet uncompression(
      [&] () -> void {
        uncompress(cellDescription1);
      },
      [&] () -> void {
        uncompress(cellDescription2);
      },
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	  peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
      true
    );

    double* solution1 = DataHeap::getInstance().getData(cellDescription1.getSolution()).data();
    double* solution2 = DataHeap::getInstance().getData(cellDescription2.getSolution()).data();

    ghostLayerFilling(solution1,solution2,pos2-pos1);
    ghostLayerFilling(solution2,solution1,pos1-pos2);
  }

  return;

  assertionMsg(false,"Not implemented.");
}

void exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData(
    const int                                 cellDescriptionsIndex,
    const int                                 element,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell) {
    synchroniseTimeStepping(cellDescription);
    uncompress(cellDescription);

    double* luh = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    boundaryConditions(
        luh,
        cellDescription.getOffset()+0.5*cellDescription.getSize(),
        cellDescription.getSize(),
        cellDescription.getTimeStamp(),
        cellDescription.getTimeStepSize(),
        posCell,posBoundary);
  }
}

#ifdef Parallel
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerNeighbourCommunication    = 1;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerForkOrJoinCommunication   = 1;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerMasterWorkerCommunication = 0;

void exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),
      cellDescriptionsIndex);

  if (Heap::getInstance().getData(cellDescriptionsIndex).size()>0) {
    Heap::getInstance().sendData(cellDescriptionsIndex,
                                 toRank,x,level,messageType);
  } else {
    sendEmptyCellDescriptions(toRank,messageType,x,level);
  }
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(
    const int                                     toRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::HeapEntries emptyMessage(0);
  Heap::getInstance().sendData(emptyMessage,
      toRank,x,level,messageType);
}

void exahype::solvers::FiniteVolumesSolver::mergeCellDescriptionsWithRemoteData(
    const int                                     fromRank,
    exahype::Cell&                                localCell,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  const int receivedCellDescriptionsIndex =
      Heap::getInstance().createData(0,exahype::solvers::RegisteredSolvers.size());
  Heap::getInstance().receiveData(
      receivedCellDescriptionsIndex,fromRank,x,level,messageType);

  logDebug("mergeCellDescriptionsWithRemoteData(...)","received " <<
      Heap::getInstance().getData(receivedCellDescriptionsIndex).size() <<
      " cell descriptions for cell ("
      "offset="<< x.toString() <<
      "level="<< level << ")");

  Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).clear();
  if ( !Heap::getInstance().getData(receivedCellDescriptionsIndex).empty() ) {
    resetDataHeapIndices(receivedCellDescriptionsIndex,
        multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex);
    assertion1(Heap::getInstance().isValidIndex(localCell.getCellDescriptionsIndex()),
        localCell.getCellDescriptionsIndex());
    Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).reserve(
            Heap::getInstance().getData(receivedCellDescriptionsIndex).size());

    for (auto& pReceived : Heap::getInstance().getData(receivedCellDescriptionsIndex)) {
      auto* solver = exahype::solvers::RegisteredSolvers[pReceived.getSolverNumber()];

      switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        assertionMsg(false,"Solver type not appropriate!");
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter()->ensureNecessaryMemoryIsAllocated(pReceived);
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        static_cast<FiniteVolumesSolver*>(solver)->ensureNecessaryMemoryIsAllocated(pReceived);
        break;
      }

      Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).
          push_back(pReceived);
    }
  }

  Heap::getInstance().deleteData(receivedCellDescriptionsIndex);
  assertion(!Heap::getInstance().isValidIndex(receivedCellDescriptionsIndex));
}

void exahype::solvers::FiniteVolumesSolver::resetDataHeapIndices(
    const int cellDescriptionsIndex,
    const int parentIndex) {
  for (auto& p : Heap::getInstance().getData(cellDescriptionsIndex)) {
    p.setParentIndex(parentIndex);

    // Default field data indices
    p.setPreviousSolution(-1);
    p.setSolution(-1);
  }
}

void exahype::solvers::FiniteVolumesSolver::ensureConsistencyOfParentIndex(
    CellDescription& fineGridCellDescription,
    const int coarseGridCellDescriptionsIndex,
    const int solverNumber) {
  fineGridCellDescription.setParentIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  int coarseGridElement = tryGetElement(coarseGridCellDescriptionsIndex,solverNumber);
  if (coarseGridElement!=exahype::solvers::Solver::NotFound) {
    fineGridCellDescription.setParentIndex(coarseGridCellDescriptionsIndex);
    // In this case, we are not at a master worker boundary only our parent is.
//    fineGridCellDescription.setHasToHoldDataForMasterWorkerCommunication(false);  TODO(Dominic): Introduce to FiniteVolumesCellDescription.def
  }
}

/**
 * Drop cell descriptions received from \p fromRank.
 */
void exahype::solvers::FiniteVolumesSolver::dropCellDescriptions(
    const int                                     fromRank,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(fromRank,x,level,messageType);
}

////////////////////////////////////
// MASTER <=> WORKER
////////////////////////////////////
bool exahype::solvers::FiniteVolumesSolver::prepareMasterCellDescriptionAtMasterWorkerBoundary(
      const int cellDescriptionsIndex,
      const int element) {
  return false;
}

void exahype::solvers::FiniteVolumesSolver::prepareWorkerCellDescriptionAtMasterWorkerBoundary(
      const int cellDescriptionsIndex,
      const int element) {
  // do nothing
}

void
exahype::solvers::FiniteVolumesSolver::appendMasterWorkerCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  for (int i = 0; i < exahype::MasterWorkerCommunicationMetadataPerSolver; ++i) {
    metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
  }
}

void exahype::solvers::FiniteVolumesSolver::mergeWithMasterMetadata(
      const MetadataHeap::HeapEntries& receivedMetadata,
      const int                        cellDescriptionsIndex,
      const int                        element) {
  // do nothing
}

bool exahype::solvers::FiniteVolumesSolver::mergeWithWorkerMetadata(
      const MetadataHeap::HeapEntries& receivedMetadata,
      const int                        cellDescriptionsIndex,
      const int                        element) {
  return false;
}

void exahype::solvers::FiniteVolumesSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
      element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& p = Heap::getInstance().getData(cellDescriptionsIndex)[element];

  if (p.getType()==CellDescription::Cell) {
    double* solution = DataHeap::getInstance().getData(p.getSolution()).data();

    logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << p.getSolverNumber() << " sent to rank "<<toRank<<
        ", cell: "<< x << ", level: " << level);

    DataHeap::getInstance().sendData(
        solution, getDataPerPatch()+getGhostDataPerPatch(), toRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
  }
}


void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  for(int sends=0; sends<DataMessagesPerForkOrJoinCommunication; ++sends)
    DataHeap::getInstance().sendData(
        exahype::EmptyDataHeapMessage, toRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}


void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  auto& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  assertion4(tarch::la::equals(x,cellDescription.getOffset()+0.5*cellDescription.getSize()),x,cellDescription.getOffset()+0.5*cellDescription.getSize(),level,cellDescription.getLevel());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);

  if (cellDescription.getType()==CellDescription::Cell) {
    logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
             ", cell: "<< x << ", level: " << level);

    DataHeap::getInstance().getData(cellDescription.getSolution()).clear();
    DataHeap::getInstance().receiveData(
        cellDescription.getSolution(),fromRank,x,level,
        peano::heap::MessageType::ForkOrJoinCommunication);
  }
}

void exahype::solvers::FiniteVolumesSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  for(int receives=0; receives<DataMessagesPerForkOrJoinCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::ForkOrJoinCommunication);
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void
exahype::solvers::FiniteVolumesSolver::appendNeighbourCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  const int element = tryGetElement(cellDescriptionsIndex,solverNumber);

  if (element!=exahype::solvers::Solver::NotFound)  {
    CellDescription& cellDescription =
        getCellDescription(cellDescriptionsIndex,element);
    metadata.push_back(static_cast<int>(cellDescription.getType()));
    metadata.push_back(exahype::InvalidMetadataEntry);
    metadata.push_back(exahype::InvalidMetadataEntry);
    metadata.push_back(exahype::InvalidMetadataEntry);
  } else {
    for (int i = 0; i < exahype::NeighbourCommunicationMetadataPerSolver; ++i) {
      metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::mergeWithNeighbourMetadata(
    const exahype::MetadataHeap::HeapEntries& metadata,
    const tarch::la::Vector<DIMENSIONS, int>& src,
    const tarch::la::Vector<DIMENSIONS, int>& dest,
    const int cellDescriptionsIndex,
    const int element) const {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  assertion( tarch::la::countEqualEntries(src,dest)==(DIMENSIONS-1) );

  const int direction    = tarch::la::equalsReturnIndex(src, dest);
  const int orientation  = (1 + dest(direction) - src(direction))/2;
  const int faceIndex    = 2*direction+orientation;

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (holdsFaceData(cellDescription.getType())) {
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));

    const int numberOfFaceDof = getDataPerPatchFace();
//    const int boundaryLayerToSendIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
    double* luhbnd = DataHeap::getInstance().getData(cellDescription.getExtrapolatedSolution()).data()
        + (faceIndex * numberOfFaceDof);
    const double* luh = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    boundaryLayerExtraction(luhbnd,luh,dest-src);

    logDebug(
        "sendDataToNeighbour(...)",
        "send "<<DataMessagesPerNeighbourCommunication<<" arrays to rank " <<
        toRank << " for cell="<<cellDescription.getOffset()
//        << "and face=" << faceIndex
        << " from vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(cellDescriptionsIndex) <<
        ", src=" << src << ", dest=" << dest <<
        ", size="<<DataHeap::getInstance().getData(cellDescription.getExtrapolatedSolution()).size()
//        << ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    // Send order: minMax,lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd,minMax
    DataHeap::getInstance().sendData(
        luhbnd, numberOfFaceDof, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
    // TODO(Dominic): If anarchic time stepping send the time step over too.

//    DataHeap::getInstance().deleteData(boundaryLayerToSendIndex,true);
  } else {
    sendEmptyDataToNeighbour(toRank,x,level);
  }
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        exahype::EmptyDataHeapMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const MetadataHeap::HeapEntries&             neighbourMetadata,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  if (tarch::la::countEqualEntries(src,dest)!=(DIMENSIONS-1)) {
    return; // We only consider faces; no corners.
  }

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  synchroniseTimeStepping(cellDescription);
  CellDescription::Type neighbourType =
      static_cast<CellDescription::Type>(neighbourMetadata[exahype::NeighbourCommunicationMetadataCellType].getU());

  const int direction    = tarch::la::equalsReturnIndex(src, dest);
  const int orientation  = (1 + src(direction) - dest(direction))/2;
  const int faceIndex    = 2*direction+orientation;

  // TODO(Dominic): Add to docu: We only perform a Riemann solve if a Cell is involved.
  // Solving Riemann problems at a Ancestor Ancestor boundary might lead to problems
  // if one Ancestor is just used for restriction.
  if(neighbourType==CellDescription::Type::Cell || cellDescription.getType()==CellDescription::Type::Cell){
    assertion1(holdsFaceData(neighbourType),neighbourType);
    assertion1(holdsFaceData(cellDescription.getType()),cellDescription.toString());

    assertion4(!cellDescription.getNeighbourMergePerformed(faceIndex),
        faceIndex,cellDescriptionsIndex,cellDescription.getOffset().toString(),cellDescription.getLevel());
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));

    logDebug(
        "mergeWithNeighbourData(...)", "receive "<<DataMessagesPerNeighbourCommunication<<" arrays from rank " <<
        fromRank << " for vertex x=" << x << ", level=" << level <<
        ", src type=" << cellDescription.getType() <<
        ", src=" << src << ", dest=" << dest <<
        ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    // TODO(Dominic): If anarchic time stepping, receive the time step too.
    //
    // Copy the received boundary layer into a ghost layer of the solution.
    // TODO(Dominic): Pipe it directly through the Riemann solver if
    // we only use the Godunov method and not higher-order FVM methods.
    // For methods that are higher order in time, e.g., MUSCL-Hancock, we usually need
    // corner neighbours. This is why we currently adapt a GATHER-UPDATE algorithm
    // instead of a SOLVE RIEMANN PROBLEM AT BOUNDARY-UPDATE INTERIOR scheme.
    const int numberOfFaceDof      = getDataPerPatchFace();
//    const int receivedBoundaryLayerIndex = DataHeap::getInstance().createData(0, numberOfFaceDof);
    double* luhbnd = DataHeap::getInstance().getData(cellDescription.getExtrapolatedSolution()).data()
        + (faceIndex * numberOfFaceDof);
//    double* luhbnd = DataHeap::getInstance().getData(receivedBoundaryLayerIndex).data();
//    assertion(DataHeap::getInstance().getData(receivedBoundaryLayerIndex).empty());


    // Send order: minMax,lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd,minMax
    DataHeap::getInstance().receiveData(luhbnd, numberOfFaceDof, fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

    logDebug(
        "mergeWithNeighbourData(...)", "[pre] solve Riemann problem with received data." <<
        " cellDescription=" << cellDescription.toString() <<
        ",faceIndexForCell=" << faceIndex <<
        ",normalOfExchangedFac=" << direction <<
        ",x=" << x.toString() << ", level=" << level <<
        ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    double* luh = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
    ghostLayerFillingAtBoundary(luh,luhbnd,src-dest);

//    DataHeap::getInstance().deleteData(receivedBoundaryLayerIndex,true);
  } else  {
    logDebug(
        "mergeWithNeighbourData(...)", "drop one array from rank " <<
        fromRank << " for vertex x=" << x << ", level=" << level <<
        ", src type=" << multiscalelinkedcell::indexToString(cellDescriptionsIndex) <<
        ", src=" << src << ", dest=" << dest
        << ", counter=" << cellDescription.getFaceDataExchangeCounter(faceIndex)
    );

    dropNeighbourData(fromRank,src,dest,x,level);
  }
}

void exahype::solvers::FiniteVolumesSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
    DataHeap::getInstance().receiveData(
        fromRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////

/*
 * At the time of sending data to the master,
 * we have already performed a time step update locally
 * on the rank. We thus need to communicate the
 * current min predictor time step size to the master.
 * The next min predictor time step size is
 * already reset locally to the maximum double value.
 *
 * However on the master's side, we need to
 * merge the received time step size with
 * the next min predictor time step size since
 * the master has not yet performed a time step update
 * (i.e. called TimeStepSizeComputation::endIteration()).
 */
void exahype::solvers::FiniteVolumesSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::HeapEntries timeStepDataToReduce(0,4);
  timeStepDataToReduce.push_back(_minTimeStepSize);
  timeStepDataToReduce.push_back(_meshUpdateRequest ? 1.0 : -1.0);
  timeStepDataToReduce.push_back(_minCellSize);
  timeStepDataToReduce.push_back(_maxCellSize);

  assertion1(timeStepDataToReduce.size()==4,timeStepDataToReduce.size());
  assertion1(std::isfinite(timeStepDataToReduce[0]),timeStepDataToReduce[0]);
  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
        tarch::parallel::Node::getInstance().getRank());
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending time step data: " <<
        "data[0]=" << timeStepDataToReduce[0] <<
        ",data[1]=" << timeStepDataToReduce[1] <<
        ",data[2]=" << timeStepDataToReduce[2] <<
        ",data[3]=" << timeStepDataToReduce[3]);
  }

  DataHeap::getInstance().sendData(
      timeStepDataToReduce.data(), timeStepDataToReduce.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

/**
 * At the time of the merging,
 * the workers and the master have already performed
 * at local update of the next predictor time step size
 * and of the predictor time stamp.
 * We thus need to minimise over both quantities.
 */
void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  DataHeap::HeapEntries receivedTimeStepData(4);

  if (true || tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data [pre] from rank " << workerRank);
  }

  DataHeap::getInstance().receiveData(
      receivedTimeStepData.data(),receivedTimeStepData.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(receivedTimeStepData.size()==4,receivedTimeStepData.size());
  assertion1(receivedTimeStepData[0]>=0,receivedTimeStepData[0]);
  assertion1(std::isfinite(receivedTimeStepData[0]),receivedTimeStepData[0]);
  // The master solver has not yet updated its minNextTimeStepSize.
  // Thus it does not yet equal MAX_DOUBLE.

  int index=0;
  _minNextTimeStepSize      = std::min( _minNextTimeStepSize, receivedTimeStepData[index++] );
  _nextMeshUpdateRequest |= ( receivedTimeStepData[index++] > 0 ) ? true : false;
  _nextMinCellSize          = std::min( _nextMinCellSize, receivedTimeStepData[index++] );
  _nextMaxCellSize          = std::max( _nextMaxCellSize, receivedTimeStepData[index++] );

  if (true || tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data: " <<
             "data[0]=" << receivedTimeStepData[0] <<
             ",data[1]=" << receivedTimeStepData[1] <<
             ",data[2]=" << receivedTimeStepData[2] <<
             ",data[3]=" << receivedTimeStepData[3] );

    logDebug("mergeWithWorkerData(...)","Updated time step fields: " <<
             "_minNextTimeStepSize="     << _minNextTimeStepSize <<
             "_nextMeshUpdateRequest=" << _nextMeshUpdateRequest <<
             ",_nextMinCellSize="        << _nextMinCellSize <<
             ",_nextMaxCellSize="        << _nextMaxCellSize);
  }
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // do nothing - limiter does not depend on solution from master for any tree cut
}

void exahype::solvers::FiniteVolumesSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
      element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  // do nothing - limiter does not depend on solution from master for any tree cut
}

void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerData(
    const int                                     workerRank,
    const exahype::MetadataHeap::HeapEntries&     workerMetadata,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  logDebug("mergeWithWorkerData(...)","Merge with worker data from rank "<<workerRank<<
             ", cell: "<< x << ", level: " << level);

  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  // do nothing - limiter does not depend on solution from master for any tree cut
}

void exahype::solvers::FiniteVolumesSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  logDebug("dropWorkerData(...)","Dropping worker data from rank "<<workerRank<<
               ", cell: "<< x << ", level: " << level);

  // do nothing - limiter does not depend on solution from master for any tree cut
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  std::vector<double> timeStepDataToSend(0,4);
  timeStepDataToSend.push_back(_minTimeStamp); // TODO(Dominic): Append previous time step size
  timeStepDataToSend.push_back(_minTimeStepSize);

  timeStepDataToSend.push_back(_minCellSize);
  timeStepDataToSend.push_back(_maxCellSize);

  assertion1(timeStepDataToSend.size()==4,timeStepDataToSend.size());
  assertion1(std::isfinite(timeStepDataToSend[0]),timeStepDataToSend[0]);
  assertion1(std::isfinite(timeStepDataToSend[1]),timeStepDataToSend[1]);

  if (_timeStepping==TimeStepping::Global) {
    assertionEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
        tarch::parallel::Node::getInstance().getRank());
  }

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","Broadcasting time step data: " <<
        " data[0]=" << timeStepDataToSend[0] <<
        ",data[1]=" << timeStepDataToSend[1] <<
        ",data[2]=" << timeStepDataToSend[2] <<
        ",data[3]=" << timeStepDataToSend[3]);
    logDebug("sendDataWorker(...)","_minNextTimeStepSize="<<_minNextTimeStepSize);
  }

  DataHeap::getInstance().sendData(
      timeStepDataToSend.data(), timeStepDataToSend.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithMasterData(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  std::vector<double> receivedTimeStepData(4);
  DataHeap::getInstance().receiveData(
      receivedTimeStepData.data(),receivedTimeStepData.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion1(receivedTimeStepData.size()==4,receivedTimeStepData.size());

  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
        _minNextTimeStepSize);
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithMasterData(...)","Received time step data: " <<
        "data[0]="  << receivedTimeStepData[0] <<
        ",data[1]=" << receivedTimeStepData[1] <<
        ",data[2]=" << receivedTimeStepData[2] <<
        ",data[3]=" << receivedTimeStepData[3]);
  }

  _minTimeStamp    = receivedTimeStepData[0];
  _minTimeStepSize = receivedTimeStepData[1];

  _minCellSize              = receivedTimeStepData[2];
  _maxCellSize              = receivedTimeStepData[3];
}

bool exahype::solvers::FiniteVolumesSolver::hasToSendDataToMaster(
    const int cellDescriptionsIndex,
    const int element) const {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  return false;
}

void exahype::solvers::FiniteVolumesSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level){
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
             element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  // do nothing - limiter does not depend on solution from master for any tree cut
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // do nothing - limiter does not depend on solution from master for any tree cut
}

void exahype::solvers::FiniteVolumesSolver::receiveDataFromMaster(
     const int                                    masterRank,
     std::deque<int>&                             receivedDataHeapIndices,
     const tarch::la::Vector<DIMENSIONS, double>& x,
     const int                                    level) const {
  // do nothing - limiter does not depend on solution from master for any tree cut
}

 void exahype::solvers::FiniteVolumesSolver::mergeWithMasterData(
     const MetadataHeap::HeapEntries&             masterMetadata,
     std::deque<int>&                             receivedDataHeapIndices,
     const int                                    cellDescriptionsIndex,
     const int                                    element) const {
   // do nothing - limiter does not depend on solution from master for any tree cut
 }

void exahype::solvers::FiniteVolumesSolver::dropMasterData(
     std::deque<int>& heapIndices) const {
  // do nothing - limiter does not depend on solution from master for any tree cut
}
#endif

void exahype::solvers::FiniteVolumesSolver::validateNoNansInFiniteVolumesSolution(
    CellDescription& cellDescription,const int cellDescriptionsIndex,const char* methodTrace)  const {
  #if defined(Asserts)
  double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  #endif

  dfor(i,_nodesPerCoordinateAxis+_ghostLayerWidth) {
    if (tarch::la::allSmaller(i,_nodesPerCoordinateAxis+_ghostLayerWidth)
    && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
      for (int unknown=0; unknown < _numberOfVariables; unknown++) {
        #if defined(Asserts)
        int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
        #endif // cellDescription.getTimeStepSize()==0.0 is an initial condition
        assertion7(tarch::la::equals(cellDescription.getTimeStepSize(),0.0)  || std::isfinite(solution[iScalar]),
                   cellDescription.toString(),cellDescriptionsIndex,solution[iScalar],i.toString(),
                   _nodesPerCoordinateAxis,_ghostLayerWidth,
                   methodTrace);
      }
    }
  } // Dead code elimination should get rid of this loop if Asserts is not set.
}

void exahype::solvers::FiniteVolumesSolver::printFiniteVolumesSolution(
    CellDescription& cellDescription)  const {
  #if DIMENSIONS==2
  double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();

  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    dfor(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth) {
      int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
      std::cout << std::setprecision(3) << solution[iScalar] << ",";
      if (i(0)==_nodesPerCoordinateAxis+2*_ghostLayerWidth-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #else
  double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();

  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    dfor(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth) {
      int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
      std::cout << solution[iScalar] << ",";
      if (i(0)==_nodesPerCoordinateAxis+2*_ghostLayerWidth-1 &&
          i(1)==_nodesPerCoordinateAxis+2*_ghostLayerWidth-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #endif
}

void exahype::solvers::FiniteVolumesSolver::printFiniteVolumesBoundaryLayer(const double* luhbnd)  const {
  #if DIMENSIONS==2
  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
    std::cout <<  "unknown=" << unknown << std::endl;
    for(int i=0; i<_nodesPerCoordinateAxis; ++i) {
      int iScalar = i*_numberOfVariables+unknown;
      std::cout << luhbnd[iScalar] << ",";
      if (i==_nodesPerCoordinateAxis-1) {
        std::cout << std::endl;
      }
    }
  }
  std::cout <<  "}" << std::endl;
  #else
  for (int unknown=0; unknown < _numberOfVariables; unknown++) {
      std::cout <<  "unknown=" << unknown << std::endl;
      for(int j=0; j<_nodesPerCoordinateAxis; ++j) {
      for(int i=0; i<_nodesPerCoordinateAxis; ++i) {
        int iScalar = (j*_nodesPerCoordinateAxis+i)*_numberOfVariables+unknown;
        std::cout << luhbnd[iScalar] << ",";
        if (j==_nodesPerCoordinateAxis-1 &&
            i==_nodesPerCoordinateAxis-1) {
          std::cout << std::endl;
        }
      }
      }
    }
    std::cout <<  "}" << std::endl;
  #endif
}

std::string exahype::solvers::FiniteVolumesSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::FiniteVolumesSolver::toString (std::ostream& out) const {
  out << "(";
  out << "_identifier:" << _identifier;
  out << ",";
  out << "_type:" << exahype::solvers::Solver::toString(_type);
  out << ",";
  out << "_numberOfVariables:" << _numberOfVariables;
  out << ",";
  out << "_numberOfParameters:" << _numberOfParameters;
  out << ",";
  out << "_nodesPerCoordinateAxis:" << _nodesPerCoordinateAxis;
  out << ",";
  out << "_maximumMeshSize:" << _maximumMeshSize;
  out << ",";
  out << "_timeStepping:" << exahype::solvers::Solver::toString(_timeStepping); // only solver attributes
  out << ",";
  out << "_dataPerPatchFace:" << getDataPerPatchFace();
  out << ",";
  out << "_dataPerPatchBoundary:" << getDataPerPatchBoundary();
  out << ",";
  out << "_dataPerPatch:" << getDataPerPatch();
  out << ",";
  out << "_previousMinTimeStepSize:" << _previousMinTimeStepSize;
  out << ",";
  out << "_minTimeStamp:" << _minTimeStamp;
  out << ",";
  out << "_minTimeStepSize:" << _minTimeStepSize;
  out << ",";
  out << "_nextMinTimeStepSize:" << _minNextTimeStepSize;
  out <<  ")";
}


void exahype::solvers::FiniteVolumesSolver::putUnknownsIntoByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  assertion( cellDescription.getPreviousSolutionCompressed()==-1 );
  assertion( cellDescription.getSolutionCompressed()==-1 );
  assertion( cellDescription.getExtrapolatedSolutionCompressed()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfExtrapolatedSolution;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolution() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> void  { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).data(),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy
      );},
    [&] () -> void  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getSolution() ).data(),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy
      );},
    [&]() -> void  { compressionOfExtrapolatedSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() ).data(),
      getDataPerPatchBoundary(),
      CompressionAccuracy
      );},
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
    true
  );

  assertion(1<=compressionOfPreviousSolution);
  assertion(1<=compressionOfSolution);
  assertion(1<=compressionOfExtrapolatedSolution);

  assertion(compressionOfPreviousSolution<=7);
  assertion(compressionOfSolution<=7);
  assertion(compressionOfExtrapolatedSolution<=7);

  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> void {
      cellDescription.setBytesPerDoFInPreviousSolution(compressionOfPreviousSolution);
      if (compressionOfPreviousSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        cellDescription.setPreviousSolutionCompressed( CompressedDataHeap::getInstance().createData(0,0) );
        assertion1(
          cellDescription.getPreviousSolutionCompressed()>=0,
          cellDescription.toString()
        );
        lock.free();

        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
        tearApart(numberOfEntries, cellDescription.getPreviousSolution(), cellDescription.getPreviousSolutionCompressed(), compressionOfPreviousSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getPreviousSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getPreviousSolution(), true );
        cellDescription.setPreviousSolution( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).size() * 8.0;
        lock.free();
        #endif
      }
    },
    [&]() -> void {
      cellDescription.setBytesPerDoFInSolution(compressionOfSolution);
      if (compressionOfSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        cellDescription.setSolutionCompressed(CompressedDataHeap::getInstance().createData(0,0));
        assertion1( cellDescription.getSolutionCompressed()>=0, cellDescription.getSolutionCompressed() );
        lock.free();

        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();

        tearApart(numberOfEntries, cellDescription.getSolution(), cellDescription.getSolutionCompressed(), compressionOfSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getSolution(), true );
        cellDescription.setSolution( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getSolution() ).size() * 8.0;
        lock.free();
        #endif
      }
    },
    [&]() -> void {
      cellDescription.setBytesPerDoFInExtrapolatedSolution(compressionOfExtrapolatedSolution);
      if (compressionOfExtrapolatedSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        cellDescription.setExtrapolatedSolutionCompressed( CompressedDataHeap::getInstance().createData(0,0) );
        assertion( cellDescription.getExtrapolatedSolutionCompressed()>=0 );
        lock.free();

        const int numberOfEntries = getDataPerPatchBoundary();
        tearApart(numberOfEntries, cellDescription.getExtrapolatedSolution(), cellDescription.getExtrapolatedSolutionCompressed(), compressionOfExtrapolatedSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() ).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedSolutionCompressed() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedSolution(), true );
        cellDescription.setExtrapolatedSolution( -1 );
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        PipedUncompressedBytes += DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() ).size() * 8.0;
        PipedCompressedBytes   += DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() ).size() * 8.0;
        lock.free();
        #endif
      }
    },
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
	peano::datatraversal::TaskSet::TaskType::Background,
    true
  );
}


void exahype::solvers::FiniteVolumesSolver::pullUnknownsFromByteStream(
    CellDescription& cellDescription) const {
  assertion(CompressionAccuracy>0.0);

  #if !defined(ValidateCompressedVsUncompressedData)
  const int dataPerCell     = getDataPerPatch() + getGhostDataPerPatch();
  const int dataPerBoundary = getDataPerPatchBoundary();

  {
    tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
    cellDescription.setPreviousSolution( DataHeap::getInstance().createData(
        dataPerCell,         dataPerCell,
      DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setSolution( DataHeap::getInstance().createData(
        dataPerCell,         dataPerCell,
      DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setExtrapolatedSolution( DataHeap::getInstance().createData(
        dataPerBoundary, dataPerBoundary,
        DataHeap::Allocation::UseOnlyRecycledEntries) );
    lock.free();

    if (cellDescription.getPreviousSolution()==-1) {
      waitUntilAllBackgroundJobsHaveTerminated();
      lock.lock();
      cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolution()==-1) {
      waitUntilAllBackgroundJobsHaveTerminated();
      lock.lock();
      cellDescription.setSolution( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedSolution()==-1) {
      waitUntilAllBackgroundJobsHaveTerminated();
      lock.lock();
      cellDescription.setExtrapolatedSolution( DataHeap::getInstance().createData(dataPerBoundary, dataPerBoundary ) );
      lock.free();
    }
  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolution() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressed() ),
      cellDescription.getPreviousSolutionCompressed()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressed() ),
    cellDescription.getSolutionCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionCompressed() ),
    cellDescription.getExtrapolatedSolutionCompressed()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolution() ), cellDescription.getPreviousSolution());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressed() ));
        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
        glueTogether(numberOfEntries, cellDescription.getPreviousSolution(), cellDescription.getPreviousSolutionCompressed(), cellDescription.getBytesPerDoFInPreviousSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionCompressed(), true );
        cellDescription.setPreviousSolutionCompressed( -1 );
        lock.free();
      }
    },
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getSolution() ), cellDescription.getSolution() );
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressed() ));
        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
        glueTogether(numberOfEntries, cellDescription.getSolution(), cellDescription.getSolutionCompressed(), cellDescription.getBytesPerDoFInSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getSolutionCompressed(), true );
        cellDescription.setSolutionCompressed( -1 );
        lock.free();
      }
    },
    [&]() -> void {
      if (cellDescription.getBytesPerDoFInExtrapolatedSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolution() ), cellDescription.getExtrapolatedSolution());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionCompressed() ));
        const int numberOfEntries = getDataPerPatchBoundary();
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedSolution(), cellDescription.getExtrapolatedSolutionCompressed(), cellDescription.getBytesPerDoFInExtrapolatedSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedSolutionCompressed(), true );
        cellDescription.setExtrapolatedSolutionCompressed( -1 );
        lock.free();
      }
    },
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	true
  );
}


void exahype::solvers::FiniteVolumesSolver::compress(CellDescription& cellDescription) {
  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if (SpawnCompressionAsBackgroundJob) {
      cellDescription.setCompressionState(CellDescription::CurrentlyProcessed);

      tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
      _NumberOfTriggeredTasks++;
      lock.free();

      CompressionTask myTask( *this, cellDescription );
      peano::datatraversal::TaskSet spawnedSet( myTask,peano::datatraversal::TaskSet::TaskType::Background );
    }
    else {
      determineUnknownAverages(cellDescription);
      computeHierarchicalTransform(cellDescription,-1.0);
      putUnknownsIntoByteStream(cellDescription);
      cellDescription.setCompressionState(CellDescription::Compressed);
    }
  }
}


void exahype::solvers::FiniteVolumesSolver::uncompress(CellDescription& cellDescription) const {
  #ifdef SharedMemoryParallelisation
  bool madeDecision = CompressionAccuracy<=0.0;
  bool uncompress   = false;

  while (!madeDecision) {
    tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
    madeDecision = cellDescription.getCompressionState() != CellDescription::CurrentlyProcessed;
    uncompress   = cellDescription.getCompressionState() == CellDescription::Compressed;
    if (uncompress) {
      cellDescription.setCompressionState( CellDescription::CurrentlyProcessed );
    }
    lock.free();

    tarch::multicore::jobs::processBackgroundJobs();
  }
  #else
  bool uncompress = CompressionAccuracy>0.0
      && cellDescription.getCompressionState() == CellDescription::Compressed;
  #endif

  if (uncompress) {
    pullUnknownsFromByteStream(cellDescription);
    computeHierarchicalTransform(cellDescription,1.0);

    tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
    cellDescription.setCompressionState(CellDescription::Uncompressed);
    lock.free();
  }

  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
}


void exahype::solvers::FiniteVolumesSolver::determineUnknownAverages(
  CellDescription& cellDescription) const {

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolution()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolution()), cellDescription.toString() );

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionAverages()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionAverages()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionAverages()), cellDescription.toString() );

  const int dataPerSubcell   = getNumberOfParameters()+getNumberOfVariables();
  const int subcellsPerPatch = (getDataPerPatch()+getGhostDataPerPatch())/ dataPerSubcell;
  const int subcellsPerFace  = (getDataPerPatchFace()) / dataPerSubcell;

  auto& solutionAverages             = DataHeap::getInstance().getData( cellDescription.getSolutionAverages() );
  auto& previousSolutionAverage      = DataHeap::getInstance().getData( cellDescription.getPreviousSolutionAverages() );
  auto& extrapolatedSolutionAverages = DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolutionAverages() );

  // patch data
  kernels::idx2 idx_patchData    (subcellsPerPatch,dataPerSubcell);
  for (int i=0; i<subcellsPerPatch; i++) {
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      solutionAverages[variableNumber]        += DataHeap::getInstance().getData( cellDescription.getSolution() )        [idx_patchData(i,variableNumber)];
      previousSolutionAverage[variableNumber] += DataHeap::getInstance().getData( cellDescription.getPreviousSolution() )[idx_patchData(i,variableNumber)];
    }
  }
  for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
    solutionAverages[variableNumber]        = solutionAverages[variableNumber]        / (double) subcellsPerPatch;
    previousSolutionAverage[variableNumber] = previousSolutionAverage[variableNumber] / (double) subcellsPerPatch;
  }

  // face data
  kernels::idx2 idx_faceDataAvg(DIMENSIONS_TIMES_TWO,dataPerSubcell);
  kernels::idx3 idx_faceData   (DIMENSIONS_TIMES_TWO,subcellsPerFace,dataPerSubcell);
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<subcellsPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
        extrapolatedSolutionAverages[idx_faceDataAvg(face,variableNumber)] +=
            DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() )[idx_faceData(face,i,variableNumber)];
      }
    }
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      extrapolatedSolutionAverages[idx_faceDataAvg(face,variableNumber)] =
          extrapolatedSolutionAverages[idx_faceDataAvg(face,variableNumber)] / (double) subcellsPerFace;
    }
  }
}


void exahype::solvers::FiniteVolumesSolver::computeHierarchicalTransform(
    CellDescription& cellDescription, double sign) const {
  const int dataPerSubcell   = getNumberOfParameters()+getNumberOfVariables();
  const int subcellsPerPatch = (getDataPerPatch()+getGhostDataPerPatch())/ dataPerSubcell;
  const int subcellsPerFace = (getDataPerPatchFace()) / dataPerSubcell;

  // patch data
  kernels::idx2 idx_patchData    (subcellsPerPatch,dataPerSubcell);
  for (int i=0; i<subcellsPerPatch; i++) {
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      DataHeap::getInstance().getData( cellDescription.getSolution() )
          [idx_patchData(i,variableNumber)] += sign * DataHeap::getInstance().getData( cellDescription.getSolutionAverages() )[variableNumber];

      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() )
          [idx_patchData(i,variableNumber)] += sign * DataHeap::getInstance().getData( cellDescription.getPreviousSolutionAverages() )[variableNumber];
    }
  }

  // face data
  kernels::idx2 idx_faceDataAvg(DIMENSIONS_TIMES_TWO,dataPerSubcell);
  kernels::idx3 idx_faceData   (DIMENSIONS_TIMES_TWO,subcellsPerFace,dataPerSubcell);
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<subcellsPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolution() ), cellDescription.toString() );
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionAverages() ), cellDescription.toString() );
        DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() )[idx_faceData(face,i,variableNumber)] +=
              sign * DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolutionAverages() )[idx_faceDataAvg(face,variableNumber)];
      }
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (
    cellDescription.getType()==CellDescription::Type::Cell
  ) {
    uncompress(cellDescription);
  }
}

void exahype::solvers::FiniteVolumesSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (
      cellDescription.getType()==CellDescription::Type::Cell
      &&
      CompressionAccuracy>0.0
      &&
      cellDescription.getCompressionState() == CellDescription::Uncompressed
    ) {
    compress(cellDescription);
  }
}

exahype::solvers::FiniteVolumesSolver::CompressionTask::CompressionTask(
  FiniteVolumesSolver&  solver,
  CellDescription&      cellDescription
):
  _solver(solver),
  _cellDescription(cellDescription) {
}


bool exahype::solvers::FiniteVolumesSolver::CompressionTask::operator()() {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  _cellDescription.setCompressionState(CellDescription::Compressed);
  // @todo raus
  _NumberOfTriggeredTasks--;
  assertion( _NumberOfTriggeredTasks>=0 );
  lock.free();
  return false;
}
