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
    const std::string& identifier,
    const int numberOfVariables,
    const int numberOfParameters,
    const int basisSize,
    const int ghostLayerWidth,
    const double maximumMeshSize,
    const exahype::solvers::Solver::TimeStepping timeStepping,
    std::unique_ptr<profilers::Profiler> profiler)
    : Solver(identifier, exahype::solvers::Solver::Type::FiniteVolumes,
             numberOfVariables, numberOfParameters, basisSize,
             maximumMeshSize, 0,
             timeStepping, std::move(profiler)),
            _previousMinTimeStamp( std::numeric_limits<double>::max() ),
            _previousMinTimeStepSize( std::numeric_limits<double>::max() ),
            _minTimeStamp( std::numeric_limits<double>::max() ),
            _minTimeStepSize( std::numeric_limits<double>::max() ),
            _minNextTimeStepSize( std::numeric_limits<double>::max() ),
            _ghostLayerWidth( ghostLayerWidth ),
            _meshUpdateEvent( MeshUpdateEvent::None ),
            _nextMeshUpdateEvent( MeshUpdateEvent::None ) {
  // register tags with profiler
  for (const char* tag : tags) {
    _profiler->registerTag(tag);
  }

  // TODO(WORKAROUND)
  const int dataPerFace = getBndFaceSize();
  _invalidExtrapolatedSolution.resize(dataPerFace);
  std::fill_n(_invalidExtrapolatedSolution.data(),_invalidExtrapolatedSolution.size(),-1);
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

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::FiniteVolumesSolver::getNextMeshUpdateEvent() const {
  return _nextMeshUpdateEvent;
}

void exahype::solvers::FiniteVolumesSolver::setNextMeshUpdateEvent() {
  _meshUpdateEvent         = _nextMeshUpdateEvent;
  _nextMeshUpdateEvent     = MeshUpdateEvent::None;
}

void exahype::solvers::FiniteVolumesSolver::updateNextMeshUpdateEvent(
    exahype::solvers::Solver::MeshUpdateEvent meshUpdateEvent) {
  _nextMeshUpdateEvent = mergeMeshUpdateEvents(_nextMeshUpdateEvent,meshUpdateEvent);
}

exahype::solvers::Solver::MeshUpdateEvent
exahype::solvers::FiniteVolumesSolver::getMeshUpdateEvent() const {
  return _meshUpdateEvent;
}

void exahype::solvers::FiniteVolumesSolver::overwriteMeshUpdateEvent(MeshUpdateEvent newMeshUpdateEvent) {
   _meshUpdateEvent = newMeshUpdateEvent;
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
    const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
    const std::vector<std::string>& cmdlineargs,
    const exahype::parser::ParserView& parserView) {
  _domainOffset=domainOffset;
  _domainSize=domainSize;
  std::pair<double,int> coarsestMeshInfo =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(_maximumMeshSize,boundingBoxSize[0]);
  _coarsestMeshSize  = coarsestMeshInfo.first;
  _coarsestMeshLevel = coarsestMeshInfo.second;

  _previousMinTimeStepSize = 0.0;
  _previousMinTimeStamp = timeStamp;

  _minTimeStepSize = 0.0;
  _minTimeStamp = timeStamp;

  overwriteMeshUpdateEvent(MeshUpdateEvent::InitialRefinementRequested);

  init(cmdlineargs,parserView); // call user defined initalisation
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
    const int element) const {
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

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
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

     _maxLevel     = _nextMaxLevel;
     _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
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

  _maxLevel     = _nextMaxLevel;
  _nextMaxLevel = -std::numeric_limits<int>::max(); // "-", min
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

bool exahype::solvers::FiniteVolumesSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) {
  bool result = cellDescriptionsIndex>=0;
  assertion1(!result || Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  return result;
}

int exahype::solvers::FiniteVolumesSolver::tryGetElement(
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
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

int exahype::solvers::FiniteVolumesSolver::computeWeight(const int cellDescriptionsIndex) {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    return getCellDescriptions(cellDescriptionsIndex).size();
  }
  else return 0;
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInEnterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber,
    const bool stillInRefiningMode) {
  // Fine grid cell based uniform mesh refinement.
  const int fineGridCellElement =
      tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
  if (
      fineGridCellElement==exahype::solvers::Solver::NotFound &&
      tarch::la::allSmallerEquals(fineGridVerticesEnumerator.getCellSize(),getMaximumMeshSize()) &&
      tarch::la::allGreater(coarseGridVerticesEnumerator.getCellSize(),getMaximumMeshSize())
  ) {
    addNewCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
               multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
               solverNumber);
    return true;
    // Fine grid cell based adaptive mesh refinement operations are not implemented.
  } else {
    return false;
  }
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


void exahype::solvers::FiniteVolumesSolver::ensureNoUnnecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
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

void exahype::solvers::FiniteVolumesSolver::checkDataHeapIndex(const CellDescription& cellDescription, const int arrayIndex,const std::string arrayName) {
  assertion1(DataHeap::getInstance().isValidIndex(arrayIndex),cellDescription.toString());
  if ( arrayIndex < 0 ) {
    logError("checkDataHeapIndex(...)","The data heap array 'cellDescription."<<arrayName<<"' could not be allocated! Potential reason: Not enough memory available." <<
             " CellDescription="<<cellDescription.toString());
    std::abort();
  }
}

void exahype::solvers::FiniteVolumesSolver::ensureNecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  switch (cellDescription.getType()) {
    case CellDescription::Cell:
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolution())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
        // Allocate volume data
        const int patchSize = getDataPerPatch()+getGhostDataPerPatch();

        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setSolution(        DataHeap::getInstance().createData( patchSize, patchSize ));
          cellDescription.setPreviousSolution(DataHeap::getInstance().createData( patchSize, patchSize ));
          checkDataHeapIndex(cellDescription,cellDescription.getSolution(),"getSolution()");
          checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolution(),"getPreviousSolution()");

          cellDescription.setSolutionCompressed(-1);
          cellDescription.setPreviousSolutionCompressed(-1);

          cellDescription.setSolutionAverages(
              DataHeap::getInstance().createData( getNumberOfVariables()+getNumberOfParameters(), getNumberOfVariables()+getNumberOfParameters() ) );
          cellDescription.setPreviousSolutionAverages(
              DataHeap::getInstance().createData( getNumberOfVariables()+getNumberOfParameters(), getNumberOfVariables()+getNumberOfParameters() ) );
          checkDataHeapIndex(cellDescription,cellDescription.getSolutionAverages(),"getSolutionAverages()");
          checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionAverages(),"getPreviousSolutionAverages()");

          // Zero out the solution and previous solution arrays. For our MUSCL-Hancock implementation which
          // does not take the corner neighbours into account e.g., it is important that the values in
          // the corner cells of the first ghost layer are set to zero.
          std::fill_n( DataHeap::getInstance().getData(cellDescription.getSolution()).data(),         patchSize, 0.0 );
          std::fill_n( DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data(), patchSize, 0.0 );

          // Allocate boundary data
          const int patchBoundarySize = getDataPerPatchBoundary();
          cellDescription.setExtrapolatedSolution(DataHeap::getInstance().createData( patchBoundarySize, patchBoundarySize ));
          std::fill_n( DataHeap::getInstance().getData(cellDescription.getExtrapolatedSolution()).data(), patchBoundarySize, 0.0 );
          checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedSolution(),"getExtrapolatedSolution()");

          cellDescription.setExtrapolatedSolutionCompressed(-1);
          cellDescription.setExtrapolatedSolutionAverages( DataHeap::getInstance().createData(
            (getNumberOfVariables()+getNumberOfParameters()) * 2 * DIMENSIONS,
            (getNumberOfVariables()+getNumberOfParameters()) * 2 * DIMENSIONS ) );
          checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedSolutionAverages(),"getExtrapolatedSolutionAverages()");
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

bool exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInLeaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
     const int solverNumber) {
  return false;
}

exahype::solvers::Solver::RefinementControl
exahype::solvers::FiniteVolumesSolver::eraseOrRefineAdjacentVertices(
        const int cellDescriptionsIndex,
        const int solverNumber,
        const tarch::la::Vector<DIMENSIONS, double>& cellOffset,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize,
        const bool checkThoroughly) const {
  if ( tarch::la::oneGreater(cellSize,_maximumMeshSize) ) {
    return RefinementControl::Refine;
  } else {
    return RefinementControl::Erase;
  }
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

double exahype::solvers::FiniteVolumesSolver::startNewTimeStep(CellDescription& cellDescription) {
  assertion1(cellDescription.getType()==exahype::records::FiniteVolumesCellDescription::Cell,cellDescription.toString());
  //         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo
  double* solution = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();

  double admissibleTimeStepSize = stableTimeStepSize(solution, cellDescription.getSize());
  assertion(!std::isnan(admissibleTimeStepSize));

  // n-1
  cellDescription.setPreviousTimeStamp(cellDescription.getTimeStamp());
  cellDescription.setPreviousTimeStepSize(cellDescription.getTimeStepSize());

  // n
  cellDescription.setTimeStamp(cellDescription.getTimeStamp()+cellDescription.getTimeStepSize());
  cellDescription.setTimeStepSize(admissibleTimeStepSize);

  return admissibleTimeStepSize;
}

double exahype::solvers::FiniteVolumesSolver::startNewTimeStepFused(
    CellDescription& cellDescription,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch) {
  assertion1(cellDescription.getType()==exahype::records::FiniteVolumesCellDescription::Cell,cellDescription.toString());
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
    CellDescription& cellDescription) const {
  cellDescription.setTimeStepSize(0.0);
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStep(CellDescription& cellDescription) const {
  cellDescription.setTimeStamp(cellDescription.getPreviousTimeStamp());
  cellDescription.setTimeStepSize(cellDescription.getPreviousTimeStepSize());

  cellDescription.setPreviousTimeStepSize(std::numeric_limits<double>::max());
}

void exahype::solvers::FiniteVolumesSolver::rollbackToPreviousTimeStepFused(CellDescription& cellDescription) const {
  rollbackToPreviousTimeStep(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::adjustSolutionDuringMeshRefinementBody(
    CellDescription& cellDescription,
    const bool isInitialMeshRefinement) {
  assertion(cellDescription.getType()==CellDescription::Cell);

  zeroTimeStepSizes(cellDescription); // TODO(Dominic): Still necessary?
  synchroniseTimeStepping(cellDescription);

  adjustSolution(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::adjustSolution(CellDescription& cellDescription) {
  double* solution = exahype::DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp(),
      cellDescription.getTimeStepSize());

  #ifdef Asserts
  for (int i=0; i<getDataPerPatch()+getGhostDataPerPatch(); i++) {
    assertion3(std::isfinite(solution[i]),cellDescription.toString(),"setInitialConditions(...)",i);
  }
  #endif
}

exahype::solvers::Solver::UpdateResult exahype::solvers::FiniteVolumesSolver::fusedTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isAtRemoteBoundary) {
  bool isSkeletonCell = isAtRemoteBoundary;

  if (
      !SpawnPredictionAsBackgroundJob ||
      isFirstIterationOfBatch ||
      isLastIterationOfBatch
  ) {
    CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

    updateSolution(cellDescription,cellDescriptionsIndex,isFirstIterationOfBatch);
    UpdateResult result;
    result._timeStepSize = startNewTimeStepFused(
        cellDescription,isFirstIterationOfBatch,isLastIterationOfBatch);
    return result;
  } else {
    FusedTimeStepJob fusedTimeStepJob( *this, cellDescriptionsIndex, element, isSkeletonCell );
    Solver::submitPredictionJob(fusedTimeStepJob,isSkeletonCell);
    return UpdateResult();
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::FiniteVolumesSolver::update(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary){
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  uncompress(cellDescription);

  updateSolution(cellDescription,cellDescriptionsIndex,true);
  UpdateResult result;
  result._timeStepSize = startNewTimeStep(cellDescription);

  compress(cellDescription,isAtRemoteBoundary);
  return result;
}

void exahype::solvers::FiniteVolumesSolver::compress(
    const int cellDescriptionsIndex,
    const int element,
    const bool isAtRemoteBoundary) const {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  compress(cellDescription,isAtRemoteBoundary);
}

void exahype::solvers::FiniteVolumesSolver::adjustSolutionDuringMeshRefinement(
    const int cellDescriptionsIndex,
    const int element) {
  const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
    AdjustSolutionDuringMeshRefinementJob job(*this,cellDescription,isInitialMeshRefinement);
    peano::datatraversal::TaskSet spawnedSet( job, peano::datatraversal::TaskSet::TaskType::Background  );
  } else {
    adjustSolutionDuringMeshRefinementBody(cellDescription,isInitialMeshRefinement);
  }
}

void exahype::solvers::FiniteVolumesSolver::updateSolution(
    CellDescription& cellDescription,
    const int cellDescriptionsIndex,
    const bool backupPreviousSolution) {
  double* newSolution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();
  double* solution    = DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data();
  if (backupPreviousSolution) {
    std::copy(newSolution,newSolution+getDataPerPatch()+getGhostDataPerPatch(),solution); // Copy (current solution) in old solution field.
  }

  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution[pre]");

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

  if ( !tarch::la::equals(cellDescription.getTimeStepSize(), 0.0) && tarch::la::smaller(admissibleTimeStepSize,cellDescription.getTimeStepSize()) ) {
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
  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution[post]");
}

void exahype::solvers::FiniteVolumesSolver::swapSolutionAndPreviousSolution(
    CellDescription& cellDescription) const {
  // Simply swap the heap indices
  const int previousSolution = cellDescription.getPreviousSolution();
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolution(previousSolution);
}


void exahype::solvers::FiniteVolumesSolver::prolongateFaceData(
    const int cellDescriptionsIndex,
    const int element) {
  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==CellDescription::Cell)
    return;

  assertionMsg(false,"Please implement!");
}

void exahype::solvers::FiniteVolumesSolver::restriction(
      const int cellDescriptionsIndex,
      const int element) {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::rollbackSolutionGlobally(
    const int cellDescriptionsIndex, const int solverElement,
    const bool fusedTimeStepping) const {
  // do nothing
  logError("rollbackSolutionGlobally(...)","Should have never been called");
  std::abort();
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

  synchroniseTimeStepping(cellDescription1);
  synchroniseTimeStepping(cellDescription2);

//  if (cellDescriptionsIndex1==2516 ||
//      cellDescriptionsIndex2==2516) {
//    std::cout << "cell1: "<< cellDescriptionsIndex1 << "," << cellDescription1.toString() << std::endl;
//    std::cout << "cell2: "<< cellDescriptionsIndex2 << "," << cellDescription2.toString() << std::endl;
//  }

  if (cellDescription1.getType()==CellDescription::Cell ||
      cellDescription2.getType()==CellDescription::Cell) {

    assertion1(cellDescription1.getTimeStamp()<std::numeric_limits<double>::max(),cellDescription1.toString());
    assertion1(cellDescription1.getTimeStepSize()<std::numeric_limits<double>::max(),cellDescription1.toString());
    assertion1(cellDescription2.getTimeStamp()<std::numeric_limits<double>::max(),cellDescription2.toString());
    assertion1(cellDescription2.getTimeStepSize()<std::numeric_limits<double>::max(),cellDescription2.toString());

    if ( CompressionAccuracy > 0.0 ) {
      peano::datatraversal::TaskSet uncompression(
        [&] () -> bool {
          uncompress(cellDescription1);
          return false;
        },
        [&] () -> bool {
          uncompress(cellDescription2);
          return false;
        },
        peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
        peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
        true
      );
    }

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

  synchroniseTimeStepping(cellDescription);

  if (cellDescription.getType()==CellDescription::Cell) {
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
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerForkOrJoinCommunication   = 2;
const int exahype::solvers::FiniteVolumesSolver::DataMessagesPerMasterWorkerCommunication = 0;

void exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  if ( isValidCellDescriptionIndex(cellDescriptionsIndex) ) {
    Heap::getInstance().sendData(cellDescriptionsIndex,toRank,x,level,messageType);
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

void exahype::solvers::FiniteVolumesSolver::ensureOnlyNecessaryMemoryIsAllocated(CellDescription& cellDescription) {
  auto* solver = RegisteredSolvers[cellDescription.getSolverNumber()];
  switch (solver->getType()) {
  case exahype::solvers::Solver::Type::ADERDG:
    assertionMsg(false,"Solver type not supported!");
    break;
  case exahype::solvers::Solver::Type::LimitingADERDG:
    static_cast<LimitingADERDGSolver*>(solver)->
        getLimiter()->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<LimitingADERDGSolver*>(solver)->
        getLimiter()->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  case exahype::solvers::Solver::Type::FiniteVolumes:
    static_cast<FiniteVolumesSolver*>(solver)->ensureNoUnnecessaryMemoryIsAllocated(cellDescription);
    static_cast<FiniteVolumesSolver*>(solver)->ensureNecessaryMemoryIsAllocated(cellDescription);
    break;
  }
}

void exahype::solvers::FiniteVolumesSolver::receiveCellDescriptions(
    const int                                     fromRank,
    exahype::Cell&                                localCell,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  Heap::getInstance().receiveData(
      localCell.getCellDescriptionsIndex(),fromRank,x,level,messageType);

  logDebug("mergeCellDescriptionsWithRemoteData(...)","received " <<
          Heap::getInstance().getData(localCell.getCellDescriptionsIndex()).size() <<
          " cell descriptions for cell (centre="<< x.toString() << "level="<< level << ")");

  for (auto& cellDescription : Heap::getInstance().getData(localCell.getCellDescriptionsIndex())) {
    resetIndicesAndFlagsOfReceivedCellDescription(
        cellDescription,multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex);
  }
}

void exahype::solvers::FiniteVolumesSolver::resetIndicesAndFlagsOfReceivedCellDescription(
    CellDescription& cellDescription,const int parentIndex) {
  cellDescription.setParentIndex(parentIndex);

  // Default field data indices
  cellDescription.setSolution(-1);
  cellDescription.setPreviousSolution(-1);
  cellDescription.setExtrapolatedSolution(-1);

  // compression
  cellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  cellDescription.setExtrapolatedSolutionCompressed(-1);
  cellDescription.setSolutionCompressed(-1);
  cellDescription.setPreviousSolutionCompressed(-1);

  cellDescription.setSolutionAverages(-1);
  cellDescription.setSolutionAverages(-1);
  cellDescription.setExtrapolatedSolutionAverages(-1);
  cellDescription.setBytesPerDoFInPreviousSolution(-1);
  cellDescription.setBytesPerDoFInSolution(-1);
  cellDescription.setBytesPerDoFInExtrapolatedSolution(-1);
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

void
exahype::solvers::FiniteVolumesSolver::appendMasterWorkerCommunicationMetadata(
    exahype::MetadataHeap::HeapEntries& metadata,
    const int cellDescriptionsIndex,
    const int solverNumber) const {
  for (int i = 0; i < exahype::MasterWorkerCommunicationMetadataPerSolver; ++i) {
    metadata.push_back(exahype::InvalidMetadataEntry); // implicit conversion
  }
}

///////////////////////

void exahype::solvers::FiniteVolumesSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);
  assertion1(element>=0,element);
  assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),
      element,Heap::getInstance().getData(cellDescriptionsIndex).size());

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  logDebug("sendDataToWorkerOrMasterDueToForkOrJoin(...)","solution of solver " << cellDescription.getSolverNumber() << " sent to rank "<<toRank<<
      ", cell: "<< x << ", level: " << level);

  assertion(cellDescription.getType()==CellDescription::Cell);
  assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()),
      cellDescriptionsIndex,cellDescription.toString());
  assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()),
      cellDescriptionsIndex,cellDescription.toString());
  DataHeap::getInstance().sendData(
      DataHeap::getInstance().getData(cellDescription.getSolution()).data(),
      getDataPerPatch()+getGhostDataPerPatch(), toRank, x, level, messageType);
  DataHeap::getInstance().sendData(
      DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).data(),
      getDataPerPatch()+getGhostDataPerPatch(), toRank, x, level, messageType);
}

void exahype::solvers::FiniteVolumesSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const peano::heap::MessageType&               messageType,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  auto& cellDescription = getCellDescription(cellDescriptionsIndex,element);
  #ifdef Asserts
  const tarch::la::Vector<DIMENSIONS,double> center = cellDescription.getOffset()+0.5*cellDescription.getSize();
  #endif
  assertion5(Vertex::equalUpToRelativeTolerance(x,center),x,center,level,cellDescription.getLevel(),tarch::parallel::Node::getInstance().getRank());
  assertion2(cellDescription.getLevel()==level,cellDescription.getLevel(),level);
  assertion(cellDescription.getType()==CellDescription::Cell);

  logDebug("mergeWithRemoteDataDueToForkOrJoin(...)","[solution] receive from rank "<<fromRank<<
           ", cell: "<< x << ", level: " << level);

  // allocate memory
  ensureNecessaryMemoryIsAllocated(cellDescription);
  ensureNoUnnecessaryMemoryIsAllocated(cellDescription);

  DataHeap::getInstance().getData(cellDescription.getSolution()).clear();
  DataHeap::getInstance().receiveData( cellDescription.getSolution(),
      fromRank, x, level, messageType);
  DataHeap::getInstance().getData(cellDescription.getPreviousSolution()).clear();
  DataHeap::getInstance().receiveData( cellDescription.getPreviousSolution(),
      fromRank, x, level, messageType );
}

void exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInPrepareSendToWorker(
    const int workerRank,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const int solverNumber) {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::sendDataToWorkerIfProlongating(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // do nothing
}

void exahype::solvers::FiniteVolumesSolver::receiveDataFromMasterIfProlongating(
    const int masterRank,
    const int receivedCellDescriptionsIndex,
    const int receivedElement,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
  // do nothing
}

bool exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInMergeWithWorker(
    const int localCellDescriptionsIndex,
    const int receivedCellDescriptionsIndex, const int receivedElement) {
  return false;
}

void exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInPrepareSendToMaster(
    const int masterRank,
    const int cellDescriptionsIndex, const int element,
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const int level) const {
 // do nothing
}

bool exahype::solvers::FiniteVolumesSolver::progressMeshRefinementInMergeWithMaster(
    const int worker,
    const int localCellDescriptionsIndex,
    const int localElement,
    const int coarseGridCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level,
    const bool stillInRefiningMode) {
  // do nothing
  return false;
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

  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolution()));
  assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolution()));

  const int numberOfFaceDof = getDataPerPatchFace();
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
}

void exahype::solvers::FiniteVolumesSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  // Send order: lQhbnd,lFhbnd,observablesMin,observablesMax
  // Receive order: observablesMax,observablesMin,lFhbnd,lQhbnd
  // TODO(WORKAROUND)
  #if defined(UsePeanosSymmetricBoundaryExchanger)
  DataHeap::getInstance().sendData(
      _invalidExtrapolatedSolution, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
  #else
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        exahype::EmptyDataHeapMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
  #endif
}

void exahype::solvers::FiniteVolumesSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertionEquals(tarch::la::countEqualEntries(src,dest),DIMENSIONS-1); // We only consider faces; no corners.

  CellDescription& cellDescription = getCellDescription(cellDescriptionsIndex,element);

  synchroniseTimeStepping(cellDescription);

  const int direction   = tarch::la::equalsReturnIndex(src, dest);
  const int orientation = (1 + src(direction) - dest(direction))/2;
  const int faceIndex   = 2*direction+orientation;

  assertion3(cellDescription.getNeighbourMergePerformed(faceIndex),
     faceIndex,cellDescriptionsIndex,cellDescription.toString());
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
  DataHeap::HeapEntries timeStepDataToReduce(0,1);
  timeStepDataToReduce.push_back(_minTimeStepSize);

  assertion1(timeStepDataToReduce.size()==1,timeStepDataToReduce.size());
  assertion1(std::isfinite(timeStepDataToReduce[0]),timeStepDataToReduce[0]);
  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
        tarch::parallel::Node::getInstance().getRank());
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending time step data: " <<
        "data[0]=" << timeStepDataToReduce[0]);
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
  DataHeap::HeapEntries receivedTimeStepData(1);

  if (true || tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data [pre] from rank " << workerRank);
  }

  DataHeap::getInstance().receiveData(
      receivedTimeStepData.data(),receivedTimeStepData.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);

  assertion1(receivedTimeStepData.size()==1,receivedTimeStepData.size());
  assertion1(receivedTimeStepData[0]>=0,receivedTimeStepData[0]);
  assertion1(std::isfinite(receivedTimeStepData[0]),receivedTimeStepData[0]);
  // The master solver has not yet updated its minNextTimeStepSize.
  // Thus it does not yet equal MAX_DOUBLE.

  int index=0;
  _minNextTimeStepSize = std::min( _minNextTimeStepSize, receivedTimeStepData[index++] );

  if (true || tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Receiving time step data: " <<
             "data[0]=" << receivedTimeStepData[0]);

    logDebug("mergeWithWorkerData(...)","Updated time step fields: " <<
             "_minNextTimeStepSize="     << _minNextTimeStepSize);
  }
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::FiniteVolumesSolver::sendDataToWorker(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  std::vector<double> timeStepDataToSend(0,2);
  timeStepDataToSend.push_back(_minTimeStamp); // TODO(Dominic): Append previous time step size
  timeStepDataToSend.push_back(_minTimeStepSize);

  assertion1(timeStepDataToSend.size()==2,timeStepDataToSend.size());
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
        ",data[1]=" << timeStepDataToSend[1]);
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
  std::vector<double> receivedTimeStepData(2);
  DataHeap::getInstance().receiveData(
      receivedTimeStepData.data(),receivedTimeStepData.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion1(receivedTimeStepData.size()==2,receivedTimeStepData.size());

  if (_timeStepping==TimeStepping::Global) {
    assertionNumericalEquals1(_minNextTimeStepSize,std::numeric_limits<double>::max(),
        _minNextTimeStepSize);
  }

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithMasterData(...)","Received time step data: " <<
        "data[0]="  << receivedTimeStepData[0] <<
        ",data[1]=" << receivedTimeStepData[1]);
  }

  _minTimeStamp    = receivedTimeStepData[0];
  _minTimeStepSize = receivedTimeStepData[1];
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
    [&]() -> bool  { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getPreviousSolution() ).data(),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&] () -> bool  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getSolution() ).data(),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfExtrapolatedSolution = peano::heap::findMostAgressiveCompression(
      DataHeap::getInstance().getData( cellDescription.getExtrapolatedSolution() ).data(),
      getDataPerPatchBoundary(),
      CompressionAccuracy,true
      );
      return false;
      },
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
    [&]() -> bool {
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
      return false;
    },
    [&]() -> bool {
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
      return false;
    },
    [&]() -> bool {
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
      return false;
    },
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
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
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setPreviousSolution( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolution()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setSolution( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedSolution()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
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
    [&]() -> bool {
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
      return false;
    },
    [&]() -> bool {
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
      return false;
    },
    [&]() -> bool {
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
      return false;
    },
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	peano::datatraversal::TaskSet::TaskType::IsTaskAndRunAsSoonAsPossible,
	true
  );
}


void exahype::solvers::FiniteVolumesSolver::compress(CellDescription& cellDescription,const bool isSkeletonCell) const {
  assertion1( cellDescription.getCompressionState() ==  CellDescription::Uncompressed, cellDescription.toString() );
  if (CompressionAccuracy>0.0) {
    if ( SpawnCompressionAsBackgroundJob ) {
      cellDescription.setCompressionState(CellDescription::CurrentlyProcessed);

      int& jobCounter = (isSkeletonCell) ? NumberOfSkeletonJobs: NumberOfEnclaveJobs;
      CompressionJob compressionJob( *this, cellDescription, jobCounter );
      peano::datatraversal::TaskSet spawnedSet( compressionJob,peano::datatraversal::TaskSet::TaskType::Background );
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

    peano::datatraversal::TaskSet::finishToProcessBackgroundJobs();
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

exahype::solvers::FiniteVolumesSolver::CompressionJob::CompressionJob(
  const FiniteVolumesSolver& solver,
  CellDescription&           cellDescription,
  const bool                 isSkeletonJob)
  :
  _solver(solver),
  _cellDescription(cellDescription),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}


bool exahype::solvers::FiniteVolumesSolver::CompressionJob::operator()() {
  _solver.determineUnknownAverages(_cellDescription);
  _solver.computeHierarchicalTransform(_cellDescription,-1.0);
  _solver.putUnknownsIntoByteStream(_cellDescription);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    _cellDescription.setCompressionState(CellDescription::Compressed);
    // @todo raus (TODO(Dominic): why?)
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}


exahype::solvers::FiniteVolumesSolver::FusedTimeStepJob::FusedTimeStepJob(
  FiniteVolumesSolver&     solver,
  const int                cellDescriptionsIndex,
  const int                element,
  const bool               isSkeletonJob):
  _solver(solver),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _element(element),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::FiniteVolumesSolver::FusedTimeStepJob::operator()() {
  _solver.fusedTimeStep(
      _cellDescriptionsIndex,_element,false,false,true);
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter--;
    assertion( jobCounter>=0 );
  }
  lock.free();
  return false;
}



exahype::solvers::FiniteVolumesSolver::AdjustSolutionDuringMeshRefinementJob::AdjustSolutionDuringMeshRefinementJob(
  FiniteVolumesSolver& solver,
  CellDescription&     cellDescription,
  const bool           isInitialMeshRefinement):
  _solver(solver),
  _cellDescription(cellDescription),
  _isInitialMeshRefinement(isInitialMeshRefinement)
{
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs++;
  }
  lock.free();
}

bool exahype::solvers::FiniteVolumesSolver::AdjustSolutionDuringMeshRefinementJob::operator()() {
  _solver.adjustSolutionDuringMeshRefinementBody(_cellDescription,_isInitialMeshRefinement);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
