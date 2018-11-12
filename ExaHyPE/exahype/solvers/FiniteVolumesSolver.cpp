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

void exahype::solvers::FiniteVolumesSolver::synchroniseTimeStepping(CellDescription& cellDescription) const {
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
  CellInfo cellInfo =
      fineGridCell.addNewCellDescription(
              solverNumber,
              CellDescription::Cell,
              CellDescription::None,
              fineGridVerticesEnumerator.getLevel(),
              coarseGridCellDescriptionsIndex,
              fineGridVerticesEnumerator.getCellSize(),
              fineGridVerticesEnumerator.getVertexPosition());

  const int fineGridCellElement = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[fineGridCellElement]; //TODO(Dominic): Multi-solvers: Might need to lock this?
  ensureNecessaryMemoryIsAllocated(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
    const int solverNumber,
    CellInfo& cellInfo,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
    const exahype::records::FiniteVolumesCellDescription::RefinementEvent refinementEvent,
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
    const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  assertion2(static_cast<unsigned int>(solverNumber) < solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

  CellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Background job completion monitoring (must be initialised with true)
  newCellDescription.setHasCompletedTimeStep(true);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setLevel(level);
  newCellDescription.setRefinementEvent(refinementEvent);
  // newCellDescription.setHelperCellNeedsToStoreFaceData(false); // TODO(Dominic): Add to FV cell descr.

  newCellDescription.setNeighbourMergePerformed((signed char) 0/*implicit conversion*/);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Default data field indices
  newCellDescription.setSolutionIndex(-1);
  newCellDescription.setSolution(nullptr);
  newCellDescription.setSolutionAveragesIndex(-1);
  newCellDescription.setSolutionAverages(nullptr);
  newCellDescription.setSolutionCompressedIndex(-1);
  newCellDescription.setSolutionCompressed(nullptr);
  newCellDescription.setPreviousSolutionIndex(-1);
  newCellDescription.setPreviousSolution(nullptr);
  newCellDescription.setPreviousSolutionAveragesIndex(-1);
  newCellDescription.setPreviousSolutionAverages(nullptr);
  newCellDescription.setPreviousSolutionCompressedIndex(-1);
  newCellDescription.setPreviousSolutionCompressed(nullptr);
  newCellDescription.setExtrapolatedSolutionIndex(-1);
  newCellDescription.setExtrapolatedSolution(nullptr);
  newCellDescription.setExtrapolatedSolutionAveragesIndex(-1);
  newCellDescription.setExtrapolatedSolutionAverages(nullptr);
  newCellDescription.setExtrapolatedSolutionCompressedIndex(-1);
  newCellDescription.setExtrapolatedSolutionCompressed(nullptr);

  newCellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
  cellInfo._FiniteVolumesCellDescriptions.push_back(newCellDescription);
  lock.free();
}


void exahype::solvers::FiniteVolumesSolver::ensureNoUnnecessaryMemoryIsAllocated(
    CellDescription& cellDescription) const {
  if (DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())) {
    switch (cellDescription.getType()) {
      case CellDescription::Erased: {
        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));
          assertion(DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionIndex()));

          if (cellDescription.getSolutionIndex()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getSolutionIndex());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getSolutionIndex()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getSolutionCompressedIndex());
          }
          DataHeap::getInstance().deleteData(cellDescription.getSolutionAveragesIndex());

          if (cellDescription.getPreviousSolutionIndex()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionIndex());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getPreviousSolutionIndex()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionCompressedIndex());
          }
          DataHeap::getInstance().deleteData(cellDescription.getPreviousSolutionAveragesIndex());

          if (cellDescription.getExtrapolatedSolutionIndex()>=0) {
            DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionIndex());
          }
          else {
            assertion(CompressionAccuracy>0.0);
            assertion(cellDescription.getExtrapolatedSolutionIndex()==-1);
            CompressedDataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionCompressedIndex());
          }
          DataHeap::getInstance().deleteData(cellDescription.getExtrapolatedSolutionAveragesIndex());

          cellDescription.setSolutionIndex(-1);
          cellDescription.setSolution(nullptr);
          cellDescription.setPreviousSolutionIndex(-1);
          cellDescription.setPreviousSolution(nullptr);
          cellDescription.setExtrapolatedSolutionIndex(-1);
          cellDescription.setExtrapolatedSolution(nullptr);

          cellDescription.setSolutionCompressedIndex(-1);
          cellDescription.setSolutionCompressed(nullptr);
          cellDescription.setPreviousSolutionCompressedIndex(-1);
          cellDescription.setPreviousSolutionCompressed(nullptr);
          cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
          cellDescription.setExtrapolatedSolutionCompressed(nullptr);

          cellDescription.setSolutionAveragesIndex(-1);
          cellDescription.setSolutionAverages(nullptr);
          cellDescription.setPreviousSolutionAveragesIndex(-1);
          cellDescription.setPreviousSolutionAverages(nullptr);
          cellDescription.setExtrapolatedSolutionAveragesIndex(-1);
          cellDescription.setExtrapolatedSolutionAverages(nullptr);
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
      if (!DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex())) {
        assertion(!DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
        // Allocate volume data
        const int patchSize = getDataPerPatch()+getGhostDataPerPatch();

        tarch::multicore::Lock lock(exahype::HeapSemaphore);
          cellDescription.setSolutionIndex(        DataHeap::getInstance().createData( patchSize, patchSize ));
          cellDescription.setPreviousSolutionIndex(DataHeap::getInstance().createData( patchSize, patchSize ));
          checkDataHeapIndex(cellDescription,cellDescription.getSolutionIndex(),"getSolutionIndex()");
          checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionIndex(),"getPreviousSolutionIndex()");
          cellDescription.setSolution        ( getDataHeapEntries(cellDescription.getSolutionIndex()).data() ) ;
          cellDescription.setPreviousSolution( getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).data() ) ;
          // std::fill_n(static_cast<double*>(cellDescription.getPreviousSolution()),patchSize,std::numeric_limits<double>::quiet_NaN());
          // std::fill_n(static_cast<double*>(cellDescription.getSolution()),patchSize,std::numeric_limits<double>::quiet_NaN());
          // Zero out the solution and previous solution arrays. For our MUSCL-Hancock implementation which
          // does not take the corner neighbours into account e.g., it is important that the values in
          // the corner cells of the first ghost layer are set to zero.
          std::fill_n( static_cast<double*>(cellDescription.getSolution()),         patchSize, 0.0 );
          std::fill_n( static_cast<double*>(cellDescription.getPreviousSolution()), patchSize, 0.0 );

          cellDescription.setSolutionCompressedIndex(-1);
          cellDescription.setSolutionCompressed(nullptr);
          cellDescription.setPreviousSolutionCompressedIndex(-1);
          cellDescription.setPreviousSolutionCompressed(nullptr);

          const int dataPerSubcell = getNumberOfVariables()+getNumberOfParameters();
          cellDescription.setSolutionAveragesIndex        (DataHeap::getInstance().createData( dataPerSubcell, dataPerSubcell ) );
          cellDescription.setPreviousSolutionAveragesIndex(DataHeap::getInstance().createData( dataPerSubcell, dataPerSubcell ) );
          checkDataHeapIndex(cellDescription,cellDescription.getSolutionAveragesIndex(),"getSolutionAveragesIndex()");
          checkDataHeapIndex(cellDescription,cellDescription.getPreviousSolutionAveragesIndex(),"getPreviousSolutionAveragesIndex()");
          cellDescription.setSolutionAverages        ( getDataHeapEntries(cellDescription.getSolutionAveragesIndex()).data() ) ;
          cellDescription.setPreviousSolutionAverages( getDataHeapEntries(cellDescription.getPreviousSolutionAveragesIndex()).data() ) ;
          std::fill_n(static_cast<double*>(cellDescription.getPreviousSolutionAverages()),dataPerSubcell,std::numeric_limits<double>::quiet_NaN());
          std::fill_n(static_cast<double*>(cellDescription.getSolutionAverages()),dataPerSubcell,std::numeric_limits<double>::quiet_NaN());

          // Allocate boundary data
          const int patchBoundarySize = getDataPerPatchBoundary();
          cellDescription.setExtrapolatedSolutionIndex(DataHeap::getInstance().createData( patchBoundarySize, patchBoundarySize ));
          checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedSolutionIndex(),"getExtrapolatedSolutionIndex()");
          cellDescription.setExtrapolatedSolution( getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).data() ) ;
          std::fill_n( static_cast<double*>(cellDescription.getExtrapolatedSolution()), patchBoundarySize, 0.0 );

          cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
          cellDescription.setExtrapolatedSolutionCompressed(nullptr);

          cellDescription.setExtrapolatedSolutionAveragesIndex( DataHeap::getInstance().createData(dataPerSubcell*2*DIMENSIONS, dataPerSubcell*2*DIMENSIONS ) );
          checkDataHeapIndex(cellDescription,cellDescription.getExtrapolatedSolutionAveragesIndex(),"getExtrapolatedSolutionAveragesIndex()");
          cellDescription.setExtrapolatedSolutionAverages( getDataHeapEntries(cellDescription.getExtrapolatedSolutionAveragesIndex()).data() ) ;
          std::fill_n(static_cast<double*>(cellDescription.getExtrapolatedSolutionAverages()),dataPerSubcell*2*DIMENSIONS,std::numeric_limits<double>::quiet_NaN());
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
     const int solverNumber,
     const bool stillInRefiningMode) {
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
      const int solverNumber,
      CellInfo& cellInfo) {
  // do nothing
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////

double exahype::solvers::FiniteVolumesSolver::startNewTimeStep(CellDescription& cellDescription) {
  return startNewTimeStepFused(cellDescription,true,true);
}

double exahype::solvers::FiniteVolumesSolver::startNewTimeStepFused(
    CellDescription& cellDescription,
    const bool isFirstIterationOfBatch, // TODOD(Dominic): same code
    const bool isLastIterationOfBatch) {
  assertion1(cellDescription.getType()==exahype::records::FiniteVolumesCellDescription::Cell,cellDescription.toString());
  double* solution = static_cast<double*>(cellDescription.getSolution());

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

double exahype::solvers::FiniteVolumesSolver::updateTimeStepSizes(
    const int solverNumber,CellInfo& cellInfo,const bool fused) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    if ( cellDescription.getType()==exahype::records::FiniteVolumesCellDescription::Cell ) {
      double* solution = static_cast<double*>(cellDescription.getSolution());

      double admissibleTimeStepSize = stableTimeStepSize(solution, cellDescription.getSize());

      assertion(!std::isnan(admissibleTimeStepSize));
      cellDescription.setTimeStepSize(admissibleTimeStepSize);

      return admissibleTimeStepSize;
    } else {
      return std::numeric_limits<double>::max();
    }
  } else {
    return std::numeric_limits<double>::max();
  }
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

  synchroniseTimeStepping(cellDescription);

  adjustSolution(cellDescription);
}

void exahype::solvers::FiniteVolumesSolver::adjustSolution(CellDescription& cellDescription) {
  double* solution = static_cast<double*>(cellDescription.getSolution());
  adjustSolution(
      solution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getTimeStamp(),
      cellDescription.getTimeStepSize());

  double* previousSolution = static_cast<double*>(cellDescription.getPreviousSolution());
  adjustSolution(
      previousSolution,
      cellDescription.getOffset()+0.5*cellDescription.getSize(),
      cellDescription.getSize(),
      cellDescription.getPreviousTimeStamp(),
      cellDescription.getPreviousTimeStepSize());

  #ifdef Asserts
  for (int i=0; i<getDataPerPatch()+getGhostDataPerPatch(); i++) {
    assertion3(std::isfinite(solution[i]),cellDescription.toString(),"setInitialConditions(...)",i);
  }
  #endif
}

exahype::solvers::Solver::UpdateResult exahype::solvers::FiniteVolumesSolver::updateBody(
    CellDescription& cellDescription,
    const int cellDescriptionsIndex,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isAtRemoteBoundary,
    const bool uncompressBefore) {
  if ( uncompressBefore ) { uncompress(cellDescription); }

  updateSolution(cellDescription,cellDescriptionsIndex,isFirstIterationOfBatch);
  UpdateResult result;
  result._timeStepSize = startNewTimeStepFused(cellDescription,isFirstIterationOfBatch,isLastIterationOfBatch);

  cellDescription.setHasCompletedTimeStep(true); // last step of the FV update

  compress(cellDescription,isAtRemoteBoundary);
  return result;
}

exahype::solvers::Solver::UpdateResult exahype::solvers::FiniteVolumesSolver::fusedTimeStepOrRestrict(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isFirstIterationOfBatch,
    const bool isLastIterationOfBatch,
    const bool isAtRemoteBoundary) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    bool isSkeletonCell = isAtRemoteBoundary;
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    if (
        !SpawnPredictionAsBackgroundJob ||
        isFirstIterationOfBatch ||
        isLastIterationOfBatch
    ) {
      return updateBody(
          cellDescription,cellInfo._cellDescriptionsIndex,
          isFirstIterationOfBatch,isLastIterationOfBatch,isAtRemoteBoundary,false/*uncompressBefore*/);
    } else {
      cellDescription.setHasCompletedTimeStep(false);
      peano::datatraversal::TaskSet( new FusedTimeStepJob( *this, cellDescription, cellInfo._cellDescriptionsIndex, isSkeletonCell ) );
      return UpdateResult();
    }
  } else {
    return UpdateResult();
  }
}

exahype::solvers::Solver::UpdateResult exahype::solvers::FiniteVolumesSolver::updateOrRestrict(
      const int  solverNumber,
      CellInfo&  cellInfo,
      const bool isAtRemoteBoundary){
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element!=NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    return updateBody(cellDescription,cellInfo._cellDescriptionsIndex,true,true,isAtRemoteBoundary,true/*uncompressBefore*/);
  } else {
    return UpdateResult();
  }
}

void exahype::solvers::FiniteVolumesSolver::compress(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool isAtRemoteBoundary) const {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    compress(cellDescription,isAtRemoteBoundary);
  }
}

void exahype::solvers::FiniteVolumesSolver::adjustSolutionDuringMeshRefinement(
    const int solverNumber,CellInfo& cellInfo) {
  const bool isInitialMeshRefinement = getMeshUpdateEvent()==MeshUpdateEvent::InitialRefinementRequested;
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    if ( exahype::solvers::Solver::SpawnAMRBackgroundJobs ) {
      peano::datatraversal::TaskSet( new AdjustSolutionDuringMeshRefinementJob(*this,cellDescription,isInitialMeshRefinement) );
    } else {
      adjustSolutionDuringMeshRefinementBody(cellDescription,isInitialMeshRefinement);
    }
  }
}

void exahype::solvers::FiniteVolumesSolver::updateSolution(
    CellDescription& cellDescription,
    const int cellDescriptionsIndex,
    const bool backupPreviousSolution) {
  assertion1( tarch::la::equals(cellDescription.getNeighbourMergePerformed(),(signed char) true) || ProfileUpdate,cellDescription.toString());
  if ( !tarch::la::equals(cellDescription.getNeighbourMergePerformed(),(signed char) true) && !ProfileUpdate ) {
    logError("updateSolution(...)","Not all ghost layers were copied to cell="<<cellDescription.toString());
    std::terminate();
  }

  #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
    static int counter = 0;
    static double timeStamp = 0;
    if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
      logInfo("mergeNeighboursData(...)","#updateSolution="<<counter);
      timeStamp = _minTimeStamp;
      counter=0;
    }
    counter++;
  #endif

  double* newSolution = static_cast<double*>(cellDescription.getSolution());
  double* solution    = static_cast<double*>(cellDescription.getPreviousSolution());
  if (backupPreviousSolution) {
    std::copy(newSolution,newSolution+getDataPerPatch()+getGhostDataPerPatch(),solution); // Copy (current solution) in old solution field.
  }

  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution[pre]");

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

  // only for profiling
  if ( Solver::ProfileUpdate ) { swapSolutionAndPreviousSolution(cellDescription); }

  validateNoNansInFiniteVolumesSolution(cellDescription,cellDescriptionsIndex,"updateSolution[post]");
}

void exahype::solvers::FiniteVolumesSolver::swapSolutionAndPreviousSolution(
    CellDescription& cellDescription) const {
  // Simply swap the heap indices
  const int previousSolutionIndex = cellDescription.getPreviousSolutionIndex();
  void* previousSolution          = cellDescription.getPreviousSolution(); // pointer
  cellDescription.setPreviousSolutionIndex(cellDescription.getSolutionIndex());
  cellDescription.setPreviousSolution(cellDescription.getSolution());
  cellDescription.setSolutionIndex(previousSolutionIndex);
  cellDescription.setSolution(previousSolution);
}

void exahype::solvers::FiniteVolumesSolver::rollbackSolutionGlobally(
    const int solverNumber,
    CellInfo& cellInfo,
    const bool fusedTimeStepping) const {
  // do nothing
  logError("rollbackSolutionGlobally(...)","Not implemented");
  std::abort();
}

///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////

void exahype::solvers::FiniteVolumesSolver::mergeNeighboursData(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo1,
    Solver::CellInfo&                         cellInfo2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  const int element1 = cellInfo1.indexOfFiniteVolumesCellDescription(solverNumber);
  const int element2 = cellInfo2.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element1 != Solver::NotFound && element2 != Solver::NotFound ) {
    Solver::InterfaceInfo face(pos1,pos2);
    CellDescription& cellDescription1 = cellInfo1._FiniteVolumesCellDescriptions[element1];
    CellDescription& cellDescription2 = cellInfo2._FiniteVolumesCellDescriptions[element2];

    #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
    static int counter = 0;
    static double timeStamp = 0;
    if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
     logInfo("mergeNeighboursData(...)","#riemanns="<<counter);
     timeStamp = _minTimeStamp;
     counter=0;
    }
    counter++;
    #endif

    synchroniseTimeStepping(cellDescription1);
    synchroniseTimeStepping(cellDescription2);

    waitUntilCompletedTimeStep<CellDescription>(cellDescription1,false,false);
    waitUntilCompletedTimeStep<CellDescription>(cellDescription2,false,false);

    assertion(cellDescription1.getType()==CellDescription::Cell && cellDescription2.getType()==CellDescription::Cell);

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

    double* solution1 = static_cast<double*>(cellDescription1.getSolution());
    double* solution2 = static_cast<double*>(cellDescription2.getSolution());

    ghostLayerFilling(solution1,solution2,pos2-pos1);
    ghostLayerFilling(solution2,solution1,pos1-pos2);
  }
}

void exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData(
    const int                                 solverNumber,
    Solver::CellInfo&                         cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>& posCell,
    const tarch::la::Vector<DIMENSIONS, int>& posBoundary) {
  assertion2(tarch::la::countEqualEntries(posCell,posBoundary)==(DIMENSIONS-1),posCell.toString(),posBoundary.toString());
  Solver::BoundaryFaceInfo face(posCell,posBoundary);

  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    assertion1( cellDescription.getType()==CellDescription::Cell, cellDescription.toString() );

    #if !defined(SharedMemoryParallelisation) && !defined(Parallel) && defined(Asserts)
    static int counter = 0;
    static double timeStamp = 0;
    if ( !tarch::la::equals(timeStamp,_minTimeStamp,1e-9) ) {
      logInfo("applyBoundaryConditions(...)","#boundaryConditions="<<counter);
      timeStamp = _minTimeStamp;
      counter=0;
    }
    counter++;
    #endif

    waitUntilCompletedTimeStep<CellDescription>(cellDescription,false,false);

    uncompress(cellDescription);

    double* luh = static_cast<double*>(cellDescription.getSolution());
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
  cellDescription.setSolutionIndex(-1);
  cellDescription.setSolution(nullptr);
  cellDescription.setPreviousSolutionIndex(-1);
  cellDescription.setPreviousSolution(nullptr);
  cellDescription.setExtrapolatedSolutionIndex(-1);
  cellDescription.setExtrapolatedSolution(nullptr);

  // compression
  cellDescription.setCompressionState(CellDescription::CompressionState::Uncompressed);

  cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
  cellDescription.setExtrapolatedSolutionCompressed(nullptr);
  cellDescription.setSolutionCompressedIndex(-1);
  cellDescription.setSolutionCompressed(nullptr);
  cellDescription.setPreviousSolutionCompressedIndex(-1);
  cellDescription.setPreviousSolutionCompressed(nullptr);

  cellDescription.setSolutionAveragesIndex(-1);
  cellDescription.setSolutionAverages(nullptr);
  cellDescription.setExtrapolatedSolutionAveragesIndex(-1);
  cellDescription.setExtrapolatedSolutionAverages(nullptr);

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
  assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()),
      cellDescriptionsIndex,cellDescription.toString());
  assertion2(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()),
      cellDescriptionsIndex,cellDescription.toString());
  DataHeap::getInstance().sendData(
      static_cast<double*>(cellDescription.getSolution()),
      getDataPerPatch()+getGhostDataPerPatch(), toRank, x, level, messageType);
  DataHeap::getInstance().sendData(
      static_cast<double*>(cellDescription.getPreviousSolution()),
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

  getDataHeapEntries(cellDescription.getSolutionIndex()).clear();
  DataHeap::getInstance().receiveData( cellDescription.getSolutionIndex(),
      fromRank, x, level, messageType);
  getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).clear();
  DataHeap::getInstance().receiveData( cellDescription.getPreviousSolutionIndex(),
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

void exahype::solvers::FiniteVolumesSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     solverNumber,
    Solver::CellInfo&                             cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != Solver::NotFound ) {
    Solver::BoundaryFaceInfo face(src,dest);

    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    waitUntilCompletedTimeStep<CellDescription>(cellDescription,true,true);

    const int dataPerFace = getDataPerPatchFace();
    double* luhbnd    = static_cast<double*>(cellDescription.getExtrapolatedSolution()) + (face._faceIndex * dataPerFace);
    const double* luh = static_cast<double*>(cellDescription.getSolution());
    boundaryLayerExtraction(luhbnd,luh,dest-src);

    // Send order: minMax,lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd,minMax
    logDebug("sendDataToNeighbour(...)","send "<<DataMessagesPerNeighbourCommunication<<" arrays to rank=" <<toRank << ",cell="<<cellDescription.getOffset()<<",x="<<x<<",level="<<level);

    DataHeap::getInstance().sendData(
        luhbnd, dataPerFace, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);
  }
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
    const int                                    solverNumber,
    Solver::CellInfo&                            cellInfo,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);
  if ( element != NotFound ) {
    Solver::BoundaryFaceInfo face(dest,src);
    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];
    synchroniseTimeStepping(cellDescription);

    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()));
    assertion(DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()));

    logDebug("mergeWithNeighbourData(...)", "receive "<<DataMessagesPerNeighbourCommunication<<" arrays from rank="<<fromRank<<",x="<<x<<",level="<<level);

    // TODO(Dominic): If anarchic time stepping, receive the time step too.
    //
    // Copy the received boundary layer into a ghost layer of the solution.
    // TODO(Dominic): Pipe it directly through the Riemann solver if
    // we only use the Godunov method and not higher-order FVM methods.
    // For methods that are higher order in time, e.g., MUSCL-Hancock, we usually need
    // corner neighbours. This is why we currently adapt a GATHER-UPDATE algorithm
    // instead of a SOLVE RIEMANN PROBLEM AT BOUNDARY-UPDATE INTERIOR scheme.
    const int dataPerFace = getDataPerPatchFace();
    double* luhbnd = static_cast<double*>(cellDescription.getExtrapolatedSolution()) + (face._faceIndex * dataPerFace);
 
    // Send order: minMax,lQhbnd,lFhbnd
    // Receive order: lFhbnd,lQhbnd,minMax
    DataHeap::getInstance().receiveData(
        luhbnd, dataPerFace, fromRank, x, level,peano::heap::MessageType::NeighbourCommunication);

    double* luh = static_cast<double*>(cellDescription.getSolution());
    ghostLayerFillingAtBoundary(luh,luhbnd,src-dest);
  }
}

void exahype::solvers::FiniteVolumesSolver::dropNeighbourData(
    const int                                     fromRank,
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
  double* solution = static_cast<double*>(cellDescription.getSolution());

  dfor(i,_nodesPerCoordinateAxis+_ghostLayerWidth) {
    if (tarch::la::allSmaller(i,_nodesPerCoordinateAxis+_ghostLayerWidth)
    && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
      for (int unknown=0; unknown < _numberOfVariables; unknown++) {
        int iScalar = peano::utils::dLinearisedWithoutLookup(i,_nodesPerCoordinateAxis+2*_ghostLayerWidth)*_numberOfVariables+unknown;
        // cellDescription.getTimeStepSize()==0.0 is an initial condition
        assertion7(tarch::la::equals(cellDescription.getTimeStepSize(),0.0)  || std::isfinite(solution[iScalar]),
                   cellDescription.toString(),cellDescriptionsIndex,solution[iScalar],i.toString(),
                   _nodesPerCoordinateAxis,_ghostLayerWidth,
                   methodTrace);
      }
    }
  } // Dead code elimination should get rid of this loop if Asserts is not set.
  #endif
}

void exahype::solvers::FiniteVolumesSolver::printFiniteVolumesSolution(
    CellDescription& cellDescription)  const {
  #if DIMENSIONS==2
  double* solution = static_cast<double*>(cellDescription.getSolution());

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
  double* solution = static_cast<double*>(cellDescription.getSolution());

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

  assertion( cellDescription.getPreviousSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getSolutionCompressedIndex()==-1 );
  assertion( cellDescription.getExtrapolatedSolutionCompressedIndex()==-1 );

  int compressionOfPreviousSolution;
  int compressionOfSolution;
  int compressionOfExtrapolatedSolution;

  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ));

  peano::datatraversal::TaskSet compressionFactorIdentification(
    [&]() -> bool  { compressionOfPreviousSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getPreviousSolution()),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&] () -> bool  { compressionOfSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getSolution()),
      getDataPerPatch() + getGhostDataPerPatch(),
      CompressionAccuracy,true
      );
      return false;
      },
    [&]() -> bool  { compressionOfExtrapolatedSolution = peano::heap::findMostAgressiveCompression(
      static_cast<double*>(cellDescription.getExtrapolatedSolution()),
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
        cellDescription.setPreviousSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
        assertion1(
          cellDescription.getPreviousSolutionCompressedIndex()>=0,
          cellDescription.toString()
        );
        lock.free();
        cellDescription.setPreviousSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionCompressedIndex()).data() ) );

        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
        tearApart(numberOfEntries, cellDescription.getPreviousSolutionIndex(), cellDescription.getPreviousSolutionCompressedIndex(), compressionOfPreviousSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getPreviousSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionIndex(), true );
        cellDescription.setPreviousSolutionIndex(-1);
        cellDescription.setPreviousSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += getDataHeapEntries(cellDescription.getPreviousSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInSolution(compressionOfSolution);
      if (compressionOfSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        cellDescription.setSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
        assertion1( cellDescription.getSolutionCompressedIndex()>=0, cellDescription.getSolutionCompressedIndex() );
        lock.free();
        cellDescription.setSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionCompressedIndex()).data() ) );

        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();

        tearApart(numberOfEntries, cellDescription.getSolutionIndex(), cellDescription.getSolutionCompressedIndex(), compressionOfSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getSolutionIndex(), true );
        cellDescription.setSolutionIndex(-1);
        cellDescription.setSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += getDataHeapEntries(cellDescription.getSolutionIndex()).size() * 8.0;
        lock.free();
        #endif
      }
      return false;
    },
    [&]() -> bool {
      cellDescription.setBytesPerDoFInExtrapolatedSolution(compressionOfExtrapolatedSolution);
      if (compressionOfExtrapolatedSolution<7) {
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        cellDescription.setExtrapolatedSolutionCompressedIndex( CompressedDataHeap::getInstance().createData(0,0) );
        assertion( cellDescription.getExtrapolatedSolutionCompressedIndex()>=0 );
        lock.free();
        cellDescription.setExtrapolatedSolutionCompressed( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedSolutionCompressedIndex()).data() ) );

        const int numberOfEntries = getDataPerPatchBoundary();
        tearApart(numberOfEntries, cellDescription.getExtrapolatedSolutionIndex(), cellDescription.getExtrapolatedSolutionCompressedIndex(), compressionOfExtrapolatedSolution);

        #if defined(TrackGridStatistics)
        lock.lock();
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += CompressedDataHeap::getInstance().getData( cellDescription.getExtrapolatedSolutionCompressedIndex() ).size();
        lock.free();
        #endif

        #if !defined(ValidateCompressedVsUncompressedData)
        lock.lock();
        DataHeap::getInstance().deleteData( cellDescription.getExtrapolatedSolutionIndex(), true );
        cellDescription.setExtrapolatedSolutionIndex(-1);
        cellDescription.setExtrapolatedSolution(nullptr);
        lock.free();
        #endif
      }
      else {
        #if defined(TrackGridStatistics)
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        PipedUncompressedBytes += getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).size() * 8.0;
        PipedCompressedBytes   += getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex()).size() * 8.0;
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
    cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPerCell,dataPerCell,DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setSolutionIndex( DataHeap::getInstance().createData( dataPerCell,dataPerCell,DataHeap::Allocation::UseOnlyRecycledEntries) );
    cellDescription.setExtrapolatedSolutionIndex( DataHeap::getInstance().createData( dataPerBoundary,dataPerBoundary,DataHeap::Allocation::UseOnlyRecycledEntries) );
    lock.free();

    if (cellDescription.getPreviousSolutionIndex()==-1) { // allocate new array if recycling has failed
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setPreviousSolutionIndex( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setSolutionIndex( DataHeap::getInstance().createData( dataPerCell, dataPerCell ) );
      lock.free();
    }
    if (cellDescription.getExtrapolatedSolutionIndex()==-1) {
      ensureAllJobsHaveTerminated(JobType::SkeletonJob);
      ensureAllJobsHaveTerminated(JobType::EnclaveJob);
      lock.lock();
      cellDescription.setExtrapolatedSolutionIndex( DataHeap::getInstance().createData(dataPerBoundary, dataPerBoundary ) );
      lock.free();
    }
    cellDescription.setPreviousSolution    ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getPreviousSolutionIndex()).data() ) );
    cellDescription.setSolution            ( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getSolutionIndex()).data() ) );
    cellDescription.setExtrapolatedSolution( static_cast<void*>(CompressedDataHeap::getInstance().getData(cellDescription.getExtrapolatedSolutionIndex()).data() ) );
  }
  #else
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ));
  assertion( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ));
  #endif

  assertion1(
      CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ),
      cellDescription.getPreviousSolutionCompressed()
    );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ),
    cellDescription.getSolutionCompressed()
  );
  assertion1(
    CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionCompressedIndex() ),
    cellDescription.getExtrapolatedSolutionCompressed()
  );

  peano::datatraversal::TaskSet glueTasks(
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInPreviousSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionIndex() ), cellDescription.getPreviousSolutionIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getPreviousSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
        glueTogether(numberOfEntries, cellDescription.getPreviousSolutionIndex(), cellDescription.getPreviousSolutionCompressedIndex(), cellDescription.getBytesPerDoFInPreviousSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getPreviousSolutionCompressedIndex(), true );
        cellDescription.setPreviousSolutionCompressedIndex(-1);
        cellDescription.setPreviousSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getSolutionIndex() ), cellDescription.getSolutionIndex() );
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerPatch() + getGhostDataPerPatch();
        glueTogether(numberOfEntries, cellDescription.getSolutionIndex(), cellDescription.getSolutionCompressedIndex(), cellDescription.getBytesPerDoFInSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getSolutionCompressedIndex(), true );
        cellDescription.setSolutionCompressedIndex(-1);
        cellDescription.setSolutionCompressed(nullptr);
        lock.free();
      }
      return false;
    },
    [&]() -> bool {
      if (cellDescription.getBytesPerDoFInExtrapolatedSolution()<7) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ), cellDescription.getExtrapolatedSolutionIndex());
        assertion( CompressedDataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionCompressedIndex() ));
        const int numberOfEntries = getDataPerPatchBoundary();
        glueTogether(numberOfEntries, cellDescription.getExtrapolatedSolutionIndex(), cellDescription.getExtrapolatedSolutionCompressedIndex(), cellDescription.getBytesPerDoFInExtrapolatedSolution());
        tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
        CompressedDataHeap::getInstance().deleteData( cellDescription.getExtrapolatedSolutionCompressedIndex(), true );
        cellDescription.setExtrapolatedSolutionCompressedIndex(-1);
        cellDescription.setExtrapolatedSolutionCompressed(nullptr);
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
      peano::datatraversal::TaskSet spawned( new CompressionJob( *this, cellDescription, jobCounter ));
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

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionIndex()), cellDescription.toString() );

  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getPreviousSolutionAveragesIndex()), cellDescription.toString() );
  assertion1( DataHeap::getInstance().isValidIndex(cellDescription.getExtrapolatedSolutionAveragesIndex()), cellDescription.toString() );

  const int dataPerSubcell   = getNumberOfParameters()+getNumberOfVariables();
  const int subcellsPerPatch = (getDataPerPatch()+getGhostDataPerPatch())/ dataPerSubcell;
  const int subcellsPerFace  = (getDataPerPatchFace()) / dataPerSubcell;

  auto& solutionAverages             = getDataHeapEntries( cellDescription.getSolutionAveragesIndex() );
  auto& previousSolutionAverage      = getDataHeapEntries( cellDescription.getPreviousSolutionAveragesIndex() );
  auto& extrapolatedSolutionAverages = getDataHeapEntries( cellDescription.getExtrapolatedSolutionAveragesIndex() );

  // patch data
  kernels::idx2 idx_patchData    (subcellsPerPatch,dataPerSubcell);
  for (int i=0; i<subcellsPerPatch; i++) {
    for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
      solutionAverages[variableNumber]        += getDataHeapEntries(cellDescription.getSolutionIndex())        [idx_patchData(i,variableNumber)];
      previousSolutionAverage[variableNumber] += getDataHeapEntries(cellDescription.getPreviousSolutionIndex())[idx_patchData(i,variableNumber)];
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
            getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex())[idx_faceData(face,i,variableNumber)];
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
      getDataHeapEntries(cellDescription.getSolutionIndex())
          [idx_patchData(i,variableNumber)] += sign * getDataHeapEntries(cellDescription.getSolutionAveragesIndex())[variableNumber];

      getDataHeapEntries(cellDescription.getPreviousSolutionIndex())
          [idx_patchData(i,variableNumber)] += sign * getDataHeapEntries(cellDescription.getPreviousSolutionAveragesIndex())[variableNumber];
    }
  }

  // face data
  kernels::idx2 idx_faceDataAvg(DIMENSIONS_TIMES_TWO,dataPerSubcell);
  kernels::idx3 idx_faceData   (DIMENSIONS_TIMES_TWO,subcellsPerFace,dataPerSubcell);
  for (int face=0; face<2*DIMENSIONS; face++) {
    for (int i=0; i<subcellsPerFace; i++) {
      for (int variableNumber=0; variableNumber<dataPerSubcell; variableNumber++) {
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionIndex() ), cellDescription.toString() );
        assertion1( DataHeap::getInstance().isValidIndex( cellDescription.getExtrapolatedSolutionAveragesIndex() ), cellDescription.toString() );
        getDataHeapEntries(cellDescription.getExtrapolatedSolutionIndex())[idx_faceData(face,i,variableNumber)] +=
              sign * getDataHeapEntries(cellDescription.getExtrapolatedSolutionAveragesIndex())[idx_faceDataAvg(face,variableNumber)];
      }
    }
  }
}

exahype::solvers::FiniteVolumesSolver::CompressionJob::CompressionJob(
  const FiniteVolumesSolver& solver,
  CellDescription&           cellDescription,
  const bool                 isSkeletonJob)
  :
  tarch::multicore::jobs::Job(Solver::getTaskType(isSkeletonJob),0),
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


bool exahype::solvers::FiniteVolumesSolver::CompressionJob::run() {
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
  CellDescription&         cellDescription,
  const int                cellDescriptionsIndex,
  const bool               isSkeletonJob):
  tarch::multicore::jobs::Job(Solver::getTaskType(isSkeletonJob),0),
  _solver(solver),
  _cellDescription(cellDescription),
  _cellDescriptionsIndex(cellDescriptionsIndex),
  _isSkeletonJob(isSkeletonJob) {
  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    int& jobCounter = (_isSkeletonJob) ? NumberOfSkeletonJobs : NumberOfEnclaveJobs;
    jobCounter++;
  }
  lock.free();
}

bool exahype::solvers::FiniteVolumesSolver::FusedTimeStepJob::run() {
  _solver.updateBody(_cellDescription,_cellDescriptionsIndex,false,false,_isSkeletonJob,false/*uncompressBefore*/);
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
  tarch::multicore::jobs::Job(Solver::getTaskType(false),0),
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

bool exahype::solvers::FiniteVolumesSolver::AdjustSolutionDuringMeshRefinementJob::run() {
  _solver.adjustSolutionDuringMeshRefinementBody(_cellDescription,_isInitialMeshRefinement);

  tarch::multicore::Lock lock(exahype::BackgroundJobSemaphore);
  {
    NumberOfAMRBackgroundJobs--;
    assertion( NumberOfAMRBackgroundJobs>=0 );
  }
  lock.free();
  return false;
}
