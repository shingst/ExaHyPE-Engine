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
 **/

#include "exahype/Cell.h"
#include "exahype/State.h"

#include "tarch/la/ScalarOperations.h"

#include "peano/utils/Loop.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "kernels/KernelCalls.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"


tarch::logging::Log exahype::Cell::_log("exahype::Cell");

#ifdef Parallel
#include <deque>

int exahype::Cell::ReceivedMetadataHeapIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

std::deque<int> exahype::Cell::ReceivedDataHeapIndices;
#endif


exahype::Cell::Cell() : Base() {
  // We initialise cells which are not touched by the
  // createCell(...) events of Peano's spacetree traversal automaton
  // with default ("do-nothing") values.
  _cellData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value)
: Base(value) {
  _cellData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::PersistentCell& argument) : Base(argument) {
  // This constructor is used to create a cell from persistent data.
  // Do not use it. This would overwrite persistent data.
}

void exahype::Cell::validateThatAllNeighbourMergesHaveBeenPerformed(
    const int cellDescriptionsIndex,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  auto& ADERDGPatches = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex);
  auto& FVPatches     = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex);

  // ADER-DG
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    switch (solver->getType()) {
      case solvers::Solver::Type::ADERDG: {
        const int element = solvers::Solver::indexOfCellDescription(ADERDGPatches,solverNumber);
        if (element!=solvers::Solver::NotFound) {
          const auto& patch = ADERDGPatches[element];
          bool allNeighbourMergesHaveBeenPerformed =
              (patch.getType()!=exahype::solvers::ADERDGSolver::CellDescription::Type::Cell ||
               patch.getNeighbourMergePerformed().all());

          assertion1( allNeighbourMergesHaveBeenPerformed, patch.toString() );
          if ( !allNeighbourMergesHaveBeenPerformed ) {
            logError("validateThatAllNeighbourMergesHaveBeenPerformed(...)",
                     "Failed for cellDescriptionsIndex="<<cellDescriptionsIndex<<". "<<
                     "Not all neighbour merges have been performed for ADERDGCellDescription="<<
                     patch.toString());
            std::terminate();
          }
        }
      } break;
      case solvers::Solver::Type::LimitingADERDG: {
        const int solverElement = solvers::Solver::indexOfCellDescription(ADERDGPatches,solverNumber);
        auto* limitingADERDG = static_cast<solvers::LimitingADERDGSolver*>(solver);
        if (solverElement!=solvers::Solver::NotFound) {
          const auto& solverPatch = ADERDGPatches[solverElement];
          bool allNeighbourMergesHaveBeenPerformed = true;
          if (
              solverPatch.getLevel()==limitingADERDG->getMaximumAdaptiveMeshLevel() &&
              solverPatch.getRefinementStatus() >= limitingADERDG->getSolver()->getMinimumRefinementStatusForActiveFVPatch()
          ) {
            const int limiterElement = solvers::Solver::indexOfCellDescription(FVPatches,solverNumber);
            assertion1(limiterElement!=solvers::Solver::NotFound,solverPatch.toString());
            if ( limiterElement==solvers::Solver::NotFound ) {
              logError("validateThatAllNeighbourMergesHaveBeenPerformed(...)",
                       "No limiter patch allocated for solver patch=" << solverPatch.toString());
              std::terminate();
            }
            const auto& limiterPatch = FVPatches[limiterElement];
            allNeighbourMergesHaveBeenPerformed &= limiterPatch.getNeighbourMergePerformed().all();
          } else {
            allNeighbourMergesHaveBeenPerformed &=
              (solverPatch.getType()!=exahype::solvers::ADERDGSolver::CellDescription::Type::Cell ||
               solverPatch.getNeighbourMergePerformed().all());
          }

          assertion1( allNeighbourMergesHaveBeenPerformed, solverPatch.toString() );
          if ( !allNeighbourMergesHaveBeenPerformed ) {
            logError("validateThatAllNeighbourMergesHaveBeenPerformed(...)",
                     "Failed for cellDescriptionsIndex="<<cellDescriptionsIndex<<". "<<
                     "Not all neighbour merges have been performed for ADERDGCellDescription"<<
                     solverPatch.toString());
            std::terminate();
          }
        }
      } break;
      case solvers::Solver::Type::FiniteVolumes: {
        const int element = solvers::Solver::indexOfCellDescription(FVPatches,solverNumber);
        if (element!=solvers::Solver::NotFound) {
          const auto& patch = FVPatches[element];
          bool allNeighbourMergesHaveBeenPerformed = patch.getNeighbourMergePerformed().all();

          assertion1( allNeighbourMergesHaveBeenPerformed, patch.toString() );
          if ( !allNeighbourMergesHaveBeenPerformed ) {
            logError("validateThatAllNeighbourMergesHaveBeenPerformed(...)",
                     "Failed for cellDescriptionsIndex="<<cellDescriptionsIndex<<". "<<
                     "Not all neighbour merges have been performed for ADERDGCellDescription="<<
                     patch.toString());
            std::terminate();
          }
        }
      } break;
      default:
        break;
    }
  }
}

void exahype::Cell::resetNeighbourMergeFlags(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex));

  // ADER-DG
  for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
    for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
      p.setNeighbourMergePerformed(faceIndex,false);

      #ifdef Parallel
      int listingsOfRemoteRank =
          countListingsOfRemoteRankAtInsideFace(
              faceIndex,fineGridVertices,fineGridVerticesEnumerator);
      if (listingsOfRemoteRank==0) {
        listingsOfRemoteRank = TWO_POWER_D;
      }
      p.setFaceDataExchangeCounter(faceIndex,listingsOfRemoteRank);
      assertion(p.getFaceDataExchangeCounter(faceIndex)>0);
      #endif
    }
  }

  // Finite-Volumes (loop body can be copied from ADER-DG loop)
  for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex)) {
    for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
      p.setNeighbourMergePerformed(faceIndex,false);

      #ifdef Parallel
      int listingsOfRemoteRank =
          countListingsOfRemoteRankAtInsideFace(
              faceIndex,fineGridVertices,fineGridVerticesEnumerator);
      if (listingsOfRemoteRank==0) {
        listingsOfRemoteRank = TWO_POWER_D;
      }
      p.setFaceDataExchangeCounter(faceIndex,listingsOfRemoteRank);
      assertion(p.getFaceDataExchangeCounter(faceIndex)>0); // TODO Info can be used to determine who is at boundary from vertex view
      #endif
    }
  }
}

std::bitset<DIMENSIONS_TIMES_TWO> exahype::Cell::determineInsideAndOutsideFaces(
    const exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  std::bitset<DIMENSIONS_TIMES_TWO> isInside;
  isInside.reset();

  for (int direction=0; direction<DIMENSIONS; direction++) {
    for (int orientation=0; orientation<2; orientation++) {
      const int faceIndex = 2*direction+orientation;
      isInside[faceIndex]=false;

      // Works for a single solver
      dfor2(v) // Loop over vertices.
      if (v(direction) == orientation) {
        isInside[faceIndex] =
            isInside[faceIndex] ||
            verticesAroundCell[ verticesEnumerator(v) ].isInside();
      }
      enddforx // v
    }
  }

  return isInside;
}

tarch::la::Vector<DIMENSIONS,double> exahype::Cell::computeFaceBarycentre(
    const tarch::la::Vector<DIMENSIONS,double>& cellOffset,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const int direction, const int orientation) {
  tarch::la::Vector<DIMENSIONS,double> faceBarycentre;
  for (int i=0; i<DIMENSIONS; i++) {
    faceBarycentre[i] = cellOffset[i] + 0.5 * cellSize[i];
  }
  faceBarycentre[direction] = cellOffset[direction] + orientation * cellSize[direction];

  return faceBarycentre;
}

bool exahype::Cell::isAtRemoteBoundary(
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  bool result = false;
  #ifdef Parallel
  tarch::la::Vector<DIMENSIONS,int> center(1);
  dfor2(v) // Loop over vertices.
    if (verticesAroundCell[ verticesEnumerator(v) ].isAdjacentToRemoteRank()) {
      dfor2(a) // Loop over adjacent ranks. Does also include own rank.
        result |= tarch::la::countEqualEntries(v+a,center)==DIMENSIONS-1 && // offset in one direction from center=>face neighbour
                  verticesAroundCell[ verticesEnumerator(v) ].isInside() && // exclude boundary and outside vertices
                  verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]!=
                      tarch::parallel::Node::getInstance().getRank();
      enddforx //a
    }
  enddforx // v
  #endif
  return result;
}

void exahype::Cell::setupMetaData() {
  assertion1(!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),toString());

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
    const int cellDescriptionIndex = exahype::solvers::ADERDGSolver::Heap::getInstance().createData(0, 0);
    assertion2(!exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(cellDescriptionIndex),cellDescriptionIndex,toString());
    exahype::solvers::FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionIndex,0,0);
    _cellData.setCellDescriptionsIndex(cellDescriptionIndex);
  lock.free();
}

void exahype::Cell::shutdownMetaDataAndResetCellDescriptionsIndex() {
  shutdownMetaData();
  _cellData.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

void exahype::Cell::shutdownMetaData() const {
  assertion1(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),toString());

  // use lock inside
  exahype::solvers::ADERDGSolver::eraseCellDescriptions(_cellData.getCellDescriptionsIndex());
  exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(_cellData.getCellDescriptionsIndex());

  tarch::multicore::Lock lock(exahype::HeapSemaphore);
    exahype::solvers::ADERDGSolver::Heap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());
    exahype::solvers::FiniteVolumesSolver::Heap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());
  lock.free();
}

bool exahype::Cell::isEmpty() const {
  if (_cellData.getCellDescriptionsIndex()!=multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    assertion1( exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
                _cellData.getCellDescriptionsIndex());
    assertion1( exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
                _cellData.getCellDescriptionsIndex());
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  return
      exahype::solvers::ADERDGSolver::Heap::getInstance().getData(_cellData.getCellDescriptionsIndex()).empty() &&
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(_cellData.getCellDescriptionsIndex()).empty();
}

bool exahype::Cell::isInitialised() const {
  if ( _cellData.getCellDescriptionsIndex() >= 0 ) {
    assertion1( exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
                _cellData.getCellDescriptionsIndex());
    assertion1( exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
                _cellData.getCellDescriptionsIndex());
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  return _cellData.getCellDescriptionsIndex() >= 0;
}

int exahype::Cell::getCellDescriptionsIndex() const {
  return _cellData.getCellDescriptionsIndex();
}

void exahype::Cell::setCellDescriptionsIndex(int cellDescriptionsIndex) {
  _cellData.setCellDescriptionsIndex(cellDescriptionsIndex);
}

void exahype::Cell::addNewCellDescription(
    const int solverNumber,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
    const exahype::records::FiniteVolumesCellDescription::RefinementEvent refinementEvent,
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
    const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  exahype::solvers::FiniteVolumesSolver::addNewCellDescription(
      _cellData.getCellDescriptionsIndex(),solverNumber,
      cellType,refinementEvent,
      level,parentIndex,cellSize,cellOffset);
}


void exahype::Cell::addNewCellDescription(
    const int                                     solverNumber,
    const exahype::records::ADERDGCellDescription::Type cellType,
    const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
    const int                                     level,
    const int                                     parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
    const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  exahype::solvers::ADERDGSolver::addNewCellDescription(
      _cellData.getCellDescriptionsIndex(),solverNumber,
      cellType,refinementEvent,
      level,parentIndex,cellSize,cellOffset);
}

int exahype::Cell::getNumberOfADERDGCellDescriptions() const {
  return exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


int exahype::Cell::getNumberOfFiniteVolumeCellDescriptions() const {
  return exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


#ifdef Parallel
int exahype::Cell::countListingsOfRemoteRankAtInsideFace(
    const int faceIndex,
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  int result = 0;

  const int f = faceIndex % 2;   // "0" indicates a left face, "1" indicates a right face.
  const int d = (faceIndex-f)/2; // The normal direction: 0: x, 1: y, 1: z.

  tarch::la::Vector<DIMENSIONS,int> pos(1); // This is now the center, i.e., (1,1,...,1).
  pos(d) = 2*f;                             // This is a shift from the center by one unit in direction d.

  int faceNeighbourRank = -1; // This variable is introduced to make sure that the adjacent remote rank is unique.
  // TODO(Dominic): Uniqueness is probably guaranteed by the SFC based DD.
  dfor2(v) // Loop over vertices.
  dfor2(a) // Loop over adjacent ranks. Does also include own rank.
  if (tarch::la::equals(v+a,pos)                             &&
      verticesAroundCell[ verticesEnumerator(v) ].isInside() &&
      verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]!=
          tarch::parallel::Node::getInstance().getRank()) {
    // Increment
    if (faceNeighbourRank==-1) {
      faceNeighbourRank = verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar];
    }
    if (verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]==faceNeighbourRank) {
      result++;
    }
  }
  enddforx // a
  enddforx // v

  // result corresponds either to no connection, edge connection, or whole face connection.
  // If the bounding box scaling is turned off, there might further be only a single inside point possible
  // in 3d. result is then 1.
#ifdef Asserts
  std::stringstream message;
  message << std::endl;
  bool foundHangingNode = false;
  dfor2(v)
    foundHangingNode |= verticesAroundCell[ verticesEnumerator(v) ].isHangingNode();
    message << "v="<<v.toString() <<
    ", hanging node=" << verticesAroundCell[ verticesEnumerator(v) ].isHangingNode() <<
    ", location=" <<
    exahype::Vertex::Records::toString(verticesAroundCell[ verticesEnumerator(v) ].getRecords().getInsideOutsideDomain()) <<
    ": adjacentRanks="<<verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks().toString() << std::endl;
  enddforx
#endif
  assertion5(foundHangingNode||result==0||result==TWO_POWER_D_DIVIDED_BY_TWO/4||result==TWO_POWER_D_DIVIDED_BY_TWO/2||result==TWO_POWER_D_DIVIDED_BY_TWO,
             result,pos.toString(),faceIndex,tarch::parallel::Node::getInstance().getRank(),message.str());
  return result;
}

bool exahype::Cell::hasToCommunicate( const int level ) const {
  return
      level >= exahype::solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
}

// MASTER->WORKER

void exahype::Cell::broadcastGlobalDataToWorker(
    const int                                   worker,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    solver->sendDataToWorker(worker,cellCentre,level);
  }

  for (auto& plotter : exahype::plotters::RegisteredPlotters) {
    plotter->sendDataToWorker(
        worker,cellCentre,level);
  }
}

void exahype::Cell::mergeWithGlobalDataFromMaster(
    const int                                   master,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    solver->mergeWithMasterData(master,cellCentre,level);
  }

  for (auto& plotter : exahype::plotters::RegisteredPlotters) {
    plotter->mergeWithMasterData(master,cellCentre,level);
  }
}


// WORKER->MASTER

void exahype::Cell::reduceGlobalDataToMaster(
    const int                                   master,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    solver->sendDataToMaster(master,cellCentre,level);
  }
}

void exahype::Cell::mergeWithGlobalDataFromWorker(
    const int                                   worker,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level) {
  for (auto& solver : exahype::solvers::RegisteredSolvers) {
    solver->mergeWithWorkerData(worker,cellCentre,level);
  }
}

#endif
