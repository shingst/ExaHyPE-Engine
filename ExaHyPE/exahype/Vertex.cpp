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
 
#include "exahype/Vertex.h"

#include "tarch/la/Scalar.h"

#include "tarch/multicore/Loop.h"
#include "tarch/multicore/dForRange.h"

#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/State.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

tarch::logging::Log exahype::Vertex::_log( "exahype::Vertex");

exahype::Vertex::Vertex() : Base() {
  _vertexData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance()
          .createVertexLinkMapForNewVertex() );
}

exahype::Vertex::Vertex(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
  // do nothing
}

exahype::Vertex::Vertex(const Base::PersistentVertex& argument)
    : Base(argument) {
  // do nothing
}

bool exahype::Vertex::equalUpToRelativeTolerance(
    const tarch::la::Vector<DIMENSIONS,double>& lhs,
    const tarch::la::Vector<DIMENSIONS,double>& rhs) {
  const double tolerance =
      tarch::la::NUMERICAL_ZERO_DIFFERENCE *
      std::max(
          1.0, std::max( tarch::la::maxAbs(lhs), tarch::la::maxAbs(rhs) )
  );
  return tarch::la::equals( lhs, rhs, tolerance );
}

tarch::la::Vector<TWO_POWER_D, int>
exahype::Vertex::getCellDescriptionsIndex() const {
  return _vertexData.getCellDescriptionsIndex();
}

tarch::la::Vector<DIMENSIONS,double> exahype::Vertex::computeFaceBarycentre(
    const tarch::la::Vector<DIMENSIONS,double>& x,
    const tarch::la::Vector<DIMENSIONS,double>& h,
    const int                                   normalDirection,
    const tarch::la::Vector<DIMENSIONS,int>&    cellPosition) {
  tarch::la::Vector<DIMENSIONS,double> barycentre;
  for (int d=0; d<DIMENSIONS; d++) {
    barycentre(d) = x(d) + 0.5*h(d)*(2*cellPosition(d)-1);
  }
  barycentre(normalDirection) = x(normalDirection);
  return barycentre;
}

exahype::solvers::Solver::RefinementControl exahype::Vertex::evaluateRefinementCriterion(
    const tarch::la::Vector<DIMENSIONS, double>& vertexOffset,
    const tarch::la::Vector<DIMENSIONS, double>& level,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const bool checkThoroughly) const {
  bool canErase   = true;
  bool mustRefine = false;
  dfor2(pos)
    tarch::la::Vector<DIMENSIONS,double> cellOffset(vertexOffset);
    for (int d = 0; d < DIMENSIONS; ++d) {
      cellOffset[d] += ( pos[d]-1 ) * cellSize[d];
    }

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      const int cellDescriptionsIndex = _vertexData.getCellDescriptionsIndex(posScalar);
      exahype::solvers::Solver::RefinementControl control =
          solver->eraseOrRefineAdjacentVertices(cellDescriptionsIndex,solverNumber,cellOffset,cellSize,checkThoroughly);
      canErase   &= (control==exahype::solvers::Solver::RefinementControl::Erase);
      mustRefine |= (control==exahype::solvers::Solver::RefinementControl::Refine);
    }
  enddforx

  if (mustRefine) {
    return exahype::solvers::Solver::RefinementControl::Refine;
  } else if (canErase) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  } else {
    return exahype::solvers::Solver::RefinementControl::Keep;
  }
}

void exahype::Vertex::mergeOnlyNeighboursMetadataLoopBody(
    const int pos1Scalar,
    const int pos2Scalar,
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const exahype::State::AlgorithmSection& section,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const bool                                   checkThoroughly) {
  tarch::la::Vector<DIMENSIONS,int> pos1 = delineariseIndex2(pos1Scalar);
  tarch::la::Vector<DIMENSIONS,int> pos2 = delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  bool validIndex1 = cellDescriptionsIndex1 >= 0;
  bool validIndex2 = cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1));
  assertion(cellDescriptionsIndex2 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2));

  if ( validIndex1 && validIndex2 ) {
    auto& ADERDGPatches1 = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1);
    auto& FVPatches1     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1);

    auto& ADERDGPatches2 = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2);
    auto& FVPatches2     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2);

    bool mergeNeighbours = ( !ADERDGPatches1.empty() && !ADERDGPatches2.empty() ) ||
                           ( !FVPatches1.empty() && !FVPatches2.empty() );
    if ( mergeNeighbours ) {
      for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        if (solver->isMergingMetadata(section)) {
          switch ( solver->getType() ) {
            case solvers::Solver::Type::ADERDG:
              static_cast<solvers::ADERDGSolver*>(solver)->
                mergeNeighboursMetadata(ADERDGPatches1,ADERDGPatches2,solverNumber,pos1,pos2,x,h,checkThoroughly);
              break;
            case solvers::Solver::Type::LimitingADERDG:
              static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
              mergeNeighboursMetadata(ADERDGPatches1,ADERDGPatches2,solverNumber,pos1,pos2,x,h,checkThoroughly);
              break;
            case solvers::Solver::Type::FiniteVolumes:
              // do nothing
              break;
            default:
              assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              logError("mergeWithBoundaryDataIfNotDoneYet(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              std::abort();
              break;
          }
        }
      }
    }
  }
}

void exahype::Vertex::mergeOnlyNeighboursMetadata(
    const exahype::State::AlgorithmSection& section,
    const tarch::la::Vector<DIMENSIONS,     double>& x,
    const tarch::la::Vector<DIMENSIONS,     double>& h,
    const bool                              checkThoroughly) const {
  assertion(!isHangingNode());
  assertion(isInside() || isBoundary());

  #if DIMENSIONS==2
  mergeOnlyNeighboursMetadataLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),section,x,h,checkThoroughly);
  #elif DIMENSIONS==3
  mergeOnlyNeighboursMetadataLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(0,4,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(4),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(1,5,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(5),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(2,6,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(6),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(3,7,_vertexData.getCellDescriptionsIndex(3),_vertexData.getCellDescriptionsIndex(7),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(4,5,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(5),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(4,6,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(6),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(5,7,_vertexData.getCellDescriptionsIndex(5),_vertexData.getCellDescriptionsIndex(7),section,x,h,checkThoroughly);
  mergeOnlyNeighboursMetadataLoopBody(6,7,_vertexData.getCellDescriptionsIndex(6),_vertexData.getCellDescriptionsIndex(7),section,x,h,checkThoroughly);
  #endif
}

void exahype::Vertex::validateNeighbourhood(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2) {
  solvers::Solver::InterfaceInfo face(pos1,pos2);

  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    switch (solver->getType()) {
    case exahype::solvers::Solver::Type::LimitingADERDG:
    case exahype::solvers::Solver::Type::ADERDG: {
      // Cell 1
      const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
      if (element1!=exahype::solvers::Solver::NotFound) {
        auto& p1 = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex1,element1);
        if (
            (p1.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell ||
            p1.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Ancestor)
            &&
            cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
        ) {
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<face._faceIndex1<<" next to empty cell: cell="<<p1.toString());
          std::terminate();
        }
      }
      // Cell 2
      const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
      if (element2!=exahype::solvers::Solver::NotFound) {
        auto& p2 = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex2,element2);
        if (
            (p2.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell ||
            p2.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Ancestor)
            &&
            cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
        ) {
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex2<<" is at face="<<face._faceIndex2<<" next to empty cell: cell="<<p2.toString());
          std::terminate();
        }
      }
    } break;
    case exahype::solvers::Solver::Type::FiniteVolumes: {
      // Cell 1
      const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
      if (element1!=exahype::solvers::Solver::NotFound) {
        auto& p1 = exahype::solvers::FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex1,element1);
        if (
            p1.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Cell
            &&
            cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
        ) {
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<face._faceIndex1<<" next to empty cell: cell="<<p1.toString());
          std::terminate();
        }
      }
      // Cell 2
      const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
      if (element2!=exahype::solvers::Solver::NotFound) {
        auto& p2 = exahype::solvers::FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex2,element2);
        if (
            p2.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Cell
            &&
            cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
        ) {
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex2<<" is at face="<<face._faceIndex2<<" next to empty cell: cell="<<p2.toString());
          std::terminate();
        }
      }
    } break;
    }
  }
}

void exahype::Vertex::mergeWithBoundaryData(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  tarch::la::Vector<DIMENSIONS,int> posCell     = pos1;
  tarch::la::Vector<DIMENSIONS,int> posBoundary = pos2;
  int cellDescriptionsIndex                     = cellDescriptionsIndex1;
  if ( cellDescriptionsIndex2 >= 0 ) {
    posCell     = pos2;
    posBoundary = pos1;
    cellDescriptionsIndex = cellDescriptionsIndex2;
  }

  auto& ADERDGPatches = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex);
  auto& FVPatches     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex);
  bool mergeWithBoundaryData = !ADERDGPatches.empty() || !FVPatches.empty();

  if (mergeWithBoundaryData) {
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      switch ( solver->getType() ) {
        case solvers::Solver::Type::ADERDG:
          static_cast<solvers::ADERDGSolver*>(solver)->
            mergeWithBoundaryData(ADERDGPatches,solverNumber,posCell,posBoundary);
          break;
        case solvers::Solver::Type::LimitingADERDG:
          static_cast<solvers::LimitingADERDGSolver*>(solver)->
            mergeWithBoundaryData(ADERDGPatches,FVPatches,solverNumber,posCell,posBoundary,false);
          break;
        case solvers::Solver::Type::FiniteVolumes:
          static_cast<solvers::FiniteVolumesSolver*>(solver)->
            mergeWithBoundaryData(FVPatches,solverNumber,posCell,posBoundary);
          break;
        default:
          assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          logError("mergeWithBoundaryDataIfNotDoneYet(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          std::abort();
          break;
      }
    }
  }
}

void exahype::Vertex::mergeNeighboursDataAndMetadata(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  auto& ADERDGPatches1 = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1);
  auto& FVPatches1     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1);

  auto& ADERDGPatches2 = solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2);
  auto& FVPatches2     = solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2);

  bool mergeNeighbours = ( !ADERDGPatches1.empty() && !ADERDGPatches2.empty() ) ||
                         ( !FVPatches1.empty() && !FVPatches2.empty() );

  if ( mergeNeighbours ) {
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      switch ( solver->getType() ) {
        case solvers::Solver::Type::ADERDG:
          static_cast<solvers::ADERDGSolver*>(solver)->
            mergeNeighboursData(ADERDGPatches1,ADERDGPatches2,solverNumber,pos1,pos2);
          static_cast<solvers::ADERDGSolver*>(solver)->
            mergeNeighboursMetadata(ADERDGPatches1,ADERDGPatches2,solverNumber,pos1,pos2,x,h,false);
          break;
        case solvers::Solver::Type::LimitingADERDG:
          static_cast<solvers::LimitingADERDGSolver*>(solver)->
            mergeNeighboursData(ADERDGPatches1,ADERDGPatches2,FVPatches1,FVPatches2,solverNumber,pos1,pos2,false);
          static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
              mergeNeighboursMetadata(ADERDGPatches1,ADERDGPatches2,solverNumber,pos1,pos2,x,h,false);
          break;
        case solvers::Solver::Type::FiniteVolumes:
          static_cast<solvers::FiniteVolumesSolver*>(solver)->
            mergeNeighboursData(FVPatches1,FVPatches2,solverNumber,pos1,pos2);
          break;
        default:
          assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          logError("mergeWithBoundaryDataIfNotDoneYet(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
          std::abort();
          break;
      }
    }
  }
}

tarch::la::Vector<DIMENSIONS,int> exahype::Vertex::delineariseIndex2(int index) {
  #if DIMENSIONS==2
  switch( index ) {
    case 0:
      return tarch::la::Vector<DIMENSIONS,int>(0,0);
    case 1:
      return tarch::la::Vector<DIMENSIONS,int>(1,0);
    case 2:
      return tarch::la::Vector<DIMENSIONS,int>(0,1);
    case 3:
      return tarch::la::Vector<DIMENSIONS,int>(1,1);
    default:
      logError("delineariseIndex2(index)","index must be in range [0,3]!");
      return tarch::la::Vector<DIMENSIONS,int>(-1,-1);
  }
  #elif DIMENSIONS==3
  switch( index ) {
    case 0:
      return tarch::la::Vector<DIMENSIONS,int>(0,0,0);
    case 1:
      return tarch::la::Vector<DIMENSIONS,int>(1,0,0);
    case 2:
      return tarch::la::Vector<DIMENSIONS,int>(0,1,0);
    case 3:
      return tarch::la::Vector<DIMENSIONS,int>(1,1,1);
    case 4:
      return tarch::la::Vector<DIMENSIONS,int>(0,0,1);
    case 5:
      return tarch::la::Vector<DIMENSIONS,int>(1,0,1);
    case 6:
      return tarch::la::Vector<DIMENSIONS,int>(0,1,1);
    case 7:
      return tarch::la::Vector<DIMENSIONS,int>(1,1,1);
    default:
      logError("delineariseIndex2(index)","index must be in range [0,7]!");
      return tarch::la::Vector<DIMENSIONS,int>(-1,-1,-1);
  }
  #endif
}

void exahype::Vertex::mergeNeighboursLoopBody(
    const int pos1Scalar,
    const int pos2Scalar,
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS, double> x,
    const tarch::la::Vector<DIMENSIONS, double> h) {
  #if defined(Asserts) || defined (ValidateNeighbourHoodDuringNeighbourMerge)
  constexpr bool validate = true;
  #else
  constexpr bool validate = false;
  #endif
  assertion2(cellDescriptionsIndex1 < 1 || cellDescriptionsIndex1 != cellDescriptionsIndex2,cellDescriptionsIndex1,cellDescriptionsIndex2);

  tarch::la::Vector<DIMENSIONS,int> pos1 = delineariseIndex2(pos1Scalar);
  tarch::la::Vector<DIMENSIONS,int> pos2 = delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  bool validIndex1 = cellDescriptionsIndex1 >= 0;
  bool validIndex2 = cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || solvers::ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex1));
  assertion(cellDescriptionsIndex2 < 0 || solvers::ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex2));

  if ( validIndex1 && validIndex2 ) {
    mergeNeighboursDataAndMetadata(cellDescriptionsIndex1,cellDescriptionsIndex2,pos1,pos2,x,h);
  } else if (
      ((validIndex1 && !validIndex2 &&
          cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex)
          ||
          (!validIndex1 && validIndex2 &&
              cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex))
  ) {
    mergeWithBoundaryData(cellDescriptionsIndex1,cellDescriptionsIndex2,pos1,pos2,x,h);
  } else if  (
      validate
      &&
      validIndex1 != validIndex2
      &&
      cellDescriptionsIndex1!=multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex&&
      cellDescriptionsIndex2!=multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex
  ) {
    validateNeighbourhood(cellDescriptionsIndex1,cellDescriptionsIndex2,pos1,pos2);
  }
}

void exahype::Vertex::mergeNeighbours(
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  if ( tarch::la::allSmallerEquals(h,exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers()) ) {
    #if DIMENSIONS==2 and defined(SharedMemoryParallelisation) // TODO(Dominic): Comment back in if it works
    #if DIMENSIONS==2
    peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> bool {
      mergeNeighboursLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),x,h);
      return false;
    },
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    true);
    #elif DIMENSIONS==3
    peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> bool {
      mergeNeighboursLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(0,4,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(4),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(1,5,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(5),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(2,6,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(6),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(3,7,_vertexData.getCellDescriptionsIndex(3),_vertexData.getCellDescriptionsIndex(7),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(4,5,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(5),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(4,6,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(6),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(5,7,_vertexData.getCellDescriptionsIndex(5),_vertexData.getCellDescriptionsIndex(7),x,h);
      return false;
    },
    [&]() -> bool {
      mergeNeighboursLoopBody(6,7,_vertexData.getCellDescriptionsIndex(6),_vertexData.getCellDescriptionsIndex(7),x,h);
      return false;
    },
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
    true);
    #endif
    #else
    #if DIMENSIONS==2
    mergeNeighboursLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),x,h);
    mergeNeighboursLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),x,h);
    mergeNeighboursLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),x,h);
    mergeNeighboursLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),x,h);
    #elif DIMENSIONS==3
    mergeNeighboursLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),x,h);
    mergeNeighboursLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),x,h);
    mergeNeighboursLoopBody(0,4,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(4),x,h);
    mergeNeighboursLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),x,h);
    mergeNeighboursLoopBody(1,5,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(5),x,h);
    mergeNeighboursLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),x,h);
    mergeNeighboursLoopBody(2,6,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(6),x,h);
    mergeNeighboursLoopBody(3,7,_vertexData.getCellDescriptionsIndex(3),_vertexData.getCellDescriptionsIndex(7),x,h);
    mergeNeighboursLoopBody(4,5,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(5),x,h);
    mergeNeighboursLoopBody(4,6,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(6),x,h);
    mergeNeighboursLoopBody(5,7,_vertexData.getCellDescriptionsIndex(5),_vertexData.getCellDescriptionsIndex(7),x,h);
    mergeNeighboursLoopBody(6,7,_vertexData.getCellDescriptionsIndex(6),_vertexData.getCellDescriptionsIndex(7),x,h);
    #endif
    #endif
  }
}

// PARALLEL

#if Parallel
bool exahype::Vertex::hasToCommunicate( const int level) const {
  return
     isInside() &&
     level >= exahype::solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
}

bool exahype::Vertex::hasToSendMetadata(
  const int toRank,
  const tarch::la::Vector<DIMENSIONS,int>& src,
  const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return adjacentRanks[destScalar]   == toRank
         &&
         adjacentRanks[destScalar]   != tarch::parallel::Node::getGlobalMasterRank() &&
         adjacentRanks[srcScalar]    != tarch::parallel::Node::getGlobalMasterRank()
         &&
         (   // Send also when a fork/join was triggered for the current rank
             adjacentRanks[srcScalar]   == tarch::parallel::Node::getInstance().getRank() ||
             State::isForkTriggeredForRank(adjacentRanks[srcScalar]) ||
             State::isJoinTriggeredForRank(adjacentRanks[srcScalar])
         )
         &&
         tarch::la::countEqualEntries(dest, src) == (DIMENSIONS-1);
}

bool exahype::Vertex::hasToSendMetadataToNeighbour(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int srcCellDescriptionIndex = getCellDescriptionsIndex()[srcScalar];

  bool sendMetadataToNeighbour =
      exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex)
      &&
      (!exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionIndex).empty() ||
      !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionIndex).empty());

  if (sendMetadataToNeighbour) {
    solvers::Solver::BoundaryFaceInfo face(src,dest);

    if (!exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionIndex).empty()) {
      auto& patch = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionIndex)[0];
      tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
          exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),face._direction,face._orientation);
      tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
          exahype::Vertex::computeFaceBarycentre(x,h,face._direction,src);
      return equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
    }
    else if (!exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionIndex).empty()) {
      auto& patch = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionIndex)[0];
      tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch = // copy & paste from here
          exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),face._direction,face._orientation);
      tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
          exahype::Vertex::computeFaceBarycentre(x,h,face._direction,src);
      return equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
    } else {
      return false;
    }
  } {
    return false;
  }
}

void exahype::Vertex::sendOnlyMetadataToNeighbour(
    const int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    int level) const {
  if ( hasToCommunicate(level) ) {
    tarch::la::Vector<TWO_POWER_D, int> adjacentADERDGCellDescriptionsIndices = getCellDescriptionsIndex();
    dfor2(dest)
      dfor2(src)
        if ( hasToSendMetadata(toRank,src,dest) ) {
          const int srcCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(srcScalar);
          if ( hasToSendMetadataToNeighbour(src,dest,x,h) ) {
            exahype::sendNeighbourCommunicationMetadata(
                toRank,srcCellDescriptionIndex,src,dest,x,level);
          } else {
            exahype::sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(
                toRank,x,level);
          }
        }
      enddforx
    enddforx
  }
}



bool exahype::Vertex::hasToReceiveMetadata(
    const int fromRank,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const tarch::la::Vector<TWO_POWER_D,int> adjacentRanks = getAdjacentRanks();

  return
      adjacentRanks[srcScalar]    == fromRank
      &&
      adjacentRanks[srcScalar]    != tarch::parallel::Node::getGlobalMasterRank() &&
      adjacentRanks[destScalar]   != tarch::parallel::Node::getGlobalMasterRank()
      &&
      ( // Receive also when a fork/join is performed for the neighbour rank, i.e. such
        // an event was triggered in the iteration before
        adjacentRanks[destScalar] == tarch::parallel::Node::getInstance().getRank() ||
        State::isForkingRank(adjacentRanks[destScalar]) ||
        State::isJoiningRank(adjacentRanks[destScalar])
      )
      &&
      tarch::la::countEqualEntries(dest, src) == (DIMENSIONS-1);
}

bool exahype::Vertex::hasToMergeWithNeighbourMetadata(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const tarch::la::Vector<DIMENSIONS, double>& x) const {
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const int destCellDescriptionsIndex = _vertexData.getCellDescriptionsIndex(destScalar);

  // TODO(Dominic): Also have a thorough and simple version of this method

  bool mergeWithNeighbourMetadata =
      exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex) &&
      (!exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty() ||
      !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty());

  if ( mergeWithNeighbourMetadata ) {
    solvers::Solver::BoundaryFaceInfo face(dest,src);

    // ADER-DG
    bool mergeWithNeighbourData = true;
    for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) { // lookup is slow
      mergeWithNeighbourData &= p.getNeighbourMergePerformed(face._faceIndex)==false;
      p.setNeighbourMergePerformed(face._faceIndex,true);
    }

    // FV
    for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) { // lookup is slow
      mergeWithNeighbourData &= p.getNeighbourMergePerformed(face._faceIndex)==false;
      p.setNeighbourMergePerformed(face._faceIndex,true);
    }
    return mergeWithNeighbourData;
  } {
    return false;
  }
}

void exahype::Vertex::mergeOnlyWithNeighbourMetadata(
    const int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int level,
    const exahype::State::AlgorithmSection& section) const {
  if ( hasToCommunicate(level) ) {
    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;
        int destScalar = TWO_POWER_D - myDestScalar - 1;

        if ( hasToReceiveMetadata(fromRank,src,dest) ) {
          logDebug("mergeOnlyWithNeighbourMetadata(...)","[pre] rec. from rank "<<fromRank<<", x:"<<
                   x.toString() << ", level=" <<level << ", adjacentRanks: "
                   << getAdjacentRanks());

          exahype::MetadataHeap::HeapEntries receivedMetadata;
          receivedMetadata.clear();
          exahype::receiveNeighbourCommunicationMetadata(receivedMetadata,fromRank, x, level);
          assertionEquals(receivedMetadata.size(),
              exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

          if ( hasToMergeWithNeighbourMetadata(src,dest,x) ) {
            for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
              auto* solver = solvers::RegisteredSolvers[solverNumber];
              if (solver->isMergingMetadata(section)) {
                const int offset  = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;
                const int element = solver->tryGetElement(
                    getCellDescriptionsIndex()[destScalar],solverNumber);
                if (element!=exahype::solvers::Solver::NotFound) {
                  if (receivedMetadata[offset]!=exahype::InvalidMetadataEntry) {
                    MetadataHeap::HeapEntries metadataPortion(
                        receivedMetadata.begin()+offset,
                        receivedMetadata.begin()+offset+exahype::NeighbourCommunicationMetadataPerSolver);

                    solver->mergeWithNeighbourMetadata(
                        metadataPortion,
                        src, dest,
                        getCellDescriptionsIndex()[destScalar],element);
                  }
                }
              }
            }
          }
        }
      enddforx
    enddforx
  }
}

void exahype::Vertex::dropNeighbourMetadata(
    const int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int level) const {
  if ( hasToCommunicate(level) ) {
    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;
        //int destScalar = TWO_POWER_D - myDestScalar - 1;

        if (hasToReceiveMetadata(fromRank,src,dest)) {
          logDebug("dropNeighbourMetadata(...)","[pre] rec. from rank "<<fromRank<<", x:"<<
                   x.toString() << ", level=" <<level << ", adjacentRanks: "
                   << getAdjacentRanks());

          MetadataHeap::getInstance().receiveData(
              fromRank, x, level,
              peano::heap::MessageType::NeighbourCommunication);
        }
      enddforx
    enddforx
  }
}


bool exahype::Vertex::hasToSendDataToNeighbour(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int srcScalar  = peano::utils::dLinearisedWithoutLookup(src,2);
  const int srcCellDescriptionsIndex =
      _vertexData.getCellDescriptionsIndex(srcScalar);

  if ( srcCellDescriptionsIndex < 0 ) {
    return false; // !!! Make sure to consider all solver types here
  } else {
    assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex) );
    solvers::Solver::BoundaryFaceInfo face(src,dest);

    bool hasToSendDataToNeighbour = true;
    // ADER-DG
    for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionsIndex)) {
      // decrement counter beforehand
      int newCounterValue = p.getFaceDataExchangeCounter(face._faceIndex)-1;
      assertion2(newCounterValue>=0,newCounterValue,p.toString());
      assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
      p.setFaceDataExchangeCounter(face._faceIndex,newCounterValue);

      hasToSendDataToNeighbour &= p.getFaceDataExchangeCounter(face._faceIndex)==0;
    }

    // FV (copied from above)
    for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionsIndex)) {
      // decrement counter beforehand
      int newCounterValue = p.getFaceDataExchangeCounter(face._faceIndex)-1;
      assertion2(newCounterValue>=0,newCounterValue,p.toString());
      assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
      p.setFaceDataExchangeCounter(face._faceIndex,newCounterValue);

      hasToSendDataToNeighbour &= p.getFaceDataExchangeCounter(face._faceIndex)==0;
    }

    return hasToSendDataToNeighbour;
  }
}

bool exahype::Vertex::hasToMergeWithNeighbourData(
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest) const {
  const int destScalar  = peano::utils::dLinearisedWithoutLookup(dest,2);
  const int destCellDescriptionsIndex =
      _vertexData.getCellDescriptionsIndex(destScalar);

  // TODO(Dominic): Not as robust as the metadata variant

  if (
      !exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex)
      ||
      (exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty() &&
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty())
  ) {
    return false; // !!! Make sure to consider all solver types here
  } else {
    solvers::Solver::BoundaryFaceInfo face(dest,src);
  
    // ADER-DG
    bool mergeWithNeighbourData = true;
    for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) { // lookup is slow
      mergeWithNeighbourData &= p.getFaceDataExchangeCounter(face._faceIndex)==0;
  
      // now we can reset the flags
      if ( p.getFaceDataExchangeCounter(face._faceIndex)==0 ) {
        assertion(p.getNeighbourMergePerformed(face._faceIndex)==false);
        p.setFaceDataExchangeCounter(face._faceIndex,TWO_POWER_D); // TODO maybe do not do that here but in the cell? Can be used to determine which cell belongs to skeleton
        p.setNeighbourMergePerformed(face._faceIndex,true);
      }
    }
  
    // FV
    for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) { // lookup is slow
      mergeWithNeighbourData &= p.getFaceDataExchangeCounter(face._faceIndex)==0;
  
      // now we can reset the flags
      if ( p.getFaceDataExchangeCounter(face._faceIndex)==0 ) {
        assertion(p.getNeighbourMergePerformed(face._faceIndex)==false);
        p.setFaceDataExchangeCounter(face._faceIndex,TWO_POWER_D);
        p.setNeighbourMergePerformed(face._faceIndex,true);
      }
    }
    return mergeWithNeighbourData;
  }
}

void exahype::Vertex::sendEmptySolverDataToNeighbour(
    const int                                     toRank,
    const bool                                    sendMetadata,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  for ( auto* solver : exahype::solvers::RegisteredSolvers ) {
    solver->sendEmptyDataToNeighbour(toRank,x,level);
  }

  if ( sendMetadata ) {
    exahype::sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(
        toRank,x,level);
  }
}

void exahype::Vertex::sendSolverDataToNeighbour(
    const int                                    toRank,
    const bool                                   sendMetadata,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));
  assertion(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));

  for ( unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber ) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    const int element = solver->tryGetElement(srcCellDescriptionIndex,solverNumber);
    if ( element!=exahype::solvers::Solver::NotFound ) {
      solver->sendDataToNeighbour(toRank,srcCellDescriptionIndex,element,src,dest,x,level);
    } else {
      solver->sendEmptyDataToNeighbour(toRank,x,level);
    }
  }

  if ( sendMetadata ){
    exahype::sendNeighbourCommunicationMetadata(
        toRank,srcCellDescriptionIndex,src,dest,x,level);
  }
}


void exahype::Vertex::sendToNeighbour(
    int toRank,
    bool isLastIterationOfBatchOrNoBatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int level) const {
  if ( hasToCommunicate(level) ) {
    dfor2(dest)
      dfor2(src)
      if ( hasToSendMetadata(toRank,src,dest) ) {
        //#ifdef Asserts
        //logInfo("prepareSendToNeighbour(...)","to rank "<<toRank <<" vertex="<<x.toString()<<" src="<<src.toString()<<" dest="<<dest.toString());
        //#endif

        bool sendNoMetadata =
            (exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps &&
            exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps)
            ||
            (exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps
            && !isLastIterationOfBatchOrNoBatch);

        if ( hasToSendDataToNeighbour(src,dest) ) {
          sendSolverDataToNeighbour(
              toRank,!sendNoMetadata,src,dest,
              getCellDescriptionsIndex()[srcScalar],
              getCellDescriptionsIndex()[destScalar],
              x,level);
        } else {
          sendEmptySolverDataToNeighbour(
              toRank,!sendNoMetadata,src,dest,x,level);
        }
      }
      enddforx
    enddforx
  }
}

void exahype::Vertex::dropNeighbourData(
    const int fromRank,
    const int srcCellDescriptionIndex,
    const int destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int level) const {
  for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
    auto* solver = solvers::RegisteredSolvers[solverNumber];
    solver->dropNeighbourData(fromRank,src,dest,x,level);
  }
}

void exahype::Vertex::mergeWithNeighbourData(
        const int fromRank,
        const exahype::MetadataHeap::HeapEntries& receivedMetadata,
        const bool mergeWithNeighbourMetadata,
        const int srcCellDescriptionIndex,
        const int destCellDescriptionIndex,
        const tarch::la::Vector<DIMENSIONS,int>& src,
        const tarch::la::Vector<DIMENSIONS,int>& dest,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int level) const {
  for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
    auto* solver = solvers::RegisteredSolvers[solverNumber];
    const int element = solver->tryGetElement(destCellDescriptionIndex,solverNumber);
    if ( element!=exahype::solvers::Solver::NotFound ) {
      logDebug(
        "mergeWithNeighbour(...)", "receive data for solver " << solverNumber << " from " <<
        fromRank << " at vertex x=" << x << ", level=" << level <<
        ", src=" << src << ", dest=" << dest);

      solver->mergeWithNeighbourData(
        fromRank,
        destCellDescriptionIndex,element,src,dest,
        x,level);

      if ( mergeWithNeighbourMetadata ) {
          assertionEquals(receivedMetadata.size(),exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

          const int offset = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;
          exahype::MetadataHeap::HeapEntries metadataPortion(
              receivedMetadata.begin()+offset,
              receivedMetadata.begin()+offset+exahype::NeighbourCommunicationMetadataPerSolver);

          solver->mergeWithNeighbourMetadata(
                metadataPortion,
                src, dest,
                destCellDescriptionIndex,element);
      }
    } else {
      solver->dropNeighbourData(fromRank,src,dest,x,level);
    }
  }
}

void exahype::Vertex::receiveNeighbourData(
    int fromRank,
    bool mergeWithReceivedData,
    bool isFirstIterationOfBatchOrNoBatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    int level) const {
  if ( hasToCommunicate(level) ) {
    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest; // "invert" points
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

        const int destScalar = TWO_POWER_D - myDestScalar - 1; // "invert" point indices
        const int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;

        if ( hasToReceiveMetadata(fromRank,src,dest) ) {
          //#ifdef Asserts
          //logInfo("receiveNeighbourData(...)","from rank "<<fromRank <<" vertex="<<x.toString()<<" src="<<src.toString()<<" dest="<<dest.toString());
          //#endif

          bool receiveNoMetadata =
              (exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps &&
              exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps)
              ||
              (exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps
              && !isFirstIterationOfBatchOrNoBatch);

          MetadataHeap::HeapEntries receivedMetadata;
          receivedMetadata.clear();
          if ( !receiveNoMetadata ) {
            exahype::receiveNeighbourCommunicationMetadata(
                receivedMetadata, fromRank, x, level);
          }

          if( hasToMergeWithNeighbourData(src,dest) ) {
            if ( mergeWithReceivedData ) {
              mergeWithNeighbourData(
                  fromRank,
                  receivedMetadata,
                  !receiveNoMetadata,
                  getCellDescriptionsIndex()[srcScalar],
                  getCellDescriptionsIndex()[destScalar],
                  src,dest,
                  x,level);
            } else {
              dropNeighbourData(
                  fromRank,
                  getCellDescriptionsIndex()[srcScalar],
                  getCellDescriptionsIndex()[destScalar],
                  src,dest,
                  x,level);
            }
          } else {
            dropNeighbourData(
                fromRank,
                getCellDescriptionsIndex()[srcScalar],
                getCellDescriptionsIndex()[destScalar],
                src,dest,
                x,level);
          }
        }
      enddforx
    enddforx
  }
}
#endif
