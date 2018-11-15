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

constexpr int exahype::Vertex::pos1Scalar[2*(DIMENSIONS-1)*DIMENSIONS];
constexpr int exahype::Vertex::pos2Scalar[2*(DIMENSIONS-1)*DIMENSIONS];

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

bool exahype::Vertex::SpawnNeighbourMergeAsThread = false;

peano::MappingSpecification exahype::Vertex::getNeighbourMergeSpecification(const int level) {
  const int coarsestSolverLevel = solvers::Solver::getCoarsestMeshLevelOfAllSolvers();
  if ( std::abs(level)>=coarsestSolverLevel ) {
    return peano::MappingSpecification(
           peano::MappingSpecification::WholeTree,
           peano::MappingSpecification::AvoidFineGridRaces,false);
  } else {
    return peano::MappingSpecification(
          peano::MappingSpecification::Nop,
          peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
  }
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
    solvers::Solver::CellInfo cellInfo1(cellDescriptionsIndex1);
    solvers::Solver::CellInfo cellInfo2(cellDescriptionsIndex2);

    if ( !cellInfo1.empty() && !cellInfo2.empty() ) {
      for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
        if (solver->isMergingMetadata(section)) {
          switch ( solver->getType() ) {
            case solvers::Solver::Type::ADERDG:
              static_cast<solvers::ADERDGSolver*>(solver)->
                mergeNeighboursMetadata(solverNumber,cellInfo1,cellInfo2,pos1,pos2,x,h,checkThoroughly);
              break;
            case solvers::Solver::Type::LimitingADERDG:
              static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
                mergeNeighboursMetadata(solverNumber,cellInfo1,cellInfo2,pos1,pos2,x,h,checkThoroughly);
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

  for (int i=0; i<2*(DIMENSIONS-1)*(DIMENSIONS); i++) {
    mergeOnlyNeighboursMetadataLoopBody(pos1Scalar[i],pos2Scalar[i],_vertexData.getCellDescriptionsIndex(pos1Scalar[i]),_vertexData.getCellDescriptionsIndex(pos2Scalar[i]),section,x,h,checkThoroughly);
  }
}

void exahype::Vertex::validateNeighbourhood(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    solvers::Solver::CellInfo* (&cellInfos)[TWO_POWER_D],
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2) {
  tarch::la::Vector<DIMENSIONS,int> posCell  = pos1;
  tarch::la::Vector<DIMENSIONS,int> posEmpty = pos2;
  if ( cellDescriptionsIndex2 >= 0 ) {
    posCell  = pos2;
    posEmpty = pos1;
  }
  const int posCellScalar = peano::utils::dLinearised(posCell,2);
  solvers::Solver::CellInfo& cellInfo = *cellInfos[posCellScalar];
  solvers::Solver::BoundaryFaceInfo face(posCell,posEmpty);

  // ADER-DG
  for (auto& p : cellInfo._ADERDGCellDescriptions) {
    if (
        (p.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell ||
         p.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Ancestor)
        &&
        cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
    ) {
      logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<face._faceIndex<<" next to empty cell: cell="<<p.toString());
      std::terminate();
    }
  }
  // FV
  for (auto& p : cellInfo._FiniteVolumesCellDescriptions) {
    if (
        p.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Cell
        &&
        cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex
    ) {
      logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<face._faceIndex<<" next to empty cell: cell="<<p.toString());
      std::terminate();
    }
  }
}

void exahype::Vertex::mergeWithBoundaryData(
    solvers::Solver::CellInfo& cellInfo,
    const tarch::la::Vector<DIMENSIONS,int>& posCell,
    const tarch::la::Vector<DIMENSIONS,int>& posBoundary,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  solvers::Solver::BoundaryFaceInfo face(posCell,posBoundary);

  if ( !cellInfo.empty() && hasToMergeAtFace(cellInfo,face._faceIndex,true/*prefetchADERDGFace*/) ) {
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      switch ( solver->getType() ) {
        case solvers::Solver::Type::ADERDG:
          static_cast<solvers::ADERDGSolver*>(solver)->
            mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary);
          break;
        case solvers::Solver::Type::LimitingADERDG:
          static_cast<solvers::LimitingADERDGSolver*>(solver)->
            mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary,false);
          break;
        case solvers::Solver::Type::FiniteVolumes:
          static_cast<solvers::FiniteVolumesSolver*>(solver)->
            mergeWithBoundaryData(solverNumber,cellInfo,posCell,posBoundary);
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

bool exahype::Vertex::hasToMergeAtFace(
    solvers::Solver::CellInfo& cellInfo,
    const int                  faceIndex,
    const bool                 prefetchADERDGFaceData) {
  bool merge = false;
  for (auto& p : cellInfo._ADERDGCellDescriptions) {
    merge = !p.getNeighbourMergePerformed(faceIndex);
    p.setNeighbourMergePerformed(faceIndex,static_cast<char>(true));
    if ( merge && prefetchADERDGFaceData ) { solvers::ADERDGSolver::prefetchFaceData(p,faceIndex); }
  }
  for (auto& p : cellInfo._FiniteVolumesCellDescriptions) {
    merge = !p.getNeighbourMergePerformed(faceIndex);
    p.setNeighbourMergePerformed(faceIndex,static_cast<char>(true));
  }
  return merge;
}

void exahype::Vertex::mergeNeighboursDataAndMetadata(
    solvers::Solver::CellInfo& cellInfo1,
    solvers::Solver::CellInfo& cellInfo2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  if ( !cellInfo1.empty() && !cellInfo2.empty() ) {
    solvers::Solver::InterfaceInfo face(pos1,pos2);
    const bool mergeWithCell1 = hasToMergeAtFace(cellInfo1,face._faceIndex1,true/*prefetchADERDGFace*/);
    const bool mergeWithCell2 = hasToMergeAtFace(cellInfo2,face._faceIndex2,true/*prefetchADERDGFace*/);
    assertion2(mergeWithCell1==mergeWithCell2,mergeWithCell1,mergeWithCell2);
    if ( mergeWithCell1 && mergeWithCell2 ) {
      for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

        switch ( solver->getType() ) {
          case solvers::Solver::Type::ADERDG:
            static_cast<solvers::ADERDGSolver*>(solver)->
            mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
            static_cast<solvers::ADERDGSolver*>(solver)->
                mergeNeighboursMetadata(solverNumber,cellInfo1,cellInfo2,pos1,pos2,x,h,false);
            break;
          case solvers::Solver::Type::LimitingADERDG:
            static_cast<solvers::LimitingADERDGSolver*>(solver)->
            mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2,false);
            static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
                mergeNeighboursMetadata(solverNumber,cellInfo1,cellInfo2,pos1,pos2,x,h,false);
            break;
          case solvers::Solver::Type::FiniteVolumes:
            static_cast<solvers::FiniteVolumesSolver*>(solver)->
            mergeNeighboursData(solverNumber,cellInfo1,cellInfo2,pos1,pos2);
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

tarch::la::Vector<DIMENSIONS,int> exahype::Vertex::delineariseIndex2(int index) {
  #if DIMENSIONS==2
  switch( index ) {                                //x,y
    case 0: return tarch::la::Vector<DIMENSIONS,int>(0,0);
    case 1: return tarch::la::Vector<DIMENSIONS,int>(1,0);
    case 2: return tarch::la::Vector<DIMENSIONS,int>(0,1);
    case 3: return tarch::la::Vector<DIMENSIONS,int>(1,1);
    default:
      logError("delineariseIndex2(index)","index must be in range [0,3]!");
      return tarch::la::Vector<DIMENSIONS,int>(-1,-1);
  }
  #elif DIMENSIONS==3
  switch( index ) {                                //x,y,z
    case 0: return tarch::la::Vector<DIMENSIONS,int>(0,0,0);
    case 1: return tarch::la::Vector<DIMENSIONS,int>(1,0,0);
    case 2: return tarch::la::Vector<DIMENSIONS,int>(0,1,0);
    case 3: return tarch::la::Vector<DIMENSIONS,int>(1,1,0);
    case 4: return tarch::la::Vector<DIMENSIONS,int>(0,0,1);
    case 5: return tarch::la::Vector<DIMENSIONS,int>(1,0,1);
    case 6: return tarch::la::Vector<DIMENSIONS,int>(0,1,1);
    case 7: return tarch::la::Vector<DIMENSIONS,int>(1,1,1);
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
    solvers::Solver::CellInfo* (&cellInfos)[TWO_POWER_D],
    const tarch::la::Vector<DIMENSIONS, double> x,
    const tarch::la::Vector<DIMENSIONS, double> h) {
  #if defined(Asserts) || defined (ValidateNeighbourHoodDuringNeighbourMerge)
  constexpr bool validate = true;
  #else
  constexpr bool validate = false;
  #endif
  assertion2(cellDescriptionsIndex1 < 1 || cellDescriptionsIndex1 != cellDescriptionsIndex2,cellDescriptionsIndex1,cellDescriptionsIndex2);

  const tarch::la::Vector<DIMENSIONS,int> pos1 = delineariseIndex2(pos1Scalar);
  const tarch::la::Vector<DIMENSIONS,int> pos2 = delineariseIndex2(pos2Scalar);
  assertion(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1));

  bool validIndex1 = cellDescriptionsIndex1 >= 0;
  bool validIndex2 = cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || solvers::ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex1));
  assertion(cellDescriptionsIndex2 < 0 || solvers::ADERDGSolver::isValidCellDescriptionIndex(cellDescriptionsIndex2));

  if ( validIndex1 && validIndex2 ) {
    mergeNeighboursDataAndMetadata(*cellInfos[pos1Scalar],*cellInfos[pos2Scalar],pos1,pos2,x,h);
  }
  else if ( validIndex1 && cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex ) {
    mergeWithBoundaryData(*cellInfos[pos1Scalar],pos1,pos2,x,h);
  }
  else if ( validIndex2 && cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex ) {
    mergeWithBoundaryData(*cellInfos[pos2Scalar],pos2,pos1,x,h);
  }
  else if  (
      validate
      &&
      validIndex1 != validIndex2
      &&
      cellDescriptionsIndex1!=multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex&&
      cellDescriptionsIndex2!=multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex
  ) {
    validateNeighbourhood(cellDescriptionsIndex1,cellDescriptionsIndex2,cellInfos,pos1,pos2);
  }
}

void exahype::Vertex::mergeNeighbours(
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  if ( tarch::la::allSmallerEquals(h,exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers()) ) {
    // collect the cell descriptions once per vertex instead of once per face;
    // 8 instead of 12x2 lookups in 3D, 4 lookups instead of 8 in 2D.
    solvers::Solver::CellInfo* cellInfos[TWO_POWER_D] = { nullptr };
    for ( unsigned int i=0; i<TWO_POWER_D; i++ ) {
      if ( _vertexData.getCellDescriptionsIndex(i)>=0 ) {
        cellInfos[i] = new solvers::Solver::CellInfo(_vertexData.getCellDescriptionsIndex(i));
      } else {
        cellInfos[i] = nullptr;
      }
    }

    if ( SpawnNeighbourMergeAsThread ) {
      #if DIMENSIONS==2
      peano::datatraversal::TaskSet runParallelTasks(
      [&]() -> bool { mergeNeighboursLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),cellInfos,x,h); return false; },
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      peano::datatraversal::TaskSet::TaskType::IsTaskAndRunImmediately,
      true);
      #elif DIMENSIONS==3
      peano::datatraversal::TaskSet runParallelTasks(
      [&]() -> bool { mergeNeighboursLoopBody(0,1,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(1),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(0,2,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(2),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(0,4,_vertexData.getCellDescriptionsIndex(0),_vertexData.getCellDescriptionsIndex(4),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(1,3,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(3),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(1,5,_vertexData.getCellDescriptionsIndex(1),_vertexData.getCellDescriptionsIndex(5),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(2,3,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(3),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(2,6,_vertexData.getCellDescriptionsIndex(2),_vertexData.getCellDescriptionsIndex(6),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(3,7,_vertexData.getCellDescriptionsIndex(3),_vertexData.getCellDescriptionsIndex(7),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(4,5,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(5),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(4,6,_vertexData.getCellDescriptionsIndex(4),_vertexData.getCellDescriptionsIndex(6),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(5,7,_vertexData.getCellDescriptionsIndex(5),_vertexData.getCellDescriptionsIndex(7),cellInfos,x,h); return false; },
      [&]() -> bool { mergeNeighboursLoopBody(6,7,_vertexData.getCellDescriptionsIndex(6),_vertexData.getCellDescriptionsIndex(7),cellInfos,x,h); return false; },
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
    } else {
      for (int i=0; i<2*(DIMENSIONS-1)*(DIMENSIONS); i++) {
        mergeNeighboursLoopBody(pos1Scalar[i],pos2Scalar[i],_vertexData.getCellDescriptionsIndex(pos1Scalar[i]),_vertexData.getCellDescriptionsIndex(pos2Scalar[i]),cellInfos,x,h);
      }
    }

    // free the helper objects
    for ( unsigned int i=0; i<TWO_POWER_D; i++ ) {
      if ( cellInfos[i]!=nullptr ) {
        delete cellInfos[i];
        cellInfos[i]=nullptr;
      }
    }
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
  const int                                 toRank,
  const int                                 srcScalar,
  const int                                 destScalar,
  const tarch::la::Vector<TWO_POWER_D,int>& adjacentRanks) {
  return adjacentRanks[destScalar] == toRank
         &&
         adjacentRanks[destScalar] != tarch::parallel::Node::getGlobalMasterRank() &&
         adjacentRanks[srcScalar]  != tarch::parallel::Node::getGlobalMasterRank()
         &&
         (   // Send also when a fork/join was triggered for the current rank
             adjacentRanks[srcScalar]   == tarch::parallel::Node::getInstance().getRank() ||
             State::isForkTriggeredForRank(adjacentRanks[srcScalar]) ||
             State::isJoinTriggeredForRank(adjacentRanks[srcScalar])
         );
}

bool exahype::Vertex::compareGeometryInformationOfCellDescriptionsAndVertex(
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionsIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) {
  if ( exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex) ) {
    solvers::Solver::CellInfo cellInfo(srcCellDescriptionsIndex);

    if ( !cellInfo.empty() ) {
      solvers::Solver::BoundaryFaceInfo face(src,dest);

      bool geometryInformationMatches = true;

      for ( auto& patch : cellInfo._ADERDGCellDescriptions ) {
        tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
            exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),face._direction,face._orientation);
        tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
            exahype::Vertex::computeFaceBarycentre(x,h,face._direction,src);
        geometryInformationMatches &= equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
      }
      for ( auto& patch : cellInfo._FiniteVolumesCellDescriptions ) {
        tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
            exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),face._direction,face._orientation);
        tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
            exahype::Vertex::computeFaceBarycentre(x,h,face._direction,src);
        geometryInformationMatches &= equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
      }
      return geometryInformationMatches;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

void exahype::Vertex::sendOnlyMetadataToNeighbourLoopBody(
    const int                                    toRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const int                                    srcCellDescriptionsIndex,
    const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const bool                                   checkThoroughly) {
  if ( hasToSendMetadata(toRank,srcScalar,destScalar,adjacentRanks) ) {
    const tarch::la::Vector<DIMENSIONS,int> src  = delineariseIndex2(srcScalar);
    const tarch::la::Vector<DIMENSIONS,int> dest = delineariseIndex2(destScalar);
    assertion(tarch::la::countEqualEntries(src,dest)==(DIMENSIONS-1));

    if ( !checkThoroughly || compareGeometryInformationOfCellDescriptionsAndVertex(src,dest,srcCellDescriptionsIndex,x,h) ) {
      exahype::sendNeighbourCommunicationMetadata(
          toRank,srcCellDescriptionsIndex,src,dest,x,level);
    } else {
      exahype::sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(
          toRank,x,level);
    }
  }
}

void exahype::Vertex::sendOnlyMetadataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const bool                                   checkThoroughly) const {
  if ( hasToCommunicate(level) ) {
    for (unsigned int i=0; i<2*(DIMENSIONS-1)*(DIMENSIONS); i++) {
      sendOnlyMetadataToNeighbourLoopBody(toRank,pos1Scalar[i],pos2Scalar[i],_vertexData.getCellDescriptionsIndex(pos1Scalar[i]),getAdjacentRanks(),x,h,level,checkThoroughly);
      sendOnlyMetadataToNeighbourLoopBody(toRank,pos2Scalar[i],pos1Scalar[i],_vertexData.getCellDescriptionsIndex(pos2Scalar[i]),getAdjacentRanks(),x,h,level,checkThoroughly); // swap
    }
  }
}


bool exahype::Vertex::hasToReceiveMetadata(
    const int                                  fromRank,
    const int                                  srcScalar,
    const int                                  destScalar,
    const tarch::la::Vector<TWO_POWER_D, int>& adjacentRanks) {
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
      );
}


void exahype::Vertex::mergeOnlyWithNeighbourMetadataLoopBody(
    const int                                    fromRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const exahype::State::AlgorithmSection&      section,
    const bool                                   checkThoroughly) {
  if ( hasToReceiveMetadata(fromRank,srcScalar,destScalar,adjacentRanks) ) {
    logDebug("mergeOnlyWithNeighbourMetadata(...)","from rank="<<fromRank<<",x="<<x.toString()<<",level="<<level<<",adjacentRanks="<<adjacentRanks);
    const tarch::la::Vector<DIMENSIONS,int> src  = delineariseIndex2(srcScalar);
    const tarch::la::Vector<DIMENSIONS,int> dest = delineariseIndex2(destScalar);
    assertion(tarch::la::countEqualEntries(dest, src) == (DIMENSIONS-1));

    exahype::MetadataHeap::HeapEntries receivedMetadata;
    receivedMetadata.clear();
    exahype::receiveNeighbourCommunicationMetadata(receivedMetadata,fromRank, x, level);
    assertionEquals(receivedMetadata.size(),exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

    bool validIndex = destCellDescriptionIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));
    if ( validIndex ) {
      solvers::Solver::CellInfo cellInfo(destCellDescriptionIndex);

      for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
        auto* solver = solvers::RegisteredSolvers[solverNumber];
        if ( solver->isMergingMetadata(section) ) {
          const int offset = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;
          exahype::MetadataHeap::HeapEntries metadataPortion(
              receivedMetadata.begin()+offset,
              receivedMetadata.begin()+offset+exahype::NeighbourCommunicationMetadataPerSolver);

          switch ( solver->getType() ) {
            case solvers::Solver::Type::ADERDG:
              static_cast<solvers::ADERDGSolver*>(solver)->
                mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,src,dest);
              break;
            case solvers::Solver::Type::LimitingADERDG:
              static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
                mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,src,dest);
              break;
            case solvers::Solver::Type::FiniteVolumes:
              // do nothing
              break;
            default:
              assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              logError("mergeOnlyWithNeighbourMetadataLoopBody(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              std::abort();
              break;
          }
        }
      }
    }
  }
}

void exahype::Vertex::mergeOnlyWithNeighbourMetadata(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int                                    level,
    const exahype::State::AlgorithmSection&      section,
    const bool                                   checkThoroughly) const {
  if ( hasToCommunicate(level) ) {
    for (unsigned int i = 2*(DIMENSIONS-1)*(DIMENSIONS); i-- > 0;) { // dest and src is swapped & order is swapped
      mergeOnlyWithNeighbourMetadataLoopBody(fromRank,pos1Scalar[i],pos2Scalar[i],_vertexData.getCellDescriptionsIndex(pos2Scalar[i]),getAdjacentRanks(),x,h,level,section,checkThoroughly);
      mergeOnlyWithNeighbourMetadataLoopBody(fromRank,pos2Scalar[i],pos1Scalar[i],_vertexData.getCellDescriptionsIndex(pos1Scalar[i]),getAdjacentRanks(),x,h,level,section,checkThoroughly);
    }
  }
}

void exahype::Vertex::dropNeighbourMetadata(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  if ( hasToCommunicate(level) ) {
    for (unsigned int i = 2*(DIMENSIONS-1)*(DIMENSIONS); i-- > 0;) { // dest and src is swapped & order is swapped
      if ( hasToReceiveMetadata(fromRank,pos1Scalar[i],pos2Scalar[i],getAdjacentRanks()) ) { // order of conditionals important; do not combine statements
        logDebug("dropNeighbourMetadata(...)","from rank="<<fromRank<<",x="<<x.toString()<<",level="<<level<<",adjacentRanks="<<getAdjacentRanks());
        MetadataHeap::getInstance().receiveData(fromRank,x,level,peano::heap::MessageType::NeighbourCommunication);
      }
      if ( hasToReceiveMetadata(fromRank,pos2Scalar[i],pos1Scalar[i],getAdjacentRanks()) ) {
        logDebug("dropNeighbourMetadata(...)","from rank="<<fromRank<<",x="<<x.toString()<<",level="<<level<<",adjacentRanks="<<getAdjacentRanks());
        MetadataHeap::getInstance().receiveData(fromRank,x,level,peano::heap::MessageType::NeighbourCommunication);
      }
    }
  }
}

bool exahype::Vertex::hasToSendToNeighbourNow(
    solvers::Solver::CellInfo&         cellInfo,
    solvers::Solver::BoundaryFaceInfo& face) {
  bool result = true;
  for (auto& p : cellInfo._ADERDGCellDescriptions) {
    result |= hasToSendToNeighbourNow(p,face); // side effects
  }
  for (auto& p : cellInfo._FiniteVolumesCellDescriptions) {
    result |= hasToSendToNeighbourNow(p,face); // side effects
  }
  return result;
}

void exahype::Vertex::sendToNeighbourLoopBody(
  const int                                    toRank,
  const int                                    srcScalar,
  const int                                    destScalar,
  const int                                    srcCellDescriptionsIndex,
  const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
  const bool                                   isLastIterationOfBatchOrNoBatch,
  const tarch::la::Vector<DIMENSIONS, double>& x,
  const int                                    level) {
  if ( hasToSendMetadata(toRank,srcScalar,destScalar,adjacentRanks) ) {
    const tarch::la::Vector<DIMENSIONS,int> src  = delineariseIndex2(srcScalar);
    const tarch::la::Vector<DIMENSIONS,int> dest = delineariseIndex2(destScalar);
    assertion(tarch::la::countEqualEntries(src,dest)==(DIMENSIONS-1));
    logDebug("sendToNeighbourLoopBody(...)","to rank="<<toRank<<",x="<<x.toString()<<",src="<<src.toString()<<",dest="<<dest.toString());

    bool validIndex = srcCellDescriptionsIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex));

    if ( validIndex ) {
      solvers::Solver::CellInfo cellInfo(srcCellDescriptionsIndex);
      solvers::Solver::BoundaryFaceInfo face(src,dest);

      if ( hasToSendToNeighbourNow(cellInfo,face) ) {
        for ( unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber ) {
          auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
          switch ( solver->getType() ) {
            case solvers::Solver::Type::ADERDG:
              static_cast<solvers::ADERDGSolver*>(solver)->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
              break;
            case solvers::Solver::Type::LimitingADERDG:
              static_cast<solvers::LimitingADERDGSolver*>(solver)->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
              break;
            case solvers::Solver::Type::FiniteVolumes:
              static_cast<solvers::FiniteVolumesSolver*>(solver)->sendDataToNeighbour(toRank,solverNumber,cellInfo,src,dest,x,level);
              break;
            default:
              assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              logError("mergeOnlyWithNeighbourMetadataLoopBody(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
              std::abort();
              break;
          }
        }
      }
    }
    // metadata is sent and received as block
    bool sendNoMetadata = /* negative */
        (exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps &&
            exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps)
            ||
            (exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps
                && !isLastIterationOfBatchOrNoBatch);
    if ( !sendNoMetadata && validIndex ){
      sendNeighbourCommunicationMetadata(
          toRank,srcCellDescriptionsIndex,src,dest,x,level);
    } else if ( !sendNoMetadata ) {
      sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(toRank,x,level);
    }
  }
}

void exahype::Vertex::sendToNeighbour(
    int                                          toRank,
    bool                                         isLastIterationOfBatchOrNoBatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  if ( hasToCommunicate(level) ) {
    for (unsigned int i=0; i<2*(DIMENSIONS-1)*(DIMENSIONS); i++) {
      sendToNeighbourLoopBody(toRank,pos1Scalar[i],pos2Scalar[i],_vertexData.getCellDescriptionsIndex(pos1Scalar[i]),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,level);
      sendToNeighbourLoopBody(toRank,pos2Scalar[i],pos1Scalar[i],_vertexData.getCellDescriptionsIndex(pos2Scalar[i]),getAdjacentRanks(),isLastIterationOfBatchOrNoBatch,x,level); // swap
    }
  }
}

bool exahype::Vertex::hasToReceiveFromNeighbourNow(
    solvers::Solver::CellInfo&         cellInfo,
    solvers::Solver::BoundaryFaceInfo& face,
    const bool prefetchADERDGFaceData) {
  bool result = true;
  for (auto& p : cellInfo._ADERDGCellDescriptions) {
    result |= hasToReceiveFromNeighbourNow(p,face); // side effects
    if ( prefetchADERDGFaceData ) { solvers::ADERDGSolver::prefetchFaceData(p,face._faceIndex); }
  }
  for (auto& p : cellInfo._FiniteVolumesCellDescriptions) {
    result |= hasToReceiveFromNeighbourNow(p,face); // side effects
  }
  return result;
}

void exahype::Vertex::receiveNeighbourDataLoopBody(
    const int                                    fromRank,
    const int                                    srcScalar,
    const int                                    destScalar,
    const int                                    destCellDescriptionIndex,
    const bool                                   mergeWithReceivedData,
    const bool                                   receiveNeighbourMetadata,
    const tarch::la::Vector<TWO_POWER_D, int>&   adjacentRanks,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  if ( hasToReceiveMetadata(fromRank,srcScalar,destScalar,adjacentRanks) ) {
    logDebug("mergeOnlyWithNeighbourMetadata(...)","from rank="<<fromRank<<",x="<<x.toString()<<",level="<<level<<",adjacentRanks="<<adjacentRanks);
    const tarch::la::Vector<DIMENSIONS,int> src  = delineariseIndex2(srcScalar);
    const tarch::la::Vector<DIMENSIONS,int> dest = delineariseIndex2(destScalar);
    assertion(tarch::la::countEqualEntries(dest, src) == (DIMENSIONS-1));

    exahype::MetadataHeap::HeapEntries receivedMetadata;
    receivedMetadata.clear();
    if ( receiveNeighbourMetadata ) {
      exahype::receiveNeighbourCommunicationMetadata(receivedMetadata,fromRank,x,level);
      assertionEquals(receivedMetadata.size(),exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());
    }

    bool validIndex = destCellDescriptionIndex >= 0;
    assertion( !validIndex || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));
    if ( validIndex ) {
      solvers::Solver::CellInfo cellInfo(destCellDescriptionIndex);
      solvers::Solver::BoundaryFaceInfo face(dest,src); // dest and src are swapped

      if ( hasToReceiveFromNeighbourNow(cellInfo,face,true/*prefetchADERDGFaceData*/) ) {
        for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
          auto* solver = solvers::RegisteredSolvers[solverNumber];

          const int offset = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;
          exahype::MetadataHeap::HeapEntries metadataPortion(
              receivedMetadata.begin()+offset,
              receivedMetadata.begin()+offset+exahype::NeighbourCommunicationMetadataPerSolver);

          switch ( solver->getType() ) {
          case solvers::Solver::Type::ADERDG:
            if ( mergeWithReceivedData ) {
              static_cast<solvers::ADERDGSolver*>(solver)->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
              if ( receiveNeighbourMetadata ) {
                static_cast<solvers::ADERDGSolver*>(solver)->mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,src,dest);
              }
            } else {
              static_cast<solvers::ADERDGSolver*>(solver)->dropNeighbourData(fromRank,x,level);
            }
            break;
          case solvers::Solver::Type::LimitingADERDG:
            if ( mergeWithReceivedData ) {
              static_cast<solvers::LimitingADERDGSolver*>(solver)->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
              if ( receiveNeighbourMetadata ) {
                static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver()->
                    mergeWithNeighbourMetadata(solverNumber,cellInfo,metadataPortion,src,dest);
              }
            } else {
              static_cast<solvers::LimitingADERDGSolver*>(solver)->dropNeighbourData(fromRank,x,level);
            }
            break;
          case solvers::Solver::Type::FiniteVolumes:
            if ( mergeWithReceivedData ) {
              static_cast<solvers::FiniteVolumesSolver*>(solver)->mergeWithNeighbourData(fromRank,solverNumber,cellInfo,src,dest,x,level);
            } else {
              static_cast<solvers::FiniteVolumesSolver*>(solver)->dropNeighbourData(fromRank,x,level);
            }
            break;
          default:
            assertionMsg(false,"Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
            logError("mergeOnlyWithNeighbourMetadataLoopBody(...)","Unrecognised solver type: "<<solvers::Solver::toString(solver->getType()));
            std::abort();
            break;
          }
        }
      }
    }
  }
}

void exahype::Vertex::receiveNeighbourData(
    const int                                    fromRank,
    const bool                                   mergeWithReceivedData,
    const bool                                   isFirstIterationOfBatchOrNoBatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  if ( hasToCommunicate(level) ) {
    bool receiveNoMetadata =
              (exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps &&
              exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps)
              ||
              (exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps
              && !isFirstIterationOfBatchOrNoBatch);
    // todo template these loops
    for (unsigned int i = 2*(DIMENSIONS-1)*(DIMENSIONS); i-- > 0;) { // dest and src is swapped & order is swapped
      receiveNeighbourDataLoopBody(fromRank,pos1Scalar[i],pos2Scalar[i],_vertexData.getCellDescriptionsIndex(pos2Scalar[i]),mergeWithReceivedData,!receiveNoMetadata,getAdjacentRanks(),x,level);
      receiveNeighbourDataLoopBody(fromRank,pos2Scalar[i],pos1Scalar[i],_vertexData.getCellDescriptionsIndex(pos1Scalar[i]),mergeWithReceivedData,!receiveNoMetadata,getAdjacentRanks(),x,level);
    }
  }
}
#endif
