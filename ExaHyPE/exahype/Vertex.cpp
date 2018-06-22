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

#include "peano/utils/Loop.h"
#include "peano/grid/Checkpoint.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/State.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

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
  double scaling =
      std::max(
          1.0, std::max( tarch::la::maxAbs(lhs), tarch::la::maxAbs(rhs) )
  );
  return tarch::la::equals( lhs, rhs, scaling*tarch::la::NUMERICAL_ZERO_DIFFERENCE );
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
    tarch::la::Vector<DIMENSIONS,double> cellOffset(get);
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

void exahype::Vertex::mergeOnlyNeighboursMetadata(
    const exahype::State::AlgorithmSection& section,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  assertion(!isHangingNode());
  assertion(isInside() || isBoundary());

  dfor2(pos1)
    dfor2(pos2)
      if ( determineInterfaceType(pos1,pos1Scalar,pos2,pos2Scalar,x,h,false)==InterfaceType::Interior ) { // Implies that we have two valid indices on the correct level
        for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
          auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
          if (solver->isMergingMetadata(section)) {
            const int cellDescriptionsIndex1 = _vertexData.getCellDescriptionsIndex(pos1Scalar);
            const int cellDescriptionsIndex2 = _vertexData.getCellDescriptionsIndex(pos2Scalar);
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            if (element2>=0 && element1>=0) {
              solver->mergeNeighboursMetadata(
                  cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
            }
          }
        }
      }
    enddforx
  enddforx
}

bool exahype::Vertex::hasToMergeNeighbours(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  assertion(!isHangingNode());

  if ( cellDescriptionsIndex1!=cellDescriptionsIndex2 ) { // happened once
    assertion1(cellDescriptionsIndex1!=cellDescriptionsIndex2,cellDescriptionsIndex1);
    assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1),
        cellDescriptionsIndex1);
    assertion1(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2),
        cellDescriptionsIndex2);

    const int direction    = tarch::la::equalsReturnIndex(pos1, pos2);
    const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
    const int orientation2 = 1-orientation1;

    const int faceIndex1 = 2*direction+orientation1;
    const int faceIndex2 = 2*direction+orientation2;

    bool mergeNeighbours =
        !exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1).empty() ||
        !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1).empty() ||
        !exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2).empty() ||
        !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2).empty();

    // cell 1
    tarch::la::Vector<DIMENSIONS,double> baryCentreFromPatch1;
    for (auto& p1 : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
      baryCentreFromPatch1 =
          exahype::Cell::computeFaceBarycentre(p1.getOffset(),p1.getSize(),direction,orientation1);
      mergeNeighbours &= !p1.getNeighbourMergePerformed(faceIndex1);

      // now we can set the flag
      p1.setNeighbourMergePerformed(faceIndex1,true);
    }
    for (auto& p1 : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
      baryCentreFromPatch1 =
          exahype::Cell::computeFaceBarycentre(p1.getOffset(),p1.getSize(),direction,orientation1);
      mergeNeighbours &= !p1.getNeighbourMergePerformed(faceIndex1);

      // now we can set the flag
      p1.setNeighbourMergePerformed(faceIndex1,true);
    }

    // cell 2
    tarch::la::Vector<DIMENSIONS,double> baryCentreFromPatch2;
    for (auto& p2 : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
      baryCentreFromPatch2 =
          exahype::Cell::computeFaceBarycentre(p2.getOffset(),p2.getSize(),direction,orientation2);
      mergeNeighbours &= !p2.getNeighbourMergePerformed(faceIndex2);
      // assertion(p2.getNeighbourMergePerformed(faceIndex2) || mergeNeighbours);
      // TODO(Dominic): This is not always true during the mesh refinement iterations with MPI turned on
      // and more than 2 ranks.

      // now we can set the flag
      p2.setNeighbourMergePerformed(faceIndex2,true);
    }
    for (auto& p2 : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
      baryCentreFromPatch2 =
          exahype::Cell::computeFaceBarycentre(p2.getOffset(),p2.getSize(),direction,orientation2);
      mergeNeighbours &= !p2.getNeighbourMergePerformed(faceIndex2);
      // assertion(p2.getNeighbourMergePerformed(faceIndex2) || mergeNeighbours);

      // now we can set the flag
      p2.setNeighbourMergePerformed(faceIndex2,true);
    }

    tarch::la::Vector<DIMENSIONS,double> baryCentreFromVertex =
              exahype::Vertex::computeFaceBarycentre(x,h,direction,pos2);

    mergeNeighbours &= // ensure the barycentres match
        equalUpToRelativeTolerance(baryCentreFromPatch1,baryCentreFromPatch2) &&
        equalUpToRelativeTolerance(baryCentreFromPatch1,baryCentreFromVertex);

    return mergeNeighbours;
  } else  {
    return false;
  }
}

bool exahype::Vertex::hasToMergeWithBoundaryData(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  const int direction    = tarch::la::equalsReturnIndex(pos1, pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;
  const int faceIndex1 = 2*direction+orientation1;
  const int faceIndex2 = 2*direction+orientation2;

  if ( cellDescriptionsIndex1 >= 0 ) {
    bool mergeWithBoundaryData =
        !exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1).empty() ||
        !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1).empty();

    tarch::la::Vector<DIMENSIONS,double> baryCentreFromPatch1;
    for (auto& p1 : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
      mergeWithBoundaryData &=
          !p1.getNeighbourMergePerformed(faceIndex1);
      baryCentreFromPatch1 = exahype::Cell::computeFaceBarycentre(p1.getOffset(),p1.getSize(),direction,orientation1);

      // now we can set the flag
      p1.setNeighbourMergePerformed(faceIndex1,true);
    }
    for (auto& p1 : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex1)) {
      mergeWithBoundaryData &=
          !p1.getNeighbourMergePerformed(faceIndex1);
      baryCentreFromPatch1 = exahype::Cell::computeFaceBarycentre(p1.getOffset(),p1.getSize(),direction,orientation1);

      // now we can set the flag
      p1.setNeighbourMergePerformed(faceIndex1,true);
    }

    tarch::la::Vector<DIMENSIONS,double> baryCentreFromVertex =
        exahype::Vertex::computeFaceBarycentre(x,h,direction,pos1);
    mergeWithBoundaryData &= equalUpToRelativeTolerance(baryCentreFromPatch1,baryCentreFromVertex);
    return mergeWithBoundaryData;
  } else if ( cellDescriptionsIndex2 >= 0 ) {
    bool mergeWithBoundaryData =
        !exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2).empty() ||
        !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2).empty();

    tarch::la::Vector<DIMENSIONS,double> baryCentreFromPatch2;
    for (auto& p2 : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
      mergeWithBoundaryData &=
          !p2.getNeighbourMergePerformed(faceIndex2);
      baryCentreFromPatch2 = exahype::Cell::computeFaceBarycentre(p2.getOffset(),p2.getSize(),direction,orientation2);

      // now we can set the flag
      p2.setNeighbourMergePerformed(faceIndex2,true);
    }
    for (auto& p2 : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex2)) {
      mergeWithBoundaryData &=
          !p2.getNeighbourMergePerformed(faceIndex2);
      baryCentreFromPatch2 = exahype::Cell::computeFaceBarycentre(p2.getOffset(),p2.getSize(),direction,orientation2);

      // now we can set the flag
      p2.setNeighbourMergePerformed(faceIndex2,true);
    }

    tarch::la::Vector<DIMENSIONS,double> baryCentreFromVertex =
        exahype::Vertex::computeFaceBarycentre(x,h,direction,pos2);
    mergeWithBoundaryData &= equalUpToRelativeTolerance(baryCentreFromPatch2,baryCentreFromVertex);
    return mergeWithBoundaryData;
  } else {
    return false;
  }
}

void exahype::Vertex::validateNeighbourhood(
    const int cellDescriptionsIndex1,
    const int cellDescriptionsIndex2,
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const tarch::la::Vector<DIMENSIONS,int>& pos2) const {
  const int direction    = tarch::la::equalsReturnIndex(pos1, pos2);
  const int orientation1 = (1 + pos2(direction) - pos1(direction))/2;
  const int orientation2 = 1-orientation1;

  const int faceIndex1 = 2*direction+orientation1;
  const int faceIndex2 = 2*direction+orientation2;

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
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<faceIndex1<<" next to empty cell: cell="<<p1.toString());
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
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex2<<" is at face="<<faceIndex2<<" next to empty cell: cell="<<p2.toString());
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
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex1<<" is at face="<<faceIndex1<<" next to empty cell: cell="<<p1.toString());
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
          logError("validateNeighbourhood(...)","cell at index="<<cellDescriptionsIndex2<<" is at face="<<faceIndex2<<" next to empty cell: cell="<<p2.toString());
          std::terminate();
        }
      }
    } break;
    }
  }
}

exahype::Vertex::InterfaceType exahype::Vertex::determineInterfaceType(
    const tarch::la::Vector<DIMENSIONS,int>& pos1,
    const int pos1Scalar,
    const tarch::la::Vector<DIMENSIONS,int>& pos2,
    const int pos2Scalar,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const bool validate) const {
  const int cellDescriptionsIndex1 = _vertexData.getCellDescriptionsIndex(pos1Scalar);
  const int cellDescriptionsIndex2 = _vertexData.getCellDescriptionsIndex(pos2Scalar);

  const bool isFace = tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1);

  bool validIndex1 =
      isFace &&
      cellDescriptionsIndex1 >= 0;
  assertion(cellDescriptionsIndex1 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex1));

  bool validIndex2 =
      isFace &&
      cellDescriptionsIndex2 >= 0;
  assertion(cellDescriptionsIndex2 < 0 || exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(cellDescriptionsIndex2));

  if (
      validIndex1 && validIndex2
      &&
      hasToMergeNeighbours(cellDescriptionsIndex1,cellDescriptionsIndex2,pos1,pos2,x,h)
  ) {
    return InterfaceType::Interior;
  } else if (
      ((validIndex1 && !validIndex2 &&
      cellDescriptionsIndex2==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex)
      ||
      (validIndex2 && !validIndex1 &&
      cellDescriptionsIndex1==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex))
      &&
      hasToMergeWithBoundaryData(cellDescriptionsIndex1,cellDescriptionsIndex2,pos1,pos2,x,h)
  ) {
    return InterfaceType::Boundary;
  } else if  (
      validate
      &&
      validIndex1 != validIndex2
      &&
      cellDescriptionsIndex1!=multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex&&
      cellDescriptionsIndex2!=multiscalelinkedcell::HangingVertexBookkeeper::RemoteAdjacencyIndex
  ) {
    validateNeighbourhood(cellDescriptionsIndex1,cellDescriptionsIndex2,pos1,pos2);
    return InterfaceType::None;
  } else {
    return InterfaceType::None;
  }
}

void exahype::Vertex::mergeNeighboursDataAndMetadata(
    const tarch::la::Vector<DIMENSIONS,int>&  pos1,
    const int pos1Scalar,
    const tarch::la::Vector<DIMENSIONS,int>&  pos2,
    const int pos2Scalar) const {
  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    const int cellDescriptionsIndex1 = _vertexData.getCellDescriptionsIndex(pos1Scalar);
    const int cellDescriptionsIndex2 = _vertexData.getCellDescriptionsIndex(pos2Scalar);
    const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
    const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
    if (element2>=0 && element1>=0) {
      solver->mergeNeighbours(
          cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
      solver->mergeNeighboursMetadata(
          cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
    }
  }
}

void exahype::Vertex::mergeWithBoundaryData(
    const tarch::la::Vector<DIMENSIONS,int>&  pos1,
    const int pos1Scalar,
    const tarch::la::Vector<DIMENSIONS,int>&  pos2,
    const int pos2Scalar) const {
  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    const int cellDescriptionsIndex1 = _vertexData.getCellDescriptionsIndex(pos1Scalar);
    const int cellDescriptionsIndex2 = _vertexData.getCellDescriptionsIndex(pos2Scalar);
    int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
    int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
    assertion4((element1==exahype::solvers::Solver::NotFound &&
                element2==exahype::solvers::Solver::NotFound)
               || (element1 >= 0 && element2==exahype::solvers::Solver::NotFound)
               || (element2 >= 0 && element1==exahype::solvers::Solver::NotFound),
               cellDescriptionsIndex1,cellDescriptionsIndex2,element1,element2);

    if (element1 >= 0) {
      solver->mergeWithBoundaryData(cellDescriptionsIndex1,element1,pos1,pos2);
    }
    if (element2 >= 0){
      solver->mergeWithBoundaryData(cellDescriptionsIndex2,element2,pos2,pos1);
    }
  }
}


void exahype::Vertex::mergeNeighbours(
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  if ( tarch::la::allSmallerEquals(h,exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers()) ) {
    dfor2(pos1)
      dfor2(pos2)
        InterfaceType interfaceType = determineInterfaceType(pos1,pos1Scalar,pos2,pos2Scalar,x,h,true);

        if ( interfaceType==InterfaceType::Interior ) { // Assumes that we have two valid indices
          mergeNeighboursDataAndMetadata(pos1,pos1Scalar,pos2,pos2Scalar);
        } else if ( interfaceType==InterfaceType::Boundary ) {
          mergeWithBoundaryData(pos1,pos1Scalar,pos2,pos2Scalar);
        }
      enddforx
    enddforx
  }
}

// PARALLEL

#if Parallel
bool exahype::Vertex::hasToCommunicate(
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  if (!isInside()) {
    return false;
  } else if (tarch::la::oneGreater(h,exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers())) {
    return  false;
  } else {
    return true;
  }
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
    const int direction   = tarch::la::equalsReturnIndex(src, dest);
    const int orientation = (1 + dest(direction) - src(direction))/2;

    if (!exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionIndex).empty()) {
      auto& patch = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionIndex)[0];
      tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
          exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),direction,orientation);
      tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
          exahype::Vertex::computeFaceBarycentre(x,h,direction,src);
      return equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
    }
    else if (!exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionIndex).empty()) {
      auto& patch = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionIndex)[0];
      tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch = // copy & paste from here
          exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),direction,orientation);
      tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
          exahype::Vertex::computeFaceBarycentre(x,h,direction,src);
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
  if ( hasToCommunicate(h) ) {
    tarch::la::Vector<TWO_POWER_D, int> adjacentADERDGCellDescriptionsIndices = getCellDescriptionsIndex();
    dfor2(dest)
      dfor2(src)
        if (hasToSendMetadata(toRank,src,dest)) {
          const int srcCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(srcScalar);
          if (hasToSendMetadataToNeighbour(src,dest,x,h)) {
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
      (   // Receive also when a fork/join is performed for the neighbour rank, i.e. such
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
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h) const {
  const int destScalar = peano::utils::dLinearisedWithoutLookup(dest,2);
  const int destCellDescriptionIndex = getCellDescriptionsIndex()[destScalar];

  bool mergeWithNeighbourMetadata =
      exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex) &&
      (!exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionIndex).empty() ||
      !exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionIndex).empty());

  if (mergeWithNeighbourMetadata) {
    const int direction   = tarch::la::equalsReturnIndex(src, dest);
    const int orientation = (1 + src(direction) - dest(direction))/2;

    if (!exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionIndex).empty()) {
      auto& patch = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionIndex)[0];
      tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch =
          exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),direction,orientation);
      tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
          exahype::Vertex::computeFaceBarycentre(x,h,direction,dest);
      return equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
    }
    else if (!exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionIndex).empty()) {
      auto& patch = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionIndex)[0];
      tarch::la::Vector<DIMENSIONS,double> barycentreFromPatch = // copy & paste from here
          exahype::Cell::computeFaceBarycentre(patch.getOffset(),patch.getSize(),direction,orientation);
      tarch::la::Vector<DIMENSIONS,double> barycentreFromVertex =
          exahype::Vertex::computeFaceBarycentre(x,h,direction,dest);
      return equalUpToRelativeTolerance(barycentreFromPatch,barycentreFromVertex);
    } else {
      return false;
    }
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
  if ( hasToCommunicate(h) ) {
    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;
        int destScalar = TWO_POWER_D - myDestScalar - 1;

        if (hasToReceiveMetadata(fromRank,src,dest)) {
          logDebug("mergeOnlyWithNeighbourMetadata(...)","[pre] rec. from rank "<<fromRank<<", x:"<<
                   x.toString() << ", level=" <<level << ", adjacentRanks: "
                   << getAdjacentRanks());

          const int receivedMetadataIndex =
              exahype::receiveNeighbourCommunicationMetadata(fromRank, x, level);
          exahype::MetadataHeap::HeapEntries& receivedMetadata =
              MetadataHeap::getInstance().getData(receivedMetadataIndex);
          assertionEquals(receivedMetadata.size(),
              exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

          if (hasToMergeWithNeighbourMetadata(src,dest,x,h)) {
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
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int level) const {
  if ( hasToCommunicate(h) ) {
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
    assertion(!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionsIndex) );

    const int direction   = tarch::la::equalsReturnIndex(src, dest);
    const int orientation = (1 + dest(direction) - src(direction))/2;
    const int faceIndex   = 2*direction+orientation;

    bool hasToSendDataToNeighbour = true;
    // ADER-DG
    for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(srcCellDescriptionsIndex)) {
      // decrement counter beforehand
      int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
      assertion2(newCounterValue>=0,newCounterValue,p.toString());
      assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
      p.setFaceDataExchangeCounter(faceIndex,newCounterValue);

      hasToSendDataToNeighbour &= p.getFaceDataExchangeCounter(faceIndex)==0;
    }

    // FV (copied from above)
    for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(srcCellDescriptionsIndex)) {
      // decrement counter beforehand
      int newCounterValue = p.getFaceDataExchangeCounter(faceIndex)-1;
      assertion2(newCounterValue>=0,newCounterValue,p.toString());
      assertion1(newCounterValue<TWO_POWER_D,newCounterValue);
      p.setFaceDataExchangeCounter(faceIndex,newCounterValue);

      hasToSendDataToNeighbour &= p.getFaceDataExchangeCounter(faceIndex)==0;
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

  if (
      !exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionsIndex)
      ||
      (exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty() &&
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionsIndex).empty())
  ) {
    return false; // !!! Make sure to consider all solver types here
  } else {
    const int direction   = tarch::la::equalsReturnIndex(src, dest);
    const int orientation = (1 + src(direction) - dest(direction))/2;
    const int faceIndex   = 2*direction+orientation;
  
    // ADER-DG
    bool mergeWithNeighbourData = true;
    for (auto& p : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) { // lookup is slow
      mergeWithNeighbourData &= p.getFaceDataExchangeCounter(faceIndex)==0;
  
      // now we can reset the flags
      if ( p.getFaceDataExchangeCounter(faceIndex)==0 ) {
        assertion(p.getNeighbourMergePerformed(faceIndex)==false);
        p.setFaceDataExchangeCounter(faceIndex,TWO_POWER_D); // TODO maybe do not do that here but in the cell? Can be used to determine which cell belongs to skeleton
        p.setNeighbourMergePerformed(faceIndex,true);
      }
    }
  
    // FV
    for (auto& p : exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(destCellDescriptionsIndex)) { // lookup is slow
      mergeWithNeighbourData &= p.getFaceDataExchangeCounter(faceIndex)==0;
  
      // now we can reset the flags
      if ( p.getFaceDataExchangeCounter(faceIndex)==0 ) {
        assertion(p.getNeighbourMergePerformed(faceIndex)==false);
        p.setFaceDataExchangeCounter(faceIndex,TWO_POWER_D);
        p.setNeighbourMergePerformed(faceIndex,true);
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
    const tarch::la::Vector<DIMENSIONS, double>& h,
    const int level) const {
  if ( hasToCommunicate(h) ) {
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
        const int receivedMetadataIndex,
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

      if ( receivedMetadataIndex != InvalidMetadataIndex ) {
          exahype::MetadataHeap::HeapEntries& receivedMetadata =
              MetadataHeap::getInstance().getData(receivedMetadataIndex);
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
    const tarch::la::Vector<DIMENSIONS, double>& h,
    int level) const {
  if ( hasToCommunicate(h) ) {
    dfor2(myDest)
      dfor2(mySrc)
        tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest; // "invert" points
        tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

        int destScalar = TWO_POWER_D - myDestScalar - 1; // "invert" point indices
        int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;

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

          int receivedMetadataIndex = InvalidMetadataIndex;
          if ( receiveNoMetadata==false ) {
            receivedMetadataIndex = exahype::receiveNeighbourCommunicationMetadata(fromRank, x, level);
          }

          if( hasToMergeWithNeighbourData(src,dest) ) {
            if ( mergeWithReceivedData ) {
              mergeWithNeighbourData(
                  fromRank,
                  receivedMetadataIndex,
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
