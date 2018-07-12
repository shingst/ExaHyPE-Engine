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

#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "exahype/solvers/Solver.h"

#include "tarch/parallel/NodePool.h"

#include <limits>

bool exahype::State::VirtuallyExpandBoundingBox = false;

#ifdef Parallel
bool exahype::State::BroadcastInThisIteration = true;
bool exahype::State::ReduceInThisIteration    = false;
#endif

bool exahype::State::isFirstIterationOfBatchOrNoBatch() const {
  return _stateData.getTotalNumberOfBatchIterations()==1 || _stateData.getBatchIteration()==0;
}

bool exahype::State::isSecondIterationOfBatchOrNoBatch() const {
  return _stateData.getTotalNumberOfBatchIterations()==1 || _stateData.getBatchIteration()==1;
}

bool exahype::State::isLastIterationOfBatchOrNoBatch() const {
  return _stateData.getTotalNumberOfBatchIterations()==1 || _stateData.getBatchIteration()==_stateData.getTotalNumberOfBatchIterations()-1;
}

bool exahype::State::isSecondToLastIterationOfBatchOrNoBatch() const {
  return _stateData.getTotalNumberOfBatchIterations()==1 || _stateData.getBatchIteration()==_stateData.getTotalNumberOfBatchIterations()-2;
}

exahype::State::State() : Base() {
  // @todo Guidebook

  _stateData.setMaxRefinementLevelAllowed(0);
}

exahype::State::State(const Base::PersistentState& argument) : Base(argument) {
  // do nothing
}

void exahype::State::setVerticalExchangeOfSolverDataRequired(bool state) {
  _stateData.setVerticalExchangeOfSolverDataRequired(state);
}

bool exahype::State::getVerticalExchangeOfSolverDataRequired() const {
  return _stateData.getVerticalExchangeOfSolverDataRequired();
}

/**
 * \see exahype/State.def
 */
void exahype::State::setAllSolversAttainedStableStateInPreviousIteration(const bool state) {
  _stateData.setAllSolversAttainedStableStateInPreviousIteration(state);
}
/**
 * \see exahype/State.def
 */
bool exahype::State::getAllSolversAttainedStableStateInPreviousIteration() const {
  return _stateData.getAllSolversAttainedStableStateInPreviousIteration();
}

/**
 * \see exahype/State.def
 */
void exahype::State::setMeshRefinementHasConverged(const bool state) {
  _stateData.setMeshRefinementHasConverged(state);
}
/**
 * \see exahype/State.def
 */
bool exahype::State::getMeshRefinementHasConverged() const {
  return _stateData.getMeshRefinementHasConverged();
}

void exahype::State::mergeWithMaster(const exahype::State& anotherState) {
  _stateData.setVerticalExchangeOfSolverDataRequired(
      _stateData.getVerticalExchangeOfSolverDataRequired() ||
      anotherState._stateData.getVerticalExchangeOfSolverDataRequired());
  _stateData.setAllSolversAttainedStableStateInPreviousIteration(
      _stateData.getAllSolversAttainedStableStateInPreviousIteration() &&
      anotherState._stateData.getAllSolversAttainedStableStateInPreviousIteration());
}

void exahype::State::writeToCheckpoint(
    peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>& checkpoint) const {
  // do nothing
}

void exahype::State::readFromCheckpoint(
    const peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>& checkpoint) {
  // do nothing
}

void exahype::State::endedGridConstructionIteration(int finestGridLevelPossible) {
  const bool idleNodesLeft =
      tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0;
  const bool nodePoolHasGivenOutRankSizeLastQuery =
      tarch::parallel::NodePool::getInstance().hasGivenOutRankSizeLastQuery();

  #ifdef Debug
  std::cout <<  "!getHasChangedVertexOrCellState=" << !_stateData.getHasChangedVertexOrCellState() << std::endl;
  std::cout <<  "!getHasRefined=" << !_stateData.getHasRefined() << std::endl;
  std::cout <<  "!getHasErased=" << !_stateData.getHasErased()  << std::endl;
  std::cout <<  "!getHasTriggeredRefinementForNextIteration=" << !_stateData.getHasTriggeredRefinementForNextIteration() << std::endl;
  std::cout <<  "!getHasTriggeredEraseForNextIteration=" << !_stateData.getHasTriggeredEraseForNextIteration() << std::endl;
  #ifdef Parallel
  std::cout <<  "!getCouldNotEraseDueToDecompositionFlag=" << !_stateData.getCouldNotEraseDueToDecompositionFlag() << std::endl;
  #endif
  std::cout << "isGridStationary()=" << isGridStationary() << std::endl;
  std::cout << "getMaxRefinementLevelAllowed()=" << _stateData.getMaxRefinementLevelAllowed() << std::endl;
  #endif

  // No more nodes left. Start to enforce refinement
  if ( !idleNodesLeft
      && _stateData.getMaxRefinementLevelAllowed()>=0
      && !nodePoolHasGivenOutRankSizeLastQuery) {
    _stateData.setMaxRefinementLevelAllowed(-1);
  }
  // Seems that max permitted level has exceeded max grid level. We may assume
  // that there are more MPI ranks than available trees.
  else if (isGridStationary()
      && _stateData.getMaxRefinementLevelAllowed()>finestGridLevelPossible
      && _stateData.getMaxRefinementLevelAllowed()>=0
      && !nodePoolHasGivenOutRankSizeLastQuery) {
    _stateData.setMaxRefinementLevelAllowed( -1 );
  }
  // Reset counter by two. Some LB has happened and we might wanna
  // give the whole system two sweeps to recover from this LB, i.e. to
  // set up all partitions properly and recompute all LB metrics.
  else if (nodePoolHasGivenOutRankSizeLastQuery
      && _stateData.getMaxRefinementLevelAllowed()>=2) {
    _stateData.setMaxRefinementLevelAllowed(
        _stateData.getMaxRefinementLevelAllowed()-2);
  }
  // Refinement is enforced. So we decrease counter. Once we underrun -2, grid
  // construction can terminate as all enforced refined instructions went
  // through.
  else if (_stateData.getMaxRefinementLevelAllowed()<=-1
      && !nodePoolHasGivenOutRankSizeLastQuery
      && isGridStationary()) {
    _stateData.setMaxRefinementLevelAllowed(
        _stateData.getMaxRefinementLevelAllowed()-1 );
  }
  // Nothing has changed in this grid iteration in the grid and we haven't
  // given out new workers. So increase the permitted maximum grid level by
  // one and give another try whether the grid adds more vertices.
  else if (
      (!nodePoolHasGivenOutRankSizeLastQuery)
      && isGridStationary()
      && (_stateData.getMaxRefinementLevelAllowed()>=0)
  ) {
    _stateData.setMaxRefinementLevelAllowed(
        _stateData.getMaxRefinementLevelAllowed()+1);
  }
}


exahype::State::RefinementAnswer exahype::State::mayRefine(bool isCreationalEvent, int level) const
{
#ifdef Parallel
  if (
      _stateData.getMaxRefinementLevelAllowed()<=-2
      &&
      isCreationalEvent
      &&
      !isInvolvedInJoinOrFork() // A Peano assertion was triggered
  ) {
    return RefinementAnswer::EnforceRefinement;
  }
  else if ( _stateData.getMaxRefinementLevelAllowed()<0 ) {
    return RefinementAnswer::Refine;
  }
  else if (
      _stateData.getMaxRefinementLevelAllowed()>level
      &&
      !isCreationalEvent
      &&
      mayForkDueToLoadBalancing()
  ) {
    return RefinementAnswer::Refine;
  }
  else {
    return RefinementAnswer::DontRefineYet;
  }
#else
  return RefinementAnswer::Refine;
#endif
}


bool exahype::State::continueToConstructGrid() const {
#ifdef Parallel
  return _stateData.getMaxRefinementLevelAllowed()>=-3 && _stateData.getMeshRefinementHasConverged();
#else
  return !isGridBalanced() && _stateData.getMeshRefinementHasConverged();
#endif
}
