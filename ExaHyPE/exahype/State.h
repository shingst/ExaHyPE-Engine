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
 
#ifndef _EXAHYPE_STATE_H_
#define _EXAHYPE_STATE_H_

#include "exahype/records/State.h"
#include "peano/grid/State.h"

#include <vector>
#include <memory>

#include "peano/grid/Checkpoint.h"

namespace exahype {
class State;

/**
 * Forward declaration
 */
class Vertex;
/**
 * Forward declaration
 */
class Cell;

namespace repositories {
  /**
   * Forward declaration
   */
  class RepositoryArrayStack;
    class RepositorySTDStack;
  }
}



/**
 * Blueprint for solver state.
 *
 * This file has originally been created by the PDT and may be manually extended
 *to
 * the needs of your application. We do not recommend to remove anything!
 */
class exahype::State : public peano::grid::State<exahype::records::State> {
 private:
  typedef class peano::grid::State<exahype::records::State> Base;

  /**
   * Needed for checkpointing.
   */
  friend class exahype::repositories::RepositoryArrayStack;
  friend class exahype::repositories::RepositorySTDStack;

  void writeToCheckpoint(
      peano::grid::Checkpoint<Vertex, Cell>& checkpoint) const;
  void readFromCheckpoint(
      const peano::grid::Checkpoint<Vertex, Cell>& checkpoint);

 public:
  /**
   * This enum is used to select certain solvers
   * in mappings like PredictionRerun, MeshRefinement etc.
   */
  enum class AlgorithmSection {
    /*
     * The runner is currently
     * performing a normal ADER-DG / FV / ... time step.
     */
    TimeStepping,


    /*
     * Currently performing mesh refinement. Only
     * relevant for merging metadata.
     */
    MeshRefinement,

    /*
     * Currently performing limiter status spreading. Only
     * relevant for merging metadata.
     */
    LimiterStatusSpreading,

    /**
     * In this section, the runner overlaps the
     * operations that must be performed after the
     * mesh refinement with operations that
     * must be performed for a local or global
     * recomputation.
     *
     * This marks the end point of the side branch.
     * Triggers a send for all registered solvers.
     */
    PredictionOrLocalRecomputationAllSend,

    /**
     * In this section, all solver have to drop
     * their messages. Then, the ADER-DG solvers which
     * have violated the CFL condition with their
     * estimated time step size are required
     * to reurn the prediction.
     * Finally, all solvers send out again their
     * face data.
     */
    PredictionRerunAllSend
  };

  /**
   * A flag indicating that the bounding box
   * has been virtually expanded (or not).
   *
   * \note In case the bounding box has been expanded,
   * the computational domain usually shrinks
   * since only those cells are considered
   * as inside which have only inside or
   * boundary vertices.
   * A cell is considered as outside
   * if it has at least one(!) outside vertex.
   *
   * TODO(Dominic): Make private and hide in init function
   */
  static bool VirtuallyExpandBoundingBox;

  /**
   * \return true if we run no batch or if
   * we are in the last iteration of a batch.
   *
   * \note It makes only sense to query the batch state from
   * within a mapping.
   */
  bool isFirstIterationOfBatchOrNoBatch() const;

  /**
   * \return true if we run no batch or if
   * we are in the first iteration of a batch.
   *
   * \note It makes only sense to query the batch state from
   * within a mapping.
   */
  bool isLastIterationOfBatchOrNoBatch() const;

  /**
   * Default Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it.
   */
  State();

  /**
   * Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it. It is kind of a copy constructor that converts an object which
   * comprises solely persistent attributes into a full attribute. This very
   * functionality is implemented within the super type, i.e. this constructor
   * has to invoke the correponsing super type's constructor and not the super
   * type standard constructor.
   */
  State(const Base::PersistentState& argument);

  /**
   * Merge this state with another state
   *
   * @todo Clarify which stuff has to be merged
   */
  void merge(const State& anotherState);

  /**
   * Set to true if we need to exchange local solver
   * data between master and worker at least for one cell at
   * a master-worker boundary.
   *
   * These local solver data are usually restricted or prolongated degrees of freedom.
   * They must not be confused with global solver data such as, e.g.
   * admissible time step sizes.
   *
   * \see exahype::mappings::MeshRefinement
   */
  void setVerticalExchangeOfSolverDataRequired(bool state);
  /**
   * \see setVerticalExchangeOfSolverDataRequired
   */
  bool getVerticalExchangeOfSolverDataRequired() const;

  /**
   * Has to be called after the iteration!
   *
   * Please consult Peano guidebook Section 6.3.2 for details.
   */
  void endedGridConstructionIteration(int finestGridLevelPossible);

  /**
   * Please consult Peano guidebook Section 6.3.2 for details.
   */
  enum RefinementAnswer {
    DontRefineYet,
    Refine,
    EnforceRefinement
  };
  RefinementAnswer mayRefine(bool isCreationalEvent, int level) const;

  /**
   * Please consult Peano guidebook Section 6.3.2 for details.
   */
  bool continueToConstructGrid() const;
};

#endif
