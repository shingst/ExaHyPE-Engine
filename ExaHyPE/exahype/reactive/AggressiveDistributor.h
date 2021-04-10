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

#if !defined(_EXAHYPE_STEALING_AggressiveDistributor_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_AggressiveDistributor_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include "exahype/reactive/OffloadingContext.h"

#include <atomic>

namespace exahype {
  namespace reactive {
    class AggressiveDistributor;
  }
}

/*
 * The AggressiveDistributor represents a load balancing strategy
 * where the number of tasks sent away to a given rank is
 * fixed for every time step. This number of tasks is either
 * hardcoded in an input file or it is computed after the mesh
 * refinement has finished based on the number of cells that
 * a rank is responsible for.
 */
class exahype::reactive::AggressiveDistributor {
  private:

    static tarch::logging::Log _log;
    AggressiveDistributor();

    int _zeroThreshold;

    int *_initialLoadPerRank;
    int *_newLoadDistribution;
    int *_idealTasksToOffload;
    // stores how many tasks should be offloaded in every time step
    int *_tasksToOffload;
    // stores how many tasks still need to be offloaded in the current time step
    std::atomic<int> *_remainingTasksToOffload;
    // stores the number of consumers per rank that defines a weight for the load balancing
    int *_consumersPerRank;

    int *_notOffloaded;
    bool _isEnabled;

  public:
    static AggressiveDistributor& getInstance();
    virtual ~AggressiveDistributor();

    /*
     *  This operation computes a new unique load distribution (and accordingly, distribution rules)
     *  that aims to achieve the best-case load distribution (i.e. load of a rank is equal to the
     *  average load of all ranks). In order to account for the fact that on some nodes, there may
     *  be less cores/consumers available than on others, the load distribution takes into account
     *  a weight (consumersPerRank) for every rank.
     *  This operation is a global (!) operation that must be called within runGlobalStep on every rank.
     */
    void computeIdealLoadDistribution(int enclaveCells, int skeletonCells);

    void resetRemainingTasksToOffload();
    void printOffloadingStatistics();
    // return next victim rank
    bool selectVictimRank(int& victim);
 
    void updateLoadDistribution();
    void handleEmergencyOnRank(int rank);

    void enable();
    void disable();

};

#endif
