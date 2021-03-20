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

#if !defined(_EXAHYPE_OFFLOADING_STATICDISTRIBUTOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_OFFLOADING_STATICDISTRIBUTOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <vector>
#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace reactive {
    class StaticDistributor;
  }
}

/*
 * The StaticDistributor represents a load balancing strategy
 * where the number of tasks sent away to a given rank is
 * fixed for every time step. This number of tasks is either
 * hardcoded in an input file or it is computed after the mesh
 * refinement has finished based on the number of cells that
 * a rank is responsible for.
 */
class exahype::reactive::StaticDistributor {
  private:
    static tarch::logging::Log _log;
    StaticDistributor();

    // stores how many tasks should be offloaded in every time step
    int *_tasksToOffload;
    // stores how many tasks still need to be offloaded in the current time step
    std::atomic<int> *_remainingTasksToOffload;
    // stores the number of consumers per rank that defines a weight for the load balancing
    int *_consumersPerRank;

  public:
    static StaticDistributor& getInstance();
    virtual ~StaticDistributor();

    /*
     *  This operation computes a new unique load distribution (and accordingly, distribution rules)
     *  that aims to achieve the best-case load distribution (i.e. load of a rank is equal to the
     *  average load of all ranks). In order to account for the fact that on some nodes, there may
     *  be less cores/consumers available than on others, the load distribution takes into account
     *  a weight (consumersPerRank) for every rank.
     *  This operation is a global (!) operation that must be called within runGlobalStep on every rank.
     */
    void computeNewLoadDistribution(int enclaveCells, int skeletonCells);

    void resetRemainingTasksToOffload();
    void getAllVictimRanks(std::vector<int>& victims);
    // return next victim rank
    bool selectVictimRank(int& victim);

    void loadDistributionFromFile(const std::string& filename);
};

#endif
