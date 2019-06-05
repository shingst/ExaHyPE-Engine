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

#if !defined(_EXAHYPE_STEALING_DYNAMICDISTRIBUTOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_DYNAMICDISTRIBUTOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace offloading {
    class DynamicDistributor;
  }
}

/*
 * The DynamicDistributor represents a load balancing strategy
 * where it is dynamically determined whether a task needs to
 * be given away at runtime based on the PerformanceMonitor's
 * information. Compared to the StaticDistributor, this
 * strategy is less aggressive as it will give away only
 * one task at a time until a new load update has been received.
 */
class exahype::offloading::DynamicDistributor {
  private:
    static tarch::logging::Log _log;
    DynamicDistributor();

    int *_tasksToOffload;
    std::atomic<int> *_remainingTasksToOffload;
    int *_consumersPerRank;

  public:
    static DynamicDistributor& getInstance();
    virtual ~DynamicDistributor();

    void computeNewLoadDistribution(int *currentLoadSnapshot);
    bool selectVictimRank(int& victim);
};

#endif
