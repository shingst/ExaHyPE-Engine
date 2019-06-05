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

#if !defined(_EXAHYPE_STEALING_DIFFUSIVEDISTRIBUTOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_DIFFUSIVEDISTRIBUTOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/timing/Watch.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace offloading {
    class DiffusiveDistributor;
  }
}

class exahype::offloading::DiffusiveDistributor {
  private:
    static tarch::logging::Log _log;
    DiffusiveDistributor();

    // stores how many tasks should be offloaded in every time step
    int *_tasksToOffload;
    // stores how many tasks still need to be offloaded in the current time step
    std::atomic<int> *_remainingTasksToOffload;
 
    int _zeroThreshold;

  public:
    static DiffusiveDistributor& getInstance();
    virtual ~DiffusiveDistributor();

    void updateLoadDistribution();

    void updateZeroThreshold(int threshold);

    void handleEmergencyOnRank(int rank);

    // return next victim rank
    bool selectVictimRank(int& victim);

};

#endif
