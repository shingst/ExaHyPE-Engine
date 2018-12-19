#if !defined(_EXAHYPE_STEALING_DIFFUSIVEDISTRIBUTOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_DIFFUSIVEDISTRIBUTOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/timing/Watch.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace stealing {
    class DiffusiveDistributor;
  }
}

class exahype::stealing::DiffusiveDistributor {
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
