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
    std::atomic<bool> _isVictim;
    std::atomic<bool> _emergencyTriggered;

    int               _zeroThreshold;

  public:
    tarch::timing::Watch _iterationTimer;

    static DiffusiveDistributor& getInstance();
    virtual ~DiffusiveDistributor();

    void updateLoadDistribution(int load);

    void updateZeroThreshold(int threshold);

    // return next victim rank
    bool selectVictimRank(int& victim);

    void triggerVictimFlag();
    void resetVictimFlag();

    void triggerEmergency();
    void resetEmergency();

};

#endif
