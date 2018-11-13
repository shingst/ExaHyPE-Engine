#if !defined(_EXAHYPE_STEALING_PERFORMANCEMONITOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_PERFORMANCEMONITOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace stealing {
    class PerformanceMonitor;
  }
}

/*
 * The PerformanceMonitor stores and distributes on-line live performance
 * information which can be used to make effective offloading decisions.
 */
class exahype::stealing::PerformanceMonitor {
  private:
    static tarch::logging::Log _log;

    PerformanceMonitor();
    // status flag, if false then a rank has terminated locally
    bool _isStarted;
    // here, current global view on the load information is stored
    int *_currentLoadSnapshot;
    /*
     *  A double buffering scheme is used: The buffer stores intermediate
     *  "in-flight" values (for asynchronous MPI).
     */
    int *_currentLoadBuffer;
    /*
     *  local load counter of a rank, represents current load (i.e. number of tasks
     *  in queue)
     */
    std::atomic<int> _currentLoadLocal;
    /*
     * remaining load uses the additional information of how many tasks will be
     * spawned in a time step, i.e. it tells us how many tasks will still need
     * to be processed before a time step can be completed
     */
    std::atomic<int> _remainingLoadLocal;

    // stores total number of tasks spawned in every time step
    int _localLoadPerTimestep;
    // send buffer for repeated MPIIallgather
    int _currentLoadLocalBuffer;

    // the current gather request
    MPI_Request _gather_request;

    tarch::multicore::BooleanSemaphore _semaphore;

    // make progress on current gather request
    void progressGather();
    // posts a new gather request
    void postGather();

  public:
    // signals that a rank has finished computing any local work
    void stop();

    void setCurrentLoad(int load);
    // increases the current load, to be called when a new task
    // is created
    void incCurrentLoad();
    // decreases the current load
    void decCurrentLoad();

    // sets the local load per time step (needs to be called again after
    // mesh refinement)
    void setLocalLoadPerTimestep(int load);
    // getter for remaining load for the current time step
    int getRemainingLocalLoad();
    // getter for local load per time step
    int getLocalLoadPerTimestep();

    // decreases remaining load for current time step
    void decRemainingLocalLoad();

    // returns true, if every rank has finished computing during
    // an ExahyPE run
    bool isGloballyTerminated();
    void run();

    static PerformanceMonitor& getInstance();
    virtual ~PerformanceMonitor();
};

#endif
