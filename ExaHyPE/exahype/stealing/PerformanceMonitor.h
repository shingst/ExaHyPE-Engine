#if !defined(_EXAHYPE_STEALING_PERFORMANCEMONITOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_PERFORMANCEMONITOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace stealing {
    class PerformanceMonitor;
    enum class PerformanceMetric;
  }
}

enum class exahype::stealing::PerformanceMetric {
  CurrentTasks,
  RemainingTasks,
  TasksPerTimestep,
  CurrentWaitingTime
};

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
    int *_currentTasksSnapshot;
    /*
     *  A double buffering scheme is used: The buffer stores intermediate
     *  "in-flight" values (for non-blocking MPI).
     */
    int *_currentTasksReceiveBuffer;
    /*
     *  local load counter of a rank, represents current load (i.e. number of tasks
     *  in queue)
     */
    std::atomic<int> _currentTasks;
    /*
     * remaining load uses the additional information of how many tasks will be
     * spawned in a time step, i.e. it tells us how many tasks will still need
     * to be processed before a time step can be completed
     */
    std::atomic<int> _remainingTasks;

    // stores total number of tasks spawned in every time step
    int _tasksPerTimestep;
    // send buffer for repeated MPIIallgather
    int _currentTasksSendBuffer;

    // the current gather requests
    MPI_Request _gatherTasksRequest;
    MPI_Request _gatherWaitingTimesRequest;

    tarch::multicore::BooleanSemaphore _semaphore;

    // make progress on current gather requests
    void progressGather();
    // posts a new gather request
    void postGatherTasks();

    void postGatherWaitingTimes();

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
    // getter for current load snapshot
    const int *getCurrentLoadSnapshot();

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
