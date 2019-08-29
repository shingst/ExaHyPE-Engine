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

#if !defined(_EXAHYPE_STEALING_PERFORMANCEMONITOR_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_PERFORMANCEMONITOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace offloading {
    class PerformanceMonitor;
  }
}

/**
 * The PerformanceMonitor stores and distributes on-line live performance
 * information which can be used to make effective offloading decisions.
 */
class exahype::offloading::PerformanceMonitor {
  private:
    static tarch::logging::Log _log;

    PerformanceMonitor();

    /**
     * Status flag, if false then a rank has stopped
     * the PerformanceMonitor locally.
     */
    bool _isRankActive;

    bool _isDisabled;

    /**
     * Indicates if all ranks have terminated.
     */
    bool _terminatedGlobally;

    double *_currentWaitingTimesSnapshot; //length nranks*nranks
    double *_currentWaitingTimes; //lenght nranks

    double *_currentBlacklistSnapshot;
    double *_currentBlacklist;

    int *_currentTasksSnapshot;

    double *_currentFusedDataReceiveBuffer;
    double *_currentFusedDataSendBuffer;

    /**
     *  local load counter of a rank, represents current load (i.e. number of tasks
     *  in queue)
     */
    std::atomic<int> _currentTasksLocal;

    /**
     * remaining load uses the additional information of how many tasks will be
     * spawned in a time step, i.e. it tells us how many tasks will still need
     * to be processed before a time step can be completed
     */
    std::atomic<int> _remainingTasks;

    /**
     * stores total number of tasks spawned in every time step
     */
    int _tasksPerTimestep;

    /**
     * The current gather request
     */
    MPI_Request _fusedGatherRequest;

    tarch::multicore::BooleanSemaphore _semaphore;

    // make progress on current gather requests
    void progressGather();

    void postFusedRequest();

  public:
    /**
     * Signals that a rank has finished computing any local work
     */
    void stop();

    /**
     * Submits waiting time for a rank.
     * @param waitingTime
     * @param rank
     */
    void submitWaitingTimeForRank(double waitingTime, int rank);

    /**
     * Returns pointer to the current waiting times snapshot.
     * @return
     */
    const double *getWaitingTimesSnapshot();

    /**
     * Submits a new blacklist value for a given rank.
     * @param bval
     * @param rank
     */
    void submitBlacklistValueForRank(double bval, int rank);

    /**
     * Returns current blacklist snapshot.
     * @return Current global blacklist snapshot.
     */
    const double *getBlacklistSnapshot();

    /**
     * Sets the current number of tasks to a given
     * value.
     * @param numTasks
     */
    void setCurrentTasks(int numTasks);

    /**
     * increases the current load, to be called when a new task
     */
    void incCurrentTasks();

    /**
     * decreases the current load
     */
    void decCurrentTasks();

    /**
     * sets the local load per time step (needs to be called again after
     * mesh refinement
     */
    void setTasksPerTimestep(int numTasks);

    /**
     * getter for remaining load for the current time step
     */
    int getRemainingTasks();

    /**
     * getter for local load per time step
     */
    int getTasksPerTimestep();

    /**getter for current number of tasks snapshot
     *
     * @return Snapshot containing the current number of
     * tasks on every MPI rank.
     */
    const int *getCurrentTasksSnapshot();

    /**
     * decrease remaining load for current time step
     */
    void decRemainingTasks();

    /**
     * returns true, if every rank has finished computing during
     * an ExahyPE run
     */
    bool isGloballyTerminated();

    /**
     * Entry routine for the performance monitor.
     * This routine makes progress and ensures
     * that performance information is distributed
     * through the system.
     */
    void run();

    void disable();

    static PerformanceMonitor& getInstance();
    virtual ~PerformanceMonitor();
};

#endif
