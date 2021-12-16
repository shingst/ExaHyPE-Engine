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

#if !defined(_EXAHYPE_OFFLOADING_PERFORMANCEMONITOR_H_)  && defined(Parallel)
#define _EXAHYPE_OFFLOADING_PERFORMANCEMONITOR_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>
#include <vector>

namespace exahype {
  namespace reactive {
    class PerformanceMonitor;
  }
}

/**
 * The PerformanceMonitor stores and distributes on-line live performance
 * information which can be used to make effective offloading decisions.
 * It receives the raw performance data (waiting times, blacklist values) and makes
 * it available to other ranks. Global performance information can be retrieved using
 * getter methods.
 * It is important to regularly call the progress method in the application to ensure
 * that non-blocking MPI messages for keeping the global performance information view up-to-date actually
 * progress in the background.
 */
class exahype::reactive::PerformanceMonitor {
  private:
    static tarch::logging::Log _log;

    PerformanceMonitor();

    /*PerformanceMonitor(const PerformanceMonitor& other) = delete;
    PerformanceMonitor& operator=(const PerformanceMonitor& other) = delete;

    PerformanceMonitor(PerformanceMonitor&& other) = delete;
    PerformanceMonitor& operator=(PerformanceMonitor&& other) = delete;*/

    /**
     * Makes progress on outstanding Iallgather, receives a new global performance snapshot (if possible) and posts a new gather operation if necessary.
     */
    void progressGatherAndCollectNewSnapshot();

    /**
     * Copies received waiting times to snapshot, reduces rank-local blacklists and sets global termination status.
     */
    void processCompletedFusedRequestAndSetTerminationStatus();

    /**
     * Posts a new iallgather request, which requires to copy the rank-local performance data to a send buffer.
     */
    void copyToSendBufferAndPostFusedRequest();

    /**
     * Status flag, if false then a rank has stopped
     * the PerformanceMonitor locally.
     */
    bool _hasTerminatedLocally;

    /**
     * Status flag set to true if performance monitoring is disabled.
     */
    bool _isDisabled;

    /**
     * Indicates if all ranks have terminated.
     */
    bool _terminatedGlobally;

    /**
     * Array of current global view of waiting times (matrix with length nranks*nranks).
     */
    std::vector<double> _currentWaitingTimesGlobalSnapshot;

    /**
     * Stores the rank-local current waiting times (length nranks).
     */
    std::vector<double> _currentWaitingTimesLocal;

    /**
     * Stores global snapshot of blacklist (length nranks*nranks).
     */
    std::vector<double> _currentBlacklistGlobalSnapshot;

    /**
     * Stores the rank-local blacklist (length nranks)
     */
    std::vector<double> _currentBlacklistLocal;

    /**
     * Stores snapshot of current number of tasks per rank (length nranks).
     */
    std::vector<int> _currentTasksSnapshot;

    /**
     * Buffer for the currently in-flight iallgather messages.
     * Contains:
     * 1) Waiting times of local rank for other ranks (length nranks)
     * 2) Local blacklist values for each other rank (length nranks)
     * 3) Currently outstanding tasks on this rank (length 1)
     * 4) Local termination signal: -1 if local rank has terminated, 0 otherwise
     */
    std::vector<double> _currentFusedDataSendBuffer;

    /**
     * The receive buffer stores the global allgathered view of the rank local performance infos.
     * As each rank sends data of size nranks*2+2, it has length nranks*(nranks*2+2)
     */
    std::vector<double> _currentFusedDataReceiveBuffer;

    /**
     *  Local load counter of a rank, represents current load (i.e. number of tasks
     *  in queue)
     */
    std::atomic<int> _currentTasksLocal;

    /**
     * Stores total number of tasks spawned in every time step.
     */
    int _tasksPerTimestep;

    /**
     * The current gather request.
     */
    MPI_Request _fusedGatherRequest;

    /**
     * Semaphore for thread-safety.
     */
    tarch::multicore::BooleanSemaphore _semaphore;

    public:
    /**
     * Signals that a rank has finished computing any local work
     */
    void signalLocalTermination();

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
    const std::vector<double>& getWaitingTimesGlobalSnapshot() const;

    /**
     * Submits a new blacklist value for a given rank.
     * @param bval Blacklist value
     * @param rank Rank which is blacklisted
     */
    void submitBlacklistValueForRank(double bval, int rank);

    /**
     * Returns current blacklist snapshot.
     * @return Current global blacklist snapshot.
     */
    const std::vector<double>& getBlacklistGlobalSnapshot() const;

    /**
     * Sets the current number of tasks to a given
     * value.
     * @param numTasks
     */
    void setCurrentNumTasks(int numTasks);

    /**
     * Increases the current load, to be called when a new task is added.
     */
    void incCurrentNumTasks();

    /**
     * Decreases the current load.
     */
    void decCurrentNumTasks();

    /**
     * Sets the local load per time step (needs to be called again after
     * mesh refinement.
     */
    void setTasksPerTimestep(int numTasks);

    /**
     * Getter for local load per time step
     */
    int getTasksPerTimestep() const;

    /**
     * Getter for current number of tasks snapshot
     *
     * @return Snapshot containing the current number of
     * tasks on every MPI rank.
     */
    const std::vector<int>& getCurrentTasksGlobalSnapshot() const;

    /**
     * @return True, if every rank has finished computing during
     * an ExahyPE run
     */
    bool isGloballyTerminated();

    /**
     * Entry routine for the performance monitor.
     * This routine makes progress and ensures
     * that performance information is distributed
     * through the system.
     */
    void progress();

    /**
     * Disables performance monitoring.
     */
    void disable();

    static PerformanceMonitor& getInstance();
};

#endif
