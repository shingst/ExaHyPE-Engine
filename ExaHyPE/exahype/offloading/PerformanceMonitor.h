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

    bool _terminatedGlobally;

    double *_currentWaitingTimesSnapshot; //length nranks*nranks
    //int *_currentWaitingTimesReceiveBuffer; //length nrank*nranks
    //int *_currentWaitingTimesSendBuffer; //length nranks
    double *_currentWaitingTimes; //lenght nranks

    double *_currentBlacklistSnapshot;
    //double *_currentBlacklistReceiveBuffer;
    //double *_currentBlacklistSendBuffer;
    double *_currentBlacklist;

    double *_currentFusedDataReceiveBuffer;
    double *_currentFusedDataSendBuffer;

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
    MPI_Request _allreduceBlacklistRequest;
    MPI_Request _fusedGatherRequest;

    tarch::multicore::BooleanSemaphore _semaphore;

    // make progress on current gather requests
    void progressGather();
    // posts a new gather request
    void postGatherTasks();

    void postGatherWaitingTimes();

    void postAllreduceBlacklist();

    void postFusedRequest();

  public:
    // signals that a rank has finished computing any local work
    void stop();

    void submitWaitingTimeForRank(double waitingTime, int rank);
    const double *getWaitingTimesSnapshot();

    void submitBlacklistValueForRank(double bval, int rank);
    const double *getBlacklistSnapshot();

    void setCurrentTasks(int numTasks);
    // increases the current load, to be called when a new task
    // is created
    void incCurrentTasks();
    // decreases the current load
    void decCurrentTasks();

    // sets the local load per time step (needs to be called again after
    // mesh refinement)
    void setTasksPerTimestep(int numTasks);
    // getter for remaining load for the current time step
    int getRemainingTasks();
    // getter for local load per time step
    int getTasksPerTimestep();
    // getter for current load snapshot
    const int *getCurrentTasksSnapshot();

    // decreases remaining load for current time step
    void decRemainingTasks();

    // returns true, if every rank has finished computing during
    // an ExahyPE run
    bool isGloballyTerminated();
    void run();

    static PerformanceMonitor& getInstance();
    virtual ~PerformanceMonitor();
};

#endif
