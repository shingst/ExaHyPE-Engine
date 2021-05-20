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

#if !defined(_EXAHYPE_OFFLOADINGPROFILER_H_)
#define _EXAHYPE_OFFLOADINGPROFILER_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/timing/Watch.h"

#include <atomic>
#include <vector>

#include "exahype/solvers/Solver.h"

namespace exahype {
  namespace reactive {
    class OffloadingProfiler;
  }
}

/**
 * This singleton can be used to measure offloading related performance
 * characteristics during an execution run of ExaHype.
 * These statistics are gathered at runtime and printed at the
 * program exit. The OffloadingProfiler is enabled with the compile-time
 * flag -DOffloadingUseProfiler.
 */
class exahype::reactive::OffloadingProfiler {
  private:
    static tarch::logging::Log _log;
    OffloadingProfiler();
    virtual ~OffloadingProfiler();

    //counters
    std::atomic<int> *_offloadedTasksPerRankPhase;
    std::atomic<int> *_offloadedTasksPerRank;
    std::atomic<int> *_targetOffloadedTasksPerRankPhase;
    std::atomic<int> *_targetOffloadedTasksPerRank;
    std::atomic<int> *_receivedTasksPerRankPhase;
    std::atomic<int> *_receivedTasksPerRank;
    std::atomic<int> _executedTasksPhase;
    std::atomic<int> _executedTasks;
    std::atomic<int> _thresholdFailsPhase;
    std::atomic<int> _thresholdFails;
    std::atomic<int> _spawnedTasksPhase;
    std::atomic<int> _spawnedTasks;
    std::atomic<int> _performanceUpdatesPhase;
    std::atomic<int> _performanceUpdates;
    std::atomic<int> _latePerformanceUpdatesPhase;
    std::atomic<int> _latePerformanceUpdates;
    std::atomic<int> _offloadingDecisionsPhase;
    std::atomic<int> _offloadingDecisions;

    //times per phase
    std::atomic<unsigned long long> _accWaitTasksPhaseTime;
    std::atomic<unsigned long long> _accUsefulCommunicationPhaseTime;
    std::atomic<unsigned long long> _accIdleCommunicationPhaseTime;
    std::atomic<unsigned long long> _accProgressPhaseTime;
    std::atomic<unsigned long long> _accProgressRequestsPhaseTime;
    std::atomic<unsigned long long> _accPollPhaseTime;
    std::atomic<unsigned long long> _accComputationPhaseTime;
    std::atomic<unsigned long long> _accHandlingPhaseTime;
    std::atomic<unsigned long long> _accOffloadPhaseTime;

    //accumulated total (i.e. over all phases) times
    std::atomic<unsigned long long> _accWaitTasksTime;
    std::atomic<unsigned long long> _accUsefulCommunicationTime;
    std::atomic<unsigned long long> _accIdleCommunicationTime;
    std::atomic<unsigned long long> _accProgressTime;
    std::atomic<unsigned long long> _accProgressRequestsTime;
    std::atomic<unsigned long long> _accPollTime;
    std::atomic<unsigned long long> _accComputationTime;
    std::atomic<unsigned long long> _accHandlingTime;
    std::atomic<unsigned long long> _accOffloadTime;

  public:
    static OffloadingProfiler& getInstance();
    void beginPhase();
    void endPhase();

    void notifyOffloadedTask(int rank);
    void notifyTargetOffloadedTask(int ntasks, int rank);
    void notifyReceivedTask(int rank);
    void notifySpawnedTask();
    void notifyPerformanceUpdate();
    void notifyLatePerformanceUpdate();
    void notifyOffloadingDecision();
    void notifyThresholdFail();
    void beginComputation();
    void endComputation(double elapsed);

    void beginCommunication();
    void endCommunication(bool successful, double elapsed);

    void beginProgress();
    void endProgress(double elapsed);

    void beginProgressRequests();
    void endProgressRequests(double elapsed);
    
    void beginPolling();
    void endPolling(double elapsed);
    
    void beginHandling();
    void endHandling(double elapsed);

    void beginWaitForTasks();
    void endWaitForTasks(double elapsed);

    void beginOffload();
    void endOffload(double elapsed);

    void printStatistics();
};

#endif
