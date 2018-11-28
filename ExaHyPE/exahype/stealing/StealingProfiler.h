#if !defined(_EXAHYPE_STEALING_STEALINGPROFILER_H_)  
#define _EXAHYPE_STEALING_STEALINGPROFILER_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/timing/Watch.h"

#include <atomic>
#include <vector>

#include "exahype/solvers/Solver.h"

namespace exahype {
  namespace stealing {
    class StealingProfiler;
  }
}

/*
 * This class can be used to measure stealing related performance
 * characteristics during an execution run of ExaHype.
 * These statistics are gathered at runtime and printed at the
 * program exit. The StealingProfiler is enabled with the compile-time
 * flag -DStealingUseProfiler.
 */
class exahype::stealing::StealingProfiler {

  private:
    static tarch::logging::Log _log;
    StealingProfiler();
    virtual ~StealingProfiler();

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
    std::atomic<int> _stealingDecisionsPhase;
    std::atomic<int> _stealingDecisions;

    //times per phase
    std::atomic<unsigned long long> _accWaitTasksPhaseTime;
    std::atomic<unsigned long long> _accUsefulCommunicationPhaseTime;
    std::atomic<unsigned long long> _accIdleCommunicationPhaseTime;
    std::atomic<unsigned long long> _accComputationPhaseTime;
    std::atomic<unsigned long long> _accHandlingPhaseTime;
    std::atomic<unsigned long long> _accOffloadPhaseTime;

    //accumulated total (i.e. over all phases) times
    std::atomic<unsigned long long> _accWaitTasksTime;
    std::atomic<unsigned long long> _accUsefulCommunicationTime;
    std::atomic<unsigned long long> _accIdleCommunicationTime;
    std::atomic<unsigned long long> _accComputationTime;
    std::atomic<unsigned long long> _accHandlingTime;
    std::atomic<unsigned long long> _accOffloadTime;

    public:
      static StealingProfiler& getInstance();
      void beginPhase();
      void endPhase();

      void notifyOffloadedTask(int rank);
      void notifyTargetOffloadedTask(int ntasks, int rank);
      void notifyReceivedTask(int rank);
      void notifySpawnedTask();
      void notifyPerformanceUpdate();
      void notifyLatePerformanceUpdate();
      void notifyStealingDecision();
      void notifyThresholdFail();
      void beginComputation();
	  void endComputation(double elapsed);

	  void beginCommunication();
	  void endCommunication(bool successful, double elapsed);

	  void beginHandling();
	  void endHandling(double elapsed);

	  void beginWaitForTasks();
	  void endWaitForTasks(double elapsed);

	  void beginOffload();
	  void endOffload(double elapsed);

	  void printStatistics();

};

#endif
