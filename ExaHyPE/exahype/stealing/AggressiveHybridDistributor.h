#if !defined(_EXAHYPE_STEALING_AggressiveHybridDistributor_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_AggressiveHybridDistributor_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>
#include <vector>

namespace exahype {
  namespace stealing {
    class AggressiveHybridDistributor;
  }
}

/*
 * TODO
 */
class exahype::stealing::AggressiveHybridDistributor {
  private:
    static tarch::logging::Log _log;
    AggressiveHybridDistributor();

    double _temperatureDiffusion;
    double _temperatureCCP;

    double _thresholdTempAdaptation;

    bool _adaptTemperature;

    int _CCPFrequency, _CCPStepsPerPhase;

    int *_initialLoadPerRank;
    int *_newLoadDistribution;
    int *_idealTasksToOffload;
    int *_optimalTasksPerRank;
    // stores how many tasks should be offloaded in every time step
    int *_tasksToOffload;
    // stores how many tasks still need to be offloaded in the current time step
    std::atomic<int> *_remainingTasksToOffload;
    int *_emergenciesPerRank;
    std::atomic<int> *_notOffloaded;
    std::atomic<int> *_actuallyOffloaded;

    int _totalTasksOffloaded;
    int _oldTotalTasksOffloaded;

    int _incrementCurrent;
    int _incrementPrevious;

    bool _isEnabled;

    void updateLoadDistributionCCP();
    void updateLoadDistributionDiffusive();

    int determineCriticalRank();
    void determineOptimalVictim(int& optimalVictim, double& waitingTimeOptimalVictim);

  public:
    static AggressiveHybridDistributor& getInstance();
    virtual ~AggressiveHybridDistributor();

    void configure(double startTempCCP, double startTempDiffusion,
                   int CCPFrequency, int CCPStepsPerPhase,
                   bool adaptTemperature,
                   double thresholdTempAdaptation);

    /*
     *  This operation computes a new unique load distribution (and accordingly, distribution rules)
     *  that aims to achieve the best-case load distribution (i.e. load of a rank is equal to the
     *  average load of all ranks). In order to account for the fact that on some nodes, there may
     *  be less cores/consumers available than on others, the load distribution takes into account
     *  a weight (consumersPerRank) for every rank.
     *  This operation is a global (!) operation that must be called within runGlobalStep on every rank.
     */
    void computeIdealLoadDistribution(int enclaveCells, int skeletonCells);

    void resetRemainingTasksToOffload();
    void printOffloadingStatistics();
    // return next victim rank
    bool selectVictimRank(int& victim, bool& last);

    void getAllVictimRanks(std::vector<int>& victimRanks);
 
    void updateLoadDistribution();
    void handleEmergencyOnRank(int rank);

    void enable();
    void disable();

};

#endif
