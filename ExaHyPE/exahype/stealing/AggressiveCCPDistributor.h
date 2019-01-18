#if !defined(_EXAHYPE_STEALING_AggressiveCCPDistributor_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_AggressiveCCPDistributor_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>

namespace exahype {
  namespace stealing {
    class AggressiveCCPDistributor;
  }
}

/*
 * The AggressiveCCPDistributor represents a load balancing strategy
 * where the number of tasks sent away to a given rank is
 * fixed for every time step. This number of tasks is either
 * hardcoded in an input file or it is computed after the mesh
 * refinement has finished based on the number of cells that
 * a rank is responsible for.
 */
class exahype::stealing::AggressiveCCPDistributor {
  private:
    static tarch::logging::Log _log;
    AggressiveCCPDistributor();

    int _zeroThreshold;

    double _temperature;

    int *_initialLoadPerRank;
    int *_newLoadDistribution;
    int *_idealTasksToOffload;
    // stores how many tasks should be offloaded in every time step
    int *_tasksToOffload;
    // stores how many tasks still need to be offloaded in the current time step
    std::atomic<int> *_remainingTasksToOffload;
    int *_emergenciesPerRank;

    int _totalTasksOffloaded;
    int _oldTotalTasksOffloaded;

    int *_notOffloaded;
    bool _isEnabled;

  public:
    static AggressiveCCPDistributor& getInstance();
    virtual ~AggressiveCCPDistributor();

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
    bool selectVictimRank(int& victim);
 
    void updateLoadDistribution();
    void handleEmergencyOnRank(int rank);

    void enable();
    void disable();

};

#endif
