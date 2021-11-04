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

#if !defined(_EXAHYPE_OFFLOADINGMANAGER_H_)  && defined(Parallel)
#define _EXAHYPE_OFFLOADINGMANAGER_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"


#include <mpi.h>
#include <atomic>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

#if defined(UseMPIThreadSplit)
#define MAX_THREADS 48  // for SuperMUC, may be increased for other architectures
#else
#define MAX_THREADS 1
#endif

namespace exahype {
  namespace reactive {
    class ReactiveContext;
  }
}

/**
 * The ReactiveContext singleton stores the configuration of the reactive mechanisms in ExaHyPE.
 * It further manages MPI communicators for offloading- or tasksharing-related communication and the rank-local(!) blacklist
 * for potentially overloaded ranks.
 */
class exahype::reactive::ReactiveContext {
  public:
    /**
     * Enum for offloading strategy.
     */
    enum class OffloadingStrategy {
      None,
      Dynamic,
      Diffusive,
      Aggressive,
      AggressiveCCP,
      AggressiveHybrid,
      Static,
      StaticHardcoded
    };

    /**
     * Enum for resilience strategy.
     */
    enum class ResilienceStrategy {
      None,
      TaskSharing,
      TaskSharingResilienceChecks,
      TaskSharingResilienceCorrection
    };

    //todo: Should probably be done in constructor? Code assumes offloading and resilience strategies to be static.
    /**
     * Sets offloading strategy.
     * @param strategy The offloading strategy to be used.
     */
    static void setOffloadingStrategy(OffloadingStrategy strategy);

    /**
     * @return Offloading strategy.
     */
    static OffloadingStrategy getOffloadingStrategy();

    /**
     * Sets resilience strategy.
     * @param The resilience strategy to be used.
     */
    static void setResilienceStrategy(ResilienceStrategy strategy);

    /**
     * @return Resilience strategy.
     */
    static ResilienceStrategy getResilienceStrategy();

    /**
     * Sets flags that indicates whether redundant computations should be saved
     * through reactive task outcome sharing.
     */
    static void setSaveRedundantComputations(bool saveRedundantComputations);

    /**
     * @return True if reactive task outcome sharing should be used to save redundant computations.
     */
    static bool getSaveRedundantComputations();

    /**
     * Controls whether skeleton tasks should be shared or not.
     * Per default only enclaves are shared.
     */
    static void setMakeSkeletonsShareable(bool makeSkeletonsShareable);

    /**
     * @return True if skeletons should be shared.
     */
    static bool getMakeSkeletonsShareable();

    /**
     * @return True if some reactive feature is enabled (load balancing or resilience).
     */
    static bool isReactivityEnabled();

    /**
     * @return True if reactive offloading is enabled.
     */
    static bool isReactiveOffloadingEnabled();

  private:
    /**
     * The logging device.
     */
    static tarch::logging::Log _log;

    //singleton per thread
    ReactiveContext(int threadId);
    virtual ~ReactiveContext();

    /**
     * Thread id of the thread the reactive context object belongs to.
     */
    int _threadId;

    /**
     * Maximum supported tag by MPI.
     */
    int _maxSupportedTag;

    /**
     * Flag is set if this rank has become a victim rank in the current
     * time step.
     */
    std::atomic<bool> _isVictim;

    /**
     * Flag indicating whether this rank has triggered an emergency.
     */
    std::atomic<bool> _emergencyTriggered;

    /**
     * Local instance of blacklist.
     */
    double *_localBlacklist;

    /**
     * Flag indicating whether the rank has notified
     * all its victims that no more tasks will be
     * offloaded in this time step.
     */
    bool _hasNotifiedSendCompleted;

    /**
     * Communicator for communication between replicating ranks.
     */
    static int _numTeams;
    static int _interTeamRank;

    static ReactiveContext* _static_managers[MAX_THREADS];

    static MPI_Comm  _offloadingComms[MAX_THREADS];
    static MPI_Comm  _offloadingCommsMapped[MAX_THREADS];

    static MPI_Comm  _interTeamComms[MAX_THREADS];

    static OffloadingStrategy ChosenOffloadingStrategy;

    static ResilienceStrategy ChosenResilienceStrategy;

    static bool SaveRedundantComputations;
    static bool MakeSkeletonsShareable;
  public:
    //Todo: with these two function, we can clean up the interface and make some more functions private
    void initializeCommunicatorsAndTeamMetadata();

    void destroy();

    int getNextMPITag();

    void setTMPIInterTeamCommunicators(MPI_Comm comm, MPI_Comm commKey, MPI_Comm commAck);
    MPI_Comm getTMPIInterTeamCommunicatorData();

    void setTMPINumTeams(int team);
    unsigned int getTMPINumTeams();

    void setTMPITeamNumber(int interTeamRank);
    int getTMPITeamNumber();

    /**
     * Creates offloading MPI communicators.
     */
    static void createMPICommunicators();
    /**
     * Destroys offloading MPI communicators.
     */
    static void destroyMPICommunicators();

    /**
     * Returns MPI communicator used for
     * distributing load information and
     * sending/receiving tasks to/from a rank.
     */
    MPI_Comm getMPICommunicator();

    /**
     * Return MPI communicator for sending/receiving back
     * results of an offloaded task.
     */
    MPI_Comm getMPICommunicatorMapped();

    /**
     * Given the current load situation and global knowledge of the performance
     * of the other ranks, this method selects a victim rank for a
     * stealable task, i.e.,
     * a rank to which local work should be offloaded in order to
     * improve load balance.
     */
    bool selectVictimRank(int& victim, bool& last);

    void triggerVictimFlag();
    void resetVictimFlag();
    bool isVictim();

    static ReactiveContext& getInstance();
};

#endif
