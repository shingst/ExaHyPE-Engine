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
#ifndef MAX_THREADS
#define MAX_THREADS 48  // for SuperMUC, may be increased for other architectures
#endif
#else
#ifndef MAX_THREADS
#define MAX_THREADS 1
#endif
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

    /**
     * Indicates that this rank is a victim rank (i.e., it has received an offloaded task).
     */
    static void triggerVictimFlag();

    /**
     * Reset victim flag (to be called after a time step).
     */
    static void resetVictimFlag();

    /**
     * @return If this rank is a victim.
     */
    static bool isVictim();

    static void setTMPINumTeams(int team);
    static unsigned int getTMPINumTeams();

    static void setTMPITeamNumber(int interTeamRank);
    static int getTMPITeamNumber();

    static double getResilienceChecksTimeout(); 
    static void setResilienceChecksTimeout(double timeout);

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
    static int MaxSupportedTag;

    /**
     * Flag is set if this rank has become a victim rank in the current
     * time step.
     */
    static std::atomic<bool> IsVictim;

    /**
     * Flag indicating whether the rank has notified
     * all its victims that no more tasks will be
     * offloaded in this time step.
     */
    bool _hasNotifiedSendCompleted;

    /**
     * Stores the number of teaMPI teams.
     */
    static int NumTeams;

    /**
     * Stores the team number the rank belongs to.
     */
    static int Team;

    /**
     * Array of managers. With UseMPIThreadSplit, there is one reactive context per thread (with each using its own MPI communicators).
     * Per default, there is only a single reactive context (singleton).
     */
    static ReactiveContext* StaticManagers[MAX_THREADS];

    /**
     * Communicators for offloading tasks to a victim rank. There may be multiple if UseMPIThreadSplit is used.
     */
    static MPI_Comm  OffloadingComms[MAX_THREADS];

    /**
     * Communicators for receiving task outcomes from a victim rank. There may be multiple if UseMPIThreadSplit is used.
     */
    static MPI_Comm  OffloadingCommsMapped[MAX_THREADS];

    /**
     * Communicators for communication between replicating ranks.
     */
    static MPI_Comm  InterTeamComms[MAX_THREADS];

    static OffloadingStrategy ChosenOffloadingStrategy;
    static ResilienceStrategy ChosenResilienceStrategy;

    static bool SaveRedundantComputations;
    static bool MakeSkeletonsShareable;
    static double ResilienceChecksTimeout;
  public:
    /**
     * Initializes a reactive context object.
     */
    void initialize();

    /**
     * Destroys a reactive context object.
     */
    void destroy();

    /**
     * Used to get a new tag for communication related to offloading/task sharing.
     */
    int getNextMPITag();


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
    MPI_Comm getMPICommunicator() const;

    /**
     * Return MPI communicator for sending/receiving back
     * results of an offloaded task.
     */
    MPI_Comm getMPICommunicatorMapped() const;

    /**
     * Return MPI communicator for communication between teams (reactive task outcome sharing).
     */
    MPI_Comm getTMPIInterTeamCommunicatorData() const;

    /**
     * Given the current load situation and global knowledge of the performance
     * of the other ranks, this method selects a victim rank for a
     * stealable task, i.e.,
     * a rank to which local work should be offloaded in order to
     * improve load balance.
     */
    bool selectVictimRank(int& victim, bool& last);

    static ReactiveContext& getInstance();
};

#endif
