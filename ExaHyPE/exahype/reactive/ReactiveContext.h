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
 * It further manages MPI communicators for offloading- or tasksharing-related communication and the blacklist
 * for potentially overloaded ranks.
 */
class exahype::reactive::ReactiveContext {
  public:
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

    enum class ResilienceStrategy {
      None,
      TaskSharing,
      TaskSharingResilienceChecks,
      TaskSharingResilienceCorrection
    };

    static void setOffloadingStrategy(OffloadingStrategy strategy);
    static OffloadingStrategy getOffloadingStrategy();

    static void setResilienceStrategy(ResilienceStrategy strategy);
    static ResilienceStrategy getResilienceStrategy();

    static void setSaveRedundantComputations(bool saveRedundantComputations);
    static bool getSaveRedundantComputations();

    static void setMakeSkeletonsShareable(bool makeSkeletonsShareable);
    static bool getMakeSkeletonsShareable();

  private:
    /**
     * The logging device.
     */
    static tarch::logging::Log _log;

    int _threadId;

    int _maxTag;

    /**
     * Flag is set, if this rank has become a victim rank in the current
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

    ReactiveContext(int threadId);

    /**
     * Communicator for communication between replicating ranks.
     */
    //MPI_Comm _interTeamComm, _interTeamCommKey, _interTeamCommAck;
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

    bool isEnabled();
    
    bool usesOffloading();

    //Todo: with these two function, we can clean up the interface and make some more functions private
    void initializeCommunicatorsAndTeamMetadata();

    void destroy();

    int getOffloadingTag();

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

    bool isBlacklisted(int rank);
    bool isEmergencyTriggered();
    bool isEmergencyTriggeredOnRank(int rank);
    void triggerEmergencyForRank(int rank);

    void recoverBlacklistedRanks();
    void printBlacklist();

    static ReactiveContext& getInstance();
    virtual ~ReactiveContext();
};

#endif
