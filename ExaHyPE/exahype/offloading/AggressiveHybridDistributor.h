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

#if !defined(_EXAHYPE_STEALING_AggressiveHybridDistributor_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_AggressiveHybridDistributor_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include <mpi.h>
#include <atomic>
#include <vector>

namespace exahype {
  namespace offloading {
    class AggressiveHybridDistributor;
  }
}

/**
 * The aggressive hybrid distributor features a hybrid diffusion+
 * chains-on-chains partitioning approach
 * in order to adapt the non-persistent load re-distribution
 * after every time step. Task migration happens aggressively, i.e. in
 * every iteration, the new load re-distribution may be changed by >=1 tasks
 * compared to the previous one.
 *
 * Chains-on-chains partitioning can optionally be used to derive an initial
 * load re-distribution. This load re-distribution is auto-tuned by diffusion
 * taking into account the time that an MPI rank is idle waiting
 * for MPI communication and the average time it takes to
 * compute an STP task. Temperature factors can be applied in order
 * to possibly speed up convergence.
 */
class exahype::offloading::AggressiveHybridDistributor {
  private:
	/**
	 * The logging device.
	 */
    static tarch::logging::Log _log;

    AggressiveHybridDistributor();

    /**
     * Temperature used for diffusion.
     */
    double _temperatureDiffusion;
    /**
     * Temperature used for chains-on-chains partitioning.
     */
    double _temperatureCCP;

    /**
     * Threshold for ratio of current and previous task increment
     * at which temperature is either increased or decreased.
     */
    double _thresholdTempAdaptation;

    /**
     * Flag indicating if temperature should be dynamically updated.
     */
    bool _adaptTemperature;

    /**
     * The frequency of after how many time steps, a CCP step is conducted.
     */
    int _CCPFrequency;

    /**
     * Number of CCP steps conducted if the algorithm decides to
     * conduct CCP steps.
     */
    int _CCPStepsPerPhase;

    /**
     * Initial load in terms of STP jobs per rank.
     */
    int *_initialLoad;

    /**
     * Solution of the chains-on-chains partitioning.
     */
    int *_idealTasksToOffloadCCP;

    /**
     * Stores the optimal number of tasks to offload in
     * this time step. In reality, we typically
     * offload fewer tasks in order to avoid
     * overloading of the network.
     */
    int *_optimalTasks;
    /**
     * Stores how many tasks the algorithm actually attempts to offload
     * in every time step.
     */
    int *_tasksToOffload;
    /**
     * Counts how many tasks still need to be offloaded in the current time step.
     */
    std::atomic<int> *_remainingTasksToOffload;

    /**
     * Tracks the number of emergency events for every rank.
     */
    int *_emergencies;

    /**
     * Counts number of tasks that were not offloaded (per rank, for statistics).
     */
    std::atomic<int> *_tasksNotOffloaded;
    /**
     * Counts the number of tasks that were actually offloaded (per rank, for statistics).
     */
    std::atomic<int> *_tasksActuallyOffloaded;


    /**
     * Counts total number of tasks offloaded from this rank (statistics).
     */
    int _totalTasksOffloaded;
    /**
     * Stores old number of offloaded tasks (statistics).
     */
    int _oldTotalTasksOffloaded;

    /**
     * Stores how many tasks will be offloaded additionally in the current
     * step when compared to the previous one.
     */
    int _incrementCurrent;
    /**
     * Stores task increment of previous time step (see above).
     */
    int _incrementPrevious;

    /**
     * Flag indicating whether this distributor is actually used.
     */
    bool _isEnabled;

    /**
     * Conducts a CCP step.
     */
    void updateLoadDistributionCCP();
    /**
     * Conducts a diffusive load balancing step.
     */
    void updateLoadDistributionDiffusive();

    /**
     * Computes the current critical rank.
     * @return Rank number of the critical rank.
     */
    int determineCriticalRank();

    /**
     * Computes the current optimal victim rank.
     * @param optimalVictim Output parameter containing the rank number of the optimal victim.
     * @param waitingTimeOptimalVictim Output parameter containing the waiting time of the optimal victim.
     */
    void determineOptimalVictim(int& optimalVictim, double& waitingTimeOptimalVictim);

  public:
    static AggressiveHybridDistributor& getInstance();
    virtual ~AggressiveHybridDistributor();

    /**
     * Configures the aggressive hybrid distribution strategy
     * @param startTempCCP Initial temperature for the CCP
     * @param startTempDiffusion Initial temperature for the diffusion
     * @param CCPFrequency Frequency after which another CCP Phase is conducted
     * @param CCPStepsPerPhase Number of CCP steps per phase
     * @param adaptTemperature Update temperature yes/no
     * @param thresholdTempAdaptation Threshold for increment ratio at which temperature is updated
     */
    void configure(double startTempCCP, double startTempDiffusion,
                   int CCPFrequency, int CCPStepsPerPhase,
                   bool adaptTemperature,
                   double thresholdTempAdaptation);

    /**
     *  This operation computes a new unique load distribution (and accordingly, distribution rules)
     *  that aims to achieve the best-case load distribution (i.e. load of a rank is equal to the
     *  average load of all ranks). In order to keep the task re-migration as local as possible (i.e., between
     *  neighbouring ranks according to the natural ordering of MPI ranks), a so-called chains-on-chains partitioning
     *  is computed where we assume uniform cost for each task.
     *  In order to account for the fact that on some nodes, there may
     *  be less cores/consumers available than on others, the load distribution may take into account
     *  a weight (consumersPerRank) for every rank.
     *  This operation is a global (!) operation that must be called within runGlobalStep on every rank.
     */
    void computeIdealLoadDistribution(int enclaveCells, int skeletonCells);

    /**
     * Updates the remaining tasks to offload counters. This typically happens after
     * a new load re-distribution has been computed.
     */
    void resetRemainingTasksToOffload();

    /**
     * Print some useful statistics for offloading.
     */
    void printOffloadingStatistics();

    /**
     * Returns next victim rank according to the current load re-distribution rules.
     * @param victim Contains the rank number of the next victim rank (output)
     * @param last Indicates if this was the last victim that could be selected (output)
     * @return True if a victim rank was chosen.
     * @Note  ToDo: This should be moved into a super class.
     */
    bool selectVictimRank(int& victim, bool& last);

    /**
     * Inserts ranks to which
     * tasks were offloaded
     * into parameter vector victimRanks.
     * @param victimRanks Result
     */
    void getAllVictimRanks(std::vector<int>& victimRanks);
 
    /**
     * Conducts another update of the load re-distribution. This is
     * typically triggered in every time step.
     */
    void updateLoadDistribution();

    /**
     * Callback from offloading manager that is used to
     * track the number of emergencies.
     * @param rank
     * @Note ToDo: This should be inherited from a super class.
     */
    void handleEmergencyOnRank(int rank);

    /**
     * Enables distributor.
     */
    void enable();
    /**
     * Disables distributor.
     */
    void disable();

};

#endif
