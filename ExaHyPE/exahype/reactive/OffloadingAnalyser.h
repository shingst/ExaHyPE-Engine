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

#ifndef _EXAHYPE_STEALING_STEALINGANALYSER_H_
#define _EXAHYPE_STEALING_STEALINGANALYSER_H_

#include "tarch/logging/Log.h"
#include "tarch/timing/Watch.h"
#include "tarch/timing/GlidingAverageMeasurement.h"

#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/Node.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include "peano/performanceanalysis/Analyser.h"

#include <vector>
#include <atomic>

namespace exahype {
  namespace reactive {
    class OffloadingAnalyser;
  }
}

/**
 * The OffloadingAnalyser is responsible for introspecting ExaHyPE's performance by collecting
 * waiting times in Peano's vertical master-worker communication.
 * Gathered waiting times are distributed to other ranks at runtime using the PerformanceMonitor.
 * The waiting times serve as input to reactive task offloading (e.g., to find overloaded critical ranks).
 */
class exahype::reactive::OffloadingAnalyser : public peano::performanceanalysis::Analyser {
  private:
    static tarch::logging::Log     _log;

    /**
     * Flag that indicates whether the analysis is actually used.
     */
    bool _isSwitchedOn;

    /**
     * Watch for waiting times from worker.
     */
    tarch::timing::Watch           _waitForWorkerDataWatch;

    /**
     * Watch for waiting times from master.
     */
    tarch::timing::Watch           _waitForMasterDataWatch;

    /**
     * Watch for waiting times from global master.
     */
    tarch::timing::Watch           _waitForGlobalMasterDataWatch;

    /**
     * Watch for wall time per time step.
     */
    tarch::timing::Watch           _timeStepWatch;

    /**
     * Stores gliding averages of waiting times.
     */
    std::vector<tarch::timing::GlidingAverageMeasurement>    _waitAvgForOtherRank;

    /**
     * Gliding average over times per STP.
     */
    tarch::timing::GlidingAverageMeasurement     _avgTimePerSTP;

    /**
     * Gliding average over times per time step.
     */
    tarch::timing::GlidingAverageMeasurement     _avgTimePerTimeStep;

    /**
     * The zero threshold is the time below which a waiting time is considered to be equal to zero.
     */
    double _currentZeroThreshold;

    /**
     * Counts number of mesh iterations. Each second one is assumed to be the ending of a time step.
     */ 
    long _iterationCounter;

    /**
     * Accumulates waiting times for workers in order to avoid shadowing effects (callbacks from Peano come one after the other) where a long wait for a worker hides the waiting times for subsequent workers.
     */
    double _currentAccumulatedWorkerTime;

    /**
     * Estimated wall time for all pending background jobs (enclave STPs) when the code waits for a worker.
     */
    double _estWalltimeForPendingJobsAtReceiveForWorker;

    /**
     * Array of filtered waiting times where waiting times below zero threshold are set to zero.
     */
    std::vector<double> _currentFilteredWaitingTimesSnapshot;

    /**
     * After the code waits for a worker, some tasks may still be received with reactive offloading. This counter counts those tasks.
     */
    std::atomic<int> _numLateIncomingJobsWhileWaiting;

    /**
     * Filters out zero waiting times and dynamically adapts the zero threshold below which
     * waiting times are considered to be zero.
     */
    void updateZeroTresholdAndFilteredSnapshot();

  public:
    OffloadingAnalyser();
    //virtual ~OffloadingAnalyser() {};

    OffloadingAnalyser(const OffloadingAnalyser& other) = delete;
    OffloadingAnalyser& operator=(const OffloadingAnalyser& other) = delete;

    /**
     * Returns snapshot of filtered waiting times.
     */
    const std::vector<double>& getFilteredWaitingTimesSnapshot() const;

    /**
     * Returns true if there is an entry > 0 in the waiting times snapshot for each rank (each rank needs to be waiting for at least one rank)
     */
    bool isValidSnapshot() const;

    /**
     * Needs to be called whenever a STP job is received.
     * (internally used to keep track of STPs that are arrive late).
     */
    void notifyReceivedSTPJob();    

    /**
     * Prints filtered waiting times.
     */
    void printWaitingTimes() const;

    /**
     * Resets all gliding averages.
     */
    void resetWaitingTimes();

    static OffloadingAnalyser& getInstance();

    double getZeroThreshold() const;

    void setTimePerSTP(double timePerSTP);
    double getTimePerSTP() const;

    void setTimePerTimeStep(double timePerStep);
    double getTimePerTimeStep() const;

    virtual void enable(bool value) override;

    virtual void beginIteration() override;
    virtual void endIteration(double numberOfInnerLeafCells, double numberOfOuterLeafCells, double numberOfInnerCells, double numberOfOuterCells, double numberOfLocalCells, double numberOfLocalVertices) override;

    virtual void beginToSendDataToWorker();
    virtual void endToSendDataToWorker(int worker);
    virtual void beginToSendDataToMaster();
    virtual void endToSendDataToMaster();
    virtual void beginToReceiveDataFromGlobalMaster();
    virtual void endToReceiveDataFromGlobalMaster();

    virtual void beginToReceiveDataFromWorker() override;
    virtual void endToReceiveDataFromWorker(int fromRank) override;
    virtual void beginToReceiveDataFromMaster(int master) override;
    virtual void endToReceiveDataFromMaster(int master) override;

    /**
      * Nop.
      */
    virtual void addWorker(
      int                                 workerRank,
      int                                 level
    ) override;

    /**
      * Nop.
      */
    virtual void removeWorker(
      int                                 workerRank,
      int                                 level
    ) override;

    /**
      * Nop.
      */
    virtual void enterCentralElementOfEnclosingSpacetree() override;

    /**
      * Nop.
      */
    virtual void leaveCentralElementOfEnclosingSpacetree() override;

    /**
     * Nop.
     */
    virtual void dataWasNotReceivedInBackground( int fromRank, int tag, int cardinality, int pageSize ) override;

    /**
     * Nop.
     */
    virtual void beginToReleaseSynchronousHeapData() override;

    /**
     * Nop.
     */
    virtual void endToReleaseSynchronousHeapData() override;

    /**
     * Nop.
     */
    virtual void beginToPrepareAsynchronousHeapDataExchange() override;

    /**
     * Nop.
     */
    virtual void endToPrepareAsynchronousHeapDataExchange() override;

    /**
     * Nop.
     */
    virtual void beginReleaseOfJoinData() override;

    /**
     * Nop.
     */
    virtual void endReleaseOfJoinData() override;

    /**
     * Nop.
     */
    virtual void beginReleaseOfBoundaryData() override;

    /**
     * Nop.
     */
    virtual void endReleaseOfBoundaryData() override;

    /**
     * Nop.
     */
    virtual void changeConcurrencyLevel(int actualChange, int maxPossibleChange) override;

    /**
     * Nop.
     */
    virtual void minuteNumberOfBackgroundTasks(int taskCount) override;

};


#endif
