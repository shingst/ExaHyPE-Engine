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

    std::vector<tarch::timing::GlidingAverageMeasurement>    _waitForOtherRank;

    tarch::timing::GlidingAverageMeasurement     _timePerMigratablePredictionJob;

    tarch::timing::GlidingAverageMeasurement     _timePerTimeStep;

    double _currentZeroThreshold;
    int _iterationCounter;
    double _currentAccumulatedWorkerTime;

    double _estWtimeForPendingJobs;

    double *_currentFilteredWaitingTimesSnapshot;

    std::atomic<int> _lateSTPJobs;

    void updateZeroTresholdAndFilteredSnapshot();

  public:
    OffloadingAnalyser();
    virtual ~OffloadingAnalyser();

    OffloadingAnalyser(const OffloadingAnalyser& other) = delete;
    OffloadingAnalyser& operator=(const OffloadingAnalyser& other) = delete;

    const double* getFilteredWaitingTimesSnapshot();

    void notifyReceivedSTPJob();    

    void printWaitingTimes();

    static OffloadingAnalyser& getInstance();

    double getZeroThreshold();

    void setTimePerSTP(double timePerSTP);
    double getTimePerSTP();

    void setTimePerTimeStep(double timePerStep);
    double getTimePerTimeStep();

    void resetMeasurements();

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
