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
  namespace offloading {
    class OffloadingAnalyser;
  }
}

class exahype::offloading::OffloadingAnalyser : public peano::performanceanalysis::Analyser {
  private:
    static tarch::logging::Log     _log;

    bool _isSwitchedOn;

    tarch::timing::Watch           _waitForWorkerDataWatch;
    tarch::timing::Watch           _waitForMasterDataWatch;
    tarch::timing::Watch           _waitForGlobalMasterDataWatch;

    tarch::timing::Watch           _timeStepWatch;

    std::vector<tarch::timing::GlidingAverageMeasurement>    _waitForOtherRank;

    tarch::timing::GlidingAverageMeasurement     _timePerMigratablePredictionJob;

    tarch::timing::GlidingAverageMeasurement     _timePerTimeStep;

    double _currentZeroThreshold;
    int _iterationCounter;
    double _currentAccumulatedWorkerTime;

    double _estimatedWtimeForPendingJobs;

    double *_currentFilteredWaitingTimesSnapshot;

    std::atomic<int> _lateSTPJobs;

    void updateZeroTresholdAndFilteredSnapshot();

  public:
    OffloadingAnalyser();
    virtual ~OffloadingAnalyser();

    const double* getFilteredWaitingTimesSnapshot();

    void notifyReceivedSTPJob();    

    void printWaitingTimes();

    static OffloadingAnalyser& getInstance();

    double getZeroThreshold();

    void setTimePerSTP(double timePerSTP);
    double getTimePerSTP();

    void setTimePerTimeStep(double timePerStep);
    double getTimePerTimeStep();

    virtual void beginIteration();
    virtual void endIteration(double numberOfInnerLeafCells, double numberOfOuterLeafCells, double numberOfInnerCells, double numberOfOuterCells, double numberOfLocalCells, double numberOfLocalVertices);

    virtual void enterCentralElementOfEnclosingSpacetree();
    virtual void leaveCentralElementOfEnclosingSpacetree();

    virtual void addWorker(
      int                                 workerRank,
      int                                 level
    );

    virtual void removeWorker(
      int                                 workerRank,
      int                                 level
    );


    virtual void beginToSendDataToWorker();
    virtual void endToSendDataToWorker(int worker);
    virtual void beginToSendDataToMaster();
    virtual void endToSendDataToMaster();
    virtual void beginToReceiveDataFromWorker();
    virtual void endToReceiveDataFromWorker(int fromRank);
    virtual void beginToReceiveDataFromMaster(int master);
    virtual void endToReceiveDataFromMaster(int master);
    virtual void beginToReceiveDataFromGlobalMaster();
    virtual void endToReceiveDataFromGlobalMaster();

    virtual void dataWasNotReceivedInBackground( int fromRank, int tag, int cardinality, int pageSize );

    virtual void beginToReleaseSynchronousHeapData();

    virtual void endToReleaseSynchronousHeapData();

    virtual void beginToPrepareAsynchronousHeapDataExchange();

    virtual void endToPrepareAsynchronousHeapDataExchange();

    virtual void beginReleaseOfJoinData();
    virtual void endReleaseOfJoinData();

    virtual void beginReleaseOfBoundaryData();
    virtual void endReleaseOfBoundaryData();

    virtual void changeConcurrencyLevel(int actualChange, int maxPossibleChange);
    virtual void minuteNumberOfBackgroundTasks(int taskCount);

    virtual void beginProcessingBackgroundJobs();
    virtual void endProcessingBackgroundJobs();

    virtual void enable(bool value);
};


#endif
