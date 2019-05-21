// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
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
  namespace stealing {
    class StealingAnalyser;
  }
}

class exahype::stealing::StealingAnalyser : public peano::performanceanalysis::Analyser {
  private:
    static tarch::logging::Log     _log;

    bool _isSwitchedOn;


    tarch::timing::Watch           _waitForWorkerDataWatch;
    tarch::timing::Watch           _waitForMasterDataWatch;
    tarch::timing::Watch           _waitForGlobalMasterDataWatch;

    std::vector<tarch::timing::GlidingAverageMeasurement>    _waitForOtherRank;

    tarch::timing::GlidingAverageMeasurement     _timePerStealablePredictionJob;

    double _currentZeroThreshold;
    int _iterationCounter;
    double _currentAccumulatedWorkerTime;

    double _estimatedWtimeForPendingJobs;

    double *_currentFilteredWaitingTimesSnapshot;

    std::atomic<int> _lateSTPJobs;

    void updateZeroTresholdAndFilteredSnapshot();

  public:
    StealingAnalyser();
    virtual ~StealingAnalyser();

    const double* getFilteredWaitingTimesSnapshot();

    void notifyReceivedSTPJob();    

    void printWaitingTimes();

    static StealingAnalyser& getInstance();

    double getZeroThreshold();

    void setTimePerSTP(double timePerSTP);
    double getTimePerSTP();

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
