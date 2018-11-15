// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_STEALING_STEALINGANALYSER_H_
#define _EXAHYPE_STEALING_STEALINGANALYSER_H_

#include "peano/performanceanalysis/Analyser.h"

#include "tarch/logging/Log.h"
#include "tarch/timing/Watch.h"
#include "tarch/timing/GlidingAverageMeasurement.h"

#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/Node.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include <vector>


namespace exahype {
  namespace stealing {
    class StealingAnalyser;
  }
}

class exahype::stealing::StealingAnalyser: public peano::performanceanalysis::Analyser {
  private:
    static tarch::logging::Log     _log;

    bool _isSwitchedOn;


    tarch::timing::Watch           _waitForWorkerDataWatch;
    tarch::timing::Watch           _waitForMasterDataWatch;

    tarch::timing::GlidingAverageMeasurement     _traversalMeasurement;
    std::vector<tarch::timing::GlidingAverageMeasurement>    _waitForOtherRank;

    double _currentMaxWaitTime;

  public:
    StealingAnalyser();
    virtual ~StealingAnalyser();

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

    virtual void beginToReceiveDataFromWorker();
    virtual void endToReceiveDataFromWorker( int fromRank );
    virtual void beginToReceiveDataFromMaster();
    virtual void endToReceiveDataFromMaster();

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
