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

#if  defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#include "PerformanceMonitor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/stealing/StealingProfiler.h"
#include "exahype/stealing/DynamicDistributor.h"
#include "exahype/stealing/StaticDistributor.h"
#include "exahype/stealing/StealingManager.h"


#define TERMINATE_SIGNAL -1.0

tarch::logging::Log exahype::stealing::PerformanceMonitor::_log( "exahype::stealing::PerformanceMonitor" );

exahype::stealing::PerformanceMonitor::PerformanceMonitor() :
    _isStarted(true),
    _gatherTasksRequest(MPI_REQUEST_NULL),
    _gatherWaitingTimesRequest(MPI_REQUEST_NULL),
    _allreduceBlacklistRequest(MPI_REQUEST_NULL),
    _fusedGatherRequest(MPI_REQUEST_NULL),
    _currentTasks(0),
    _currentTasksSendBuffer(0),
    _remainingTasks(0),
    _tasksPerTimestep(0),
    _terminatedGlobally(false) {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //_currentTasksSnapshot = new int[nnodes];
  //_currentTasksReceiveBuffer   = new int[nnodes];

  _currentWaitingTimesSnapshot = new double[nnodes*nnodes];
  //_currentWaitingTimesReceiveBuffer = new int[nnodes*nnodes];
  //_currentWaitingTimesSendBuffer = new int[nnodes];
  _currentWaitingTimes = new double[nnodes];

  _currentBlacklistSnapshot = new double[nnodes];
  //_currentBlacklistReceiveBuffer = new double[nnodes*nnodes];
  //_currentBlacklistSendBuffer = new double[nnodes];
  _currentBlacklist = new double[nnodes];

  _currentFusedDataSendBuffer = new double[nnodes+nnodes+1];
  _currentFusedDataReceiveBuffer = new double[nnodes*(nnodes+nnodes+1)];

  //std::fill(_currentTasksSnapshot, _currentTasksSnapshot+nnodes, 0);
  //std::fill(_currentTasksReceiveBuffer, _currentTasksReceiveBuffer+nnodes, 0);

  std::fill(_currentWaitingTimesSnapshot, _currentWaitingTimesSnapshot+nnodes*nnodes, 0);
  //std::fill(_currentWaitingTimesReceiveBuffer, _currentWaitingTimesReceiveBuffer+nnodes*nnodes, 0);
  //std::fill(_currentWaitingTimesSendBuffer, _currentWaitingTimesSendBuffer+nnodes, 0);
  std::fill(_currentWaitingTimes, _currentWaitingTimes+nnodes, 0);

  std::fill(_currentBlacklistSnapshot, _currentBlacklistSnapshot+nnodes, 0);
  //std::fill(_currentBlacklistReceiveBuffer, _currentBlacklistReceiveBuffer+nnodes, 0);
  //std::fill(_currentBlacklistSendBuffer, _currentBlacklistSendBuffer+nnodes, 0);
  std::fill(_currentBlacklist, _currentBlacklist+nnodes, 0);

  std::fill(_currentFusedDataSendBuffer, _currentFusedDataSendBuffer+nnodes+nnodes+1, 0);
  std::fill(_currentFusedDataReceiveBuffer, _currentFusedDataReceiveBuffer+nnodes*(nnodes+nnodes+1), 0);
}

exahype::stealing::PerformanceMonitor::~PerformanceMonitor() {
  //delete[] _currentWaitingTimesSendBuffer;
  //delete[] _currentWaitingTimesReceiveBuffer;
  delete[] _currentWaitingTimesSnapshot;
  delete[] _currentWaitingTimes;

  delete[] _currentBlacklist;
  delete[] _currentBlacklistSnapshot;

  delete[] _currentFusedDataSendBuffer;
  delete[] _currentFusedDataReceiveBuffer;
 
  //delete[] _currentTasksSnapshot;
  //delete[] _currentTasksReceiveBuffer;
}

void exahype::stealing::PerformanceMonitor::submitWaitingTimeForRank(double waitingTime, int rank) {
  if(waitingTime>0)
    _currentWaitingTimes[rank] =  waitingTime;
    //logInfo("submitWaitingTimes", "submitting new waiting time "<<waitingTime<< " for rank "<<rank); 
}

const double *exahype::stealing::PerformanceMonitor::getWaitingTimesSnapshot() {
  return _currentWaitingTimesSnapshot;
}

void exahype::stealing::PerformanceMonitor::submitBlacklistValueForRank(double bval, int rank) {

  //logInfo("submitBlacklistValue", "new value "<<bval<<" for "<<rank);
  _currentBlacklist[rank] = bval;
}

const double *exahype::stealing::PerformanceMonitor::getBlacklistSnapshot() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //for(int j=0; j<nnodes; j++) 
  //  logInfo("getBlacklistSnapshot()"," val "<<_currentBlacklistSnapshot[j]<< " for "<< j);

  return _currentBlacklistSnapshot;
}

void exahype::stealing::PerformanceMonitor::setTasksPerTimestep(int load) {
  logInfo("setLocalLoadPerTimestep", "setting local load per timestep to "<<load);
  _tasksPerTimestep = load;
  _remainingTasks = _tasksPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getTasksPerTimestep() {
  return _tasksPerTimestep;
}

int exahype::stealing::PerformanceMonitor::getRemainingTasks() {
  return _remainingTasks;
}

const int* exahype::stealing::PerformanceMonitor::getCurrentTasksSnapshot() {
  return _currentTasksSnapshot;
}

exahype::stealing::PerformanceMonitor& exahype::stealing::PerformanceMonitor::getInstance() {
  static PerformanceMonitor perfMon;
  return perfMon;
}

void exahype::stealing::PerformanceMonitor::stop() {
    _isStarted=false;
}

void exahype::stealing::PerformanceMonitor::setCurrentTasks(int num) {
  //logInfo("performance monitor", "setting current load to "<<num);
  int myRank = tarch::parallel::Node::getInstance().getRank();
  tarch::multicore::Lock lock(_semaphore);
  _currentTasksSnapshot[myRank] = num;
  _currentTasks = num;
  lock.free();
}

void exahype::stealing::PerformanceMonitor::incCurrentTasks() {
#ifndef StealingStrategyDiffusive
  assertion(_currentTasks>=0);
  _currentTasks++;
#endif
}

void exahype::stealing::PerformanceMonitor::decCurrentTasks() {
#ifndef StealingStrategyDiffusive
  _currentTasks--;
  assertion(_currentTasks>=0);
#endif
}

void exahype::stealing::PerformanceMonitor::decRemainingTasks() {
#ifndef StealingStrategyDiffusive
  tarch::multicore::Lock lock(_semaphore);
  _remainingTasks--;
  if(_remainingTasks==0) {
    _remainingTasks=_tasksPerTimestep;
  }
  lock.free();
  assertion(_remainingTasks>=0);
#endif
}

void exahype::stealing::PerformanceMonitor::run() {
  progressGather();
}

void exahype::stealing::PerformanceMonitor::progressGather() {
  int myRank    = tarch::parallel::Node::getInstance().getRank();
  int nnodes    = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //int completed_tasks = 0;
  //int completed_waiting_times = 0;
  //int completed_blacklist = 0;

  int completed_fused = 0;
//#if defined(PerformanceAnalysisStealing)
//  double timeSinceLastGather=0;
//  static std::atomic<double> lastGather = 0;
//  static std::atomic<int> unsuccessful  = 0;
//  static std::atomic<int> successful    = 0;

//  tarch::timing::Watch watch("exahype::stealing::", "-", false,false);
//  watch.startTimer();

//  timeSinceLastGather = lastGather + MPI_Wtime();
//#endif

  tarch::multicore::Lock lock(_semaphore);

//  if( !isGloballyTerminated() && _gatherTasksRequest!=MPI_REQUEST_NULL) {
//	double time = - MPI_Wtime();
//    exahype::stealing::StealingProfiler::getInstance().beginCommunication();
//    int err = MPI_Test(&_gatherTasksRequest, &completed_tasks, MPI_STATUS_IGNORE); //assert(err==MPI_SUCCESS);
//    time += MPI_Wtime();

//    //std::string str;
//    //for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_currentTasksReceiveBuffer[i]);
//    //str+="\n";
//    //logInfo("performance monitor", "progressing "<<_gatherTasksRequest);
//    

//#if defined(PerformanceAnalysisStealing)
//    if(completed) {
//      exahype::stealing::StealingProfiler::getInstance().endCommunication(true, time);
//      successful++;
//    }
//    else {
//      exahype::stealing::StealingProfiler::getInstance().endCommunication(false, time);
//      unsuccessful++;
//    }
//    if(successful%1000==0 || unsuccessful%10000==0) {
//      logInfo("performance monitor", " successful "<<successful<<" unsuccessful "<<unsuccessful<<" ratio "<<1.0f*successful/unsuccessful);
//    }
//#endif
//  }

//  if( !isGloballyTerminated() && _gatherWaitingTimesRequest!=MPI_REQUEST_NULL) {
////    double time = - MPI_Wtime();
////    exahype::stealing::StealingProfiler::getInstance().beginCommunication();
//    //logInfo("performance monitoŕ()","progressing waiting times");
//    int err= MPI_Test(&_gatherWaitingTimesRequest, &completed_waiting_times, MPI_STATUS_IGNORE);
//   //assert(err==MPI_SUCCESS);
////    time += MPI_Wtime();_currentTasksSendBuffer
// 
//  }

//  if( !isGloballyTerminated() && _allreduceBlacklistRequest!=MPI_REQUEST_NULL) {
//     int err= MPI_Test(&_allreduceBlacklistRequest, &completed_blacklist, MPI_STATUS_IGNORE);
//  }

//#if defined(PerformanceAnalysisStealing)
//  watch.stopTimer();
//  if(watch.getCalendarTime() >= 0.00001) {
//    logInfo(
//        "performance monitor()",
//        "MPI " <<
//        "time=" << std::fixed <<
//        watch.getCalendarTime() <<
//        ", cpu time=" <<
//        watch.getCPUTime()
//    );
//  }
//  watch.startTimer();
//#endif

//  if(completed_tasks) {
//  //  logInfo("progressGather","collected new tasks snapshot");
//    stealing::StealingProfiler::getInstance().notifyPerformanceUpdate();
//    std::copy(&_currentTasksReceiveBuffer[0], &_currentTasksReceiveBuffer[nnodes], &_currentTasksSnapshot[0]);
// 
//    //std::string str;
//    //for(int i=0;i<nnodes;i++) str=str+" , "+std::to_string(_currentTasksSnapshot[i]);
//    //str+="\n";
//    //logInfo("performance monitor", str);

//#if defined(PerformanceAnalysisStealing)
//    std::string str="received new update, current load "+std::to_string(_currentTasks.load());
////    if(timeSinceLastGather>0.001) {
////      str=str+ " took too long: "+std::to_string(timeSinceLastGather);
////      stealing::StealingProfiler::getInstance().notifyLatePerformanceUpdate();
////    }

//#endif
//    if(_currentTasks.load()>0 && std::all_of(&_currentTasksSnapshot[0], &_currentTasksSnapshot[nnodes], [](int i) {return i>=0;})) {
//      exahype::stealing::DynamicDistributor::getInstance().computeNewLoadDistribution(_currentTasksSnapshot);
//      stealing::StealingProfiler::getInstance().notifyStealingDecision();
//    }
//#if defined(PerformanceAnalysisStealing)
//    watch.stopTimer();
//    if(watch.getCalendarTime() >= 0.00001) {
//      logInfo(
//          "performance monitor()",
//          "completion" <<
//          "time=" << std::fixed <<
//          watch.getCalendarTime() <<
//          ", cpu time=" <<
//          watch.getCPUTime()
//      );
//    }
//    watch.startTimer();
//#endif
//    _gatherTasksRequest = MPI_REQUEST_NULL;
//  }

//  if(completed_waiting_times) {
//    //logInfo("progressGather","collected new waiting times snapshot, request"<< _gatherWaitingTimesRequest);
//    std::copy(&_currentWaitingTimesReceiveBuffer[0], &_currentWaitingTimesReceiveBuffer[nnodes*nnodes], &_currentWaitingTimesSnapshot[0]);
//    _gatherWaitingTimesRequest = MPI_REQUEST_NULL;
//      
//    //int k = 0;
//    //for(int i=0; i<nnodes; i++) {
//    //  for(int j=0; j<nnodes; j++) {
//    //    logInfo("progressGather()","rank "<<i<<" waiting for "<<_currentWaitingTimesSnapshot[k+j]<<" for rank "<<j);
//    //  }
//    //  k+= nnodes;
//    //}
//  }

//  if(_gatherWaitingTimesRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
//    //logInfo("progressGather","post gather waiting times");
//    postGatherWaitingTimes();
//  }

//  if(_gatherTasksRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
//    //logInfo("progressGather","post gather tasks");
//    //postGatherTasks();
//#if defined(PerformanceAnalysisStealing)
//    lastGather=-MPI_Wtime();
//    watch.stopTimer();
//    if(watch.getCalendarTime() >= 0.00001) {
//      logInfo(
//          "performance monitor()",
//          "post gather" <<
//          "time=" << std::fixed <<
//          watch.getCalendarTime() <<
//          ", cpu time=" <<
//          watch.getCPUTime()
//      );
//    }
//    watch.startTimer();
//#endif
//  }

//  if(completed_blacklist) {
//    std::copy(&_currentBlacklistReceiveBuffer[0], &_currentBlacklistReceiveBuffer[nnodes], &_currentBlacklistSnapshot[0]);
//    _allreduceBlacklistRequest = MPI_REQUEST_NULL;
//      
//  }

//  if(_allreduceBlacklistRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
//    postAllreduceBlacklist();
//  }

  if( !isGloballyTerminated() && _fusedGatherRequest!=MPI_REQUEST_NULL) {
    int err = MPI_Test(&_fusedGatherRequest, &completed_fused, MPI_STATUS_IGNORE);
  }
 
  if(completed_fused) {
    bool newGlobalTerminationStatus = true;
    double *newSnapshot = new double[nnodes];
    std::fill(&newSnapshot[0], &newSnapshot[nnodes], 0);
    //logInfo("progressGather", " got new fused result" );
    for(int i=0; i<nnodes; i++) {
       //copy waiting times
       int offsetWaitingTimes = i*(nnodes+nnodes+1);
       std::copy(&_currentFusedDataReceiveBuffer[offsetWaitingTimes], &_currentFusedDataReceiveBuffer[offsetWaitingTimes+nnodes], &_currentWaitingTimesSnapshot[i*nnodes]);
       int offsetBlacklistValues = i*(nnodes+nnodes+1)+nnodes;

       //reduce blacklist values
       for(int j=0; j<nnodes; j++) {
      //   if(_currentFusedDataReceiveBuffer[offsetBlacklistValues+j]>0)
      //     logInfo("reduceBVal()"," val "<<_currentFusedDataReceiveBuffer[offsetBlacklistValues+j]<< " for "<< j);
         newSnapshot[j] += _currentFusedDataReceiveBuffer[offsetBlacklistValues+j];
       }

       //reduce termination status
       newGlobalTerminationStatus &= (_currentFusedDataReceiveBuffer[offsetBlacklistValues+nnodes]==TERMINATE_SIGNAL);
    }
    _terminatedGlobally = newGlobalTerminationStatus;
    if(_terminatedGlobally)
     logInfo("progressStealing", "received terminated"<<_terminatedGlobally);

 
    for(int j=0; j<nnodes; j++) {
       _currentBlacklistSnapshot[j] = newSnapshot[j];
       //logInfo("afterReduction()"," val "<<_currentBlacklistSnapshot[j]<< " for "<< j);
    }

    _fusedGatherRequest = MPI_REQUEST_NULL;
    delete[] newSnapshot;
  }


  if(_fusedGatherRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    postFusedRequest();
  }

  lock.free();
}

//void exahype::stealing::PerformanceMonitor::postGatherTasks() {
//  int myRank = tarch::parallel::Node::getInstance().getRank();

//  _currentTasksSendBuffer = _isStarted ? _currentTasks.load() : TERMINATE_SIGNAL;
//  
//  //if(!_isStarted)
//  //  logInfo("postGatherTasks","posted terminate signal");

//  assert(_gatherTasksRequest==MPI_REQUEST_NULL);
//  int err = MPI_Iallgather(&_currentTasksSendBuffer, 1, MPI_INTEGER, _currentTasksReceiveBuffer, 1, MPI_INTEGER, exahype::stealing::StealingManager::getInstance().getMPICommunicator(), &_gatherTasksRequest); //assert(err==MPI_SUCCESS);
//}

//void exahype::stealing::PerformanceMonitor::postGatherWaitingTimes() {
//  int nnodes    = tarch::parallel::Node::getInstance().getNumberOfNodes();
//  std::copy(&_currentWaitingTimes[0], &_currentWaitingTimes[nnodes], &_currentWaitingTimesSendBuffer[0]);


// // static int nposted = 0;
//  //logInfo("postGatherWaitingTimes","posting gather for waiting times");
//  //for(int i=0; i<nnodes; i++) {
//  //  logInfo("postGatherWaitingTimes","_currentWaitingTimesSendBuffer["<<i<<"]: "<<_currentWaitingTimesSendBuffer[i]);
//  //}
//  assert(_gatherWaitingTimesRequest==MPI_REQUEST_NULL);
//  int err = MPI_Iallgather(&_currentWaitingTimesSendBuffer[0], nnodes, MPI_INTEGER, &_currentWaitingTimesReceiveBuffer[0],
//                   nnodes, MPI_INTEGER, exahype::stealing::StealingManager::getInstance().getMPICommunicator(),
//                   &_gatherWaitingTimesRequest);// assert(err==MPI_SUCCESS);

//  //nposted++;
//  //logInfo("postGatherWaitingTimes", "nposted "<<nposted<< " request "<<_gatherWaitingTimesRequest);
//}

//void exahype::stealing::PerformanceMonitor::postAllreduceBlacklist() {
//  int nnodes  = tarch::parallel::Node::getInstance().getNumberOfNodes();
//  std::copy(&_currentBlacklist[0], &_currentBlacklist[nnodes], &_currentBlacklistSendBuffer[0]);

//  assert(_allreduceBlacklistRequest==MPI_REQUEST_NULL);
//  int err = MPI_Iallreduce(&_currentBlacklistSendBuffer[0], &_currentBlacklistReceiveBuffer[0], nnodes, MPI_DOUBLE, MPI_SUM,
//             exahype::stealing::StealingManager::getInstance().getMPICommunicator(), &_allreduceBlacklistRequest);
//}

void exahype::stealing::PerformanceMonitor::postFusedRequest() {
  //logInfo("postFusedRequest()", "performance monitor posted fused request");
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  std::copy(&_currentWaitingTimes[0], &_currentWaitingTimes[nnodes], &_currentFusedDataSendBuffer[0]);
  std::copy(&_currentBlacklist[0], &_currentBlacklist[nnodes], &_currentFusedDataSendBuffer[nnodes]);
  _currentFusedDataSendBuffer[2*nnodes]= !_isStarted ? TERMINATE_SIGNAL : 0; //0 means running   

  assert(_fusedGatherRequest==MPI_REQUEST_NULL);
 
  //for(int i=0; i< nnodes*2;i++)
  //  logInfo("postFusedrequest()", "send buffer "<<_currentFusedDataSendBuffer[i]);

  int err = MPI_Iallgather(&_currentFusedDataSendBuffer[0], 2*nnodes+1, MPI_DOUBLE, &_currentFusedDataReceiveBuffer[0],
                   2*nnodes+1, MPI_DOUBLE, exahype::stealing::StealingManager::getInstance().getMPICommunicator(),
                   &_fusedGatherRequest);// assert(err==MPI_SUCCESS);

}

bool exahype::stealing::PerformanceMonitor::isGloballyTerminated() {
/*  static bool globalTermination=false;
  if(globalTermination) return true;
  if(_isStarted) return false;

  bool result=true;
  for(int i=0; i<tarch::parallel::Node::getInstance().getNumberOfNodes(); i++)
    result &= _currentTasksReceiveBuffer[i]==TERMINATE_SIGNAL;

  globalTermination=result; */
  return _terminatedGlobally;
}

#endif
