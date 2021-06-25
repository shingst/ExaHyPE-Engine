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

#if  defined(Parallel)
#include "exahype/reactive/PerformanceMonitor.h"

#include <algorithm>
#include <numeric>

#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"
#include "tarch/timing/Watch.h"

#include "exahype/reactive/OffloadingProfiler.h"
#include "exahype/reactive/StaticDistributor.h"
#include "exahype/reactive/ReactiveContext.h"


#define TERMINATE_SIGNAL -1.0

tarch::logging::Log exahype::reactive::PerformanceMonitor::_log( "exahype::reactive::PerformanceMonitor" );

exahype::reactive::PerformanceMonitor::PerformanceMonitor() :
    _isRankActive(true),
	  _isDisabled(false),
    _terminatedGlobally(false),
    _currentTasksLocal(0),
    _remainingTasks(0),
    _tasksPerTimestep(0),
    _fusedGatherRequest(MPI_REQUEST_NULL)
 {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _currentTasksSnapshot = new int[nnodes];

  _currentWaitingTimes = new double[nnodes];
  _currentWaitingTimesSnapshot = new double[nnodes*nnodes];

  _currentBlacklist = new double[nnodes];
  _currentBlacklistSnapshot = new double[nnodes];

  _currentFusedDataSendBuffer = new double[2*nnodes+2];
  _currentFusedDataReceiveBuffer = new double[nnodes*(2*nnodes+2)];

  std::fill(_currentTasksSnapshot, _currentTasksSnapshot+nnodes, 0);

  std::fill(_currentWaitingTimesSnapshot, _currentWaitingTimesSnapshot+nnodes*nnodes, 0);
  std::fill(_currentWaitingTimes, _currentWaitingTimes+nnodes, 0);

  std::fill(_currentBlacklistSnapshot, _currentBlacklistSnapshot+nnodes, 0);
  std::fill(_currentBlacklist, _currentBlacklist+nnodes, 0);

  std::fill(_currentFusedDataSendBuffer, _currentFusedDataSendBuffer+nnodes+nnodes+1, 0);
  std::fill(_currentFusedDataReceiveBuffer, _currentFusedDataReceiveBuffer+nnodes*(nnodes+nnodes+1), 0);
}

exahype::reactive::PerformanceMonitor::~PerformanceMonitor() {

  delete[] _currentWaitingTimesSnapshot;
  delete[] _currentWaitingTimes;

  delete[] _currentFusedDataSendBuffer;
  delete[] _currentFusedDataReceiveBuffer;

  delete[] _currentBlacklist;
  delete[] _currentBlacklistSnapshot;
 
  delete[] _currentTasksSnapshot;
}

void exahype::reactive::PerformanceMonitor::submitWaitingTimeForRank(double waitingTime, int rank) {
  if(waitingTime>0)
    _currentWaitingTimes[rank] =  waitingTime;
  //logInfo("submitWaitingTimes", "submitting new waiting time "<<waitingTime<< " for rank "<<rank); 
}

const double *exahype::reactive::PerformanceMonitor::getWaitingTimesSnapshot() {
  return _currentWaitingTimesSnapshot;
}

void exahype::reactive::PerformanceMonitor::submitBlacklistValueForRank(double bval, int rank) {
  //logInfo("submitBlacklistValue", "new value "<<bval<<" for "<<rank);
  _currentBlacklist[rank] = bval;
}

const double *exahype::reactive::PerformanceMonitor::getBlacklistSnapshot() {
  //int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  //for(int j=0; j<nnodes; j++) 
  //  logInfo("getBlacklistSnapshot()"," val "<<_currentBlacklistSnapshot[j]<< " for "<< j);

  return _currentBlacklistSnapshot;
}

void exahype::reactive::PerformanceMonitor::setTasksPerTimestep(int load) {
  //logInfo("setLocalLoadPerTimestep", "setting local load per timestep to "<<load);
  _tasksPerTimestep = load;
  _remainingTasks = _tasksPerTimestep;
}

int exahype::reactive::PerformanceMonitor::getTasksPerTimestep() {
  return _tasksPerTimestep;
}

int exahype::reactive::PerformanceMonitor::getRemainingTasks() {
  return _remainingTasks;
}

const int* exahype::reactive::PerformanceMonitor::getCurrentTasksSnapshot() {
  return _currentTasksSnapshot;
}

exahype::reactive::PerformanceMonitor& exahype::reactive::PerformanceMonitor::getInstance() {
  static PerformanceMonitor perfMon;
  return perfMon;
}

void exahype::reactive::PerformanceMonitor::stop() {
  _isRankActive=false;
}

void exahype::reactive::PerformanceMonitor::setCurrentTasks(int num) {
  //logInfo("performance monitor", "setting current load to "<<num);
  tarch::multicore::Lock lock(_semaphore);
  _currentTasksLocal = num;
  lock.free();
}

void exahype::reactive::PerformanceMonitor::incCurrentTasks() {
  assertion(_currentTasksLocal>=0);
  _currentTasksLocal++;
}

void exahype::reactive::PerformanceMonitor::decCurrentTasks() {
  _currentTasksLocal--;
  assertion(_currentTasksLocal>=0);
}

void exahype::reactive::PerformanceMonitor::decRemainingTasks() {
  tarch::multicore::Lock lock(_semaphore);
  _remainingTasks--;
  if(_remainingTasks==0) {
    _remainingTasks=_tasksPerTimestep;
  }
  lock.free();
  assertion(_remainingTasks>=0);
}

void exahype::reactive::PerformanceMonitor::disable() {
	_isDisabled = true;
}

void exahype::reactive::PerformanceMonitor::run() {
  if(!_isDisabled)
    progressGather();
}

void exahype::reactive::PerformanceMonitor::progressGather() {
  int nnodes    = tarch::parallel::Node::getInstance().getNumberOfNodes();

  int completed_fused = 0;

  tarch::multicore::Lock lock(_semaphore);

  if( !isGloballyTerminated() && _fusedGatherRequest!=MPI_REQUEST_NULL) {
    int err = MPI_Test(&_fusedGatherRequest, &completed_fused, MPI_STATUS_IGNORE);
    assert(err==MPI_SUCCESS);
  }
 
  if(completed_fused) {
    bool newGlobalTerminationStatus = true;
    double *newSnapshot = new double[nnodes];
    std::fill(&newSnapshot[0], &newSnapshot[nnodes], 0);
    //logInfo("progressGather", " got new fused result" );
    for(int i=0; i<nnodes; i++) {
      //copy waiting times
      int offsetWaitingTimes = i*(nnodes+nnodes+2);
      std::copy(&_currentFusedDataReceiveBuffer[offsetWaitingTimes], &_currentFusedDataReceiveBuffer[offsetWaitingTimes+nnodes], &_currentWaitingTimesSnapshot[i*nnodes]);
      int offsetBlacklistValues = i*(nnodes+nnodes+2)+nnodes;

      //reduce blacklist values using a sum reduction
      for(int j=0; j<nnodes; j++) {
      //   if(_currentFusedDataReceiveBuffer[offsetBlacklistValues+j]>0)
      //     logInfo("reduceBVal()"," val "<<_currentFusedDataReceiveBuffer[offsetBlacklistValues+j]<< " for "<< j);
        newSnapshot[j] += _currentFusedDataReceiveBuffer[offsetBlacklistValues+j];
      }

      _currentTasksSnapshot[i] = _currentFusedDataReceiveBuffer[offsetBlacklistValues+nnodes];
      //reduce termination status
      newGlobalTerminationStatus &= (_currentFusedDataReceiveBuffer[offsetBlacklistValues+nnodes+1]==TERMINATE_SIGNAL);
    }
    _terminatedGlobally = newGlobalTerminationStatus;
    if(_terminatedGlobally)
      logInfo("progressOffloading", "received terminated"<<_terminatedGlobally);

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

void exahype::reactive::PerformanceMonitor::postFusedRequest() {
  logDebug("postFusedRequest()", "performance monitor posted fused request");

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::copy(&_currentWaitingTimes[0], &_currentWaitingTimes[nnodes], &_currentFusedDataSendBuffer[0]);
  std::copy(&_currentBlacklist[0], &_currentBlacklist[nnodes], &_currentFusedDataSendBuffer[nnodes]);
  _currentFusedDataSendBuffer[2*nnodes]= static_cast<double> (_currentTasksLocal.load());
  _currentFusedDataSendBuffer[2*nnodes+1]= !_isRankActive ? TERMINATE_SIGNAL : 0; //0 means running

  assertion(_fusedGatherRequest==MPI_REQUEST_NULL);
 
  //for(int i=0; i< nnodes*2;i++)
  //  logInfo("postFusedrequest()", "send buffer "<<_currentFusedDataSendBuffer[i]);
  MPI_Comm comm =  exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(); 
  
  if(comm!=MPI_COMM_NULL) {
    int err = MPI_Iallgather(&_currentFusedDataSendBuffer[0], 2*nnodes+2, MPI_DOUBLE, &_currentFusedDataReceiveBuffer[0],
                   2*nnodes+2, MPI_DOUBLE, comm,
                   &_fusedGatherRequest); assertion(err==MPI_SUCCESS);
    assert(err==MPI_SUCCESS);
  }
}

bool exahype::reactive::PerformanceMonitor::isGloballyTerminated() {
  return _terminatedGlobally;
}

#endif
