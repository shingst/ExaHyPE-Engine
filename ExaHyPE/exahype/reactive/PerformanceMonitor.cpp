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

#ifndef MPI_CHECK
#ifndef Asserts
#define MPI_CHECK(func, x) do { \
  ierr = (x); \
  if (ierr != MPI_SUCCESS) { \
    logError(#func, "Runtime error:"<<#x<<" returned "<<ierr<<" at " << __FILE__<< ":"<< __LINE__); \
  } \
} while (0)
#else
#define MPI_CHECK(func, x) do { \
  ierr = (x); \
  } while (0)
#endif
#endif

tarch::logging::Log exahype::reactive::PerformanceMonitor::_log( "exahype::reactive::PerformanceMonitor" );

exahype::reactive::PerformanceMonitor::PerformanceMonitor() :
    _hasTerminatedLocally(false),
	  _isDisabled(false),
    _terminatedGlobally(false),
    _currentTasksLocal(0),
    _tasksPerTimestep(0),
    _fusedGatherRequest(MPI_REQUEST_NULL)
 {

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  _currentTasksSnapshot = new int[nnodes];

  _currentWaitingTimesLocal = new double[nnodes];
  _currentWaitingTimesGlobalSnapshot = new double[nnodes*nnodes];

  _currentBlacklistLocal = new double[nnodes];
  _currentBlacklistGlobalSnapshot = new double[nnodes];

  _currentFusedDataSendBuffer = new double[2*nnodes+2];
  _currentFusedDataReceiveBuffer = new double[nnodes*(2*nnodes+2)];

  std::fill(_currentTasksSnapshot, _currentTasksSnapshot+nnodes, 0);

  std::fill(_currentWaitingTimesGlobalSnapshot, _currentWaitingTimesGlobalSnapshot+nnodes*nnodes, 0);
  std::fill(_currentWaitingTimesLocal, _currentWaitingTimesLocal+nnodes, 0);

  std::fill(_currentBlacklistGlobalSnapshot, _currentBlacklistGlobalSnapshot+nnodes, 0);
  std::fill(_currentBlacklistLocal, _currentBlacklistLocal+nnodes, 0);

  std::fill(_currentFusedDataSendBuffer, _currentFusedDataSendBuffer+nnodes+nnodes+1, 0);
  std::fill(_currentFusedDataReceiveBuffer, _currentFusedDataReceiveBuffer+nnodes*(nnodes+nnodes+1), 0);
}

exahype::reactive::PerformanceMonitor::~PerformanceMonitor() {

  delete[] _currentWaitingTimesGlobalSnapshot;
  delete[] _currentWaitingTimesLocal;

  delete[] _currentFusedDataSendBuffer;
  delete[] _currentFusedDataReceiveBuffer;

  delete[] _currentBlacklistLocal;
  delete[] _currentBlacklistGlobalSnapshot;
 
  delete[] _currentTasksSnapshot;
}

void exahype::reactive::PerformanceMonitor::submitWaitingTimeForRank(double waitingTime, int rank) {
  if(waitingTime>0)
    _currentWaitingTimesLocal[rank] =  waitingTime;
}

const double *exahype::reactive::PerformanceMonitor::getWaitingTimesGlobalSnapshot() const {
  return _currentWaitingTimesGlobalSnapshot;
}

void exahype::reactive::PerformanceMonitor::submitBlacklistValueForRank(double bval, int rank) {
  _currentBlacklistLocal[rank] = bval;
}

const double *exahype::reactive::PerformanceMonitor::getBlacklistGlobalSnapshot() const {
  return _currentBlacklistGlobalSnapshot;
}

void exahype::reactive::PerformanceMonitor::setTasksPerTimestep(int load) {
  _tasksPerTimestep = load;
}

int exahype::reactive::PerformanceMonitor::getTasksPerTimestep() const {
  return _tasksPerTimestep;
}

const int* exahype::reactive::PerformanceMonitor::getCurrentTasksGlobalSnapshot() const {
  return _currentTasksSnapshot;
}

exahype::reactive::PerformanceMonitor& exahype::reactive::PerformanceMonitor::getInstance() {
  static PerformanceMonitor perfMon;
  return perfMon;
}

void exahype::reactive::PerformanceMonitor::signalLocalTermination() {
  _hasTerminatedLocally = true;
}

void exahype::reactive::PerformanceMonitor::setCurrentNumTasks(int num) {
  tarch::multicore::Lock lock(_semaphore);
  _currentTasksLocal = num;
  lock.free();
}

void exahype::reactive::PerformanceMonitor::incCurrentNumTasks() {
  assertion(_currentTasksLocal>=0);
  _currentTasksLocal++;
}

void exahype::reactive::PerformanceMonitor::decCurrentNumTasks() {
  _currentTasksLocal--;
  assertion(_currentTasksLocal>=0);
}

void exahype::reactive::PerformanceMonitor::disable() {
	_isDisabled = true;
}

void exahype::reactive::PerformanceMonitor::progress() {
  if(!_isDisabled)
    progressGatherAndCollectNewSnapshot();
}

void exahype::reactive::PerformanceMonitor::progressGatherAndCollectNewSnapshot() {
  int nnodes    = tarch::parallel::Node::getInstance().getNumberOfNodes();

  int reqCompleted = 0;
  int ierr = MPI_SUCCESS;

  tarch::multicore::Lock lock(_semaphore);

  if( !isGloballyTerminated() && _fusedGatherRequest!=MPI_REQUEST_NULL) {
    MPI_CHECK("progressGatherAndCollectNewSnapshot", MPI_Test(&_fusedGatherRequest, &reqCompleted, MPI_STATUS_IGNORE));
  }
 
  if(reqCompleted) {
    bool newGlobalTerminationStatus = true;
    double *newSnapshot = new double[nnodes];
    std::fill(&newSnapshot[0], &newSnapshot[nnodes], 0);

    for(int i=0; i<nnodes; i++) {
      //copy waiting times
      int offsetWaitingTimes = i*(nnodes+nnodes+2);
      std::copy(&_currentFusedDataReceiveBuffer[offsetWaitingTimes], &_currentFusedDataReceiveBuffer[offsetWaitingTimes+nnodes], &_currentWaitingTimesGlobalSnapshot[i*nnodes]);
      int offsetBlacklistValues = i*(nnodes+nnodes+2)+nnodes;

      //reduce blacklist values using a sum reduction
      for(int j=0; j<nnodes; j++) {
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
      _currentBlacklistGlobalSnapshot[j] = newSnapshot[j];
    }

    _fusedGatherRequest = MPI_REQUEST_NULL;
    delete[] newSnapshot;
  }

  if(_fusedGatherRequest==MPI_REQUEST_NULL && !isGloballyTerminated()) {
    copyToSendBufferAndPostFusedRequest();
  }

  lock.free();
}

void exahype::reactive::PerformanceMonitor::copyToSendBufferAndPostFusedRequest() {
  logDebug("postFusedRequest()", "performance monitor posted fused request");

  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  int ierr = MPI_SUCCESS;

  std::copy(&_currentWaitingTimesLocal[0], &_currentWaitingTimesLocal[nnodes], &_currentFusedDataSendBuffer[0]);
  std::copy(&_currentBlacklistLocal[0], &_currentBlacklistLocal[nnodes], &_currentFusedDataSendBuffer[nnodes]);
  _currentFusedDataSendBuffer[2*nnodes]= static_cast<double> (_currentTasksLocal.load());
  _currentFusedDataSendBuffer[2*nnodes+1]= _hasTerminatedLocally ? TERMINATE_SIGNAL : 0; //0 means running

  assertion(_fusedGatherRequest==MPI_REQUEST_NULL);
 
  MPI_Comm comm =  exahype::reactive::ReactiveContext::getInstance().getMPICommunicator(); 
  
  if(comm!=MPI_COMM_NULL) {
    MPI_CHECK("copyToSendBufferAndPostFusedRequest", MPI_Iallgather(&_currentFusedDataSendBuffer[0], 2*nnodes+2, MPI_DOUBLE, &_currentFusedDataReceiveBuffer[0],
                   2*nnodes+2, MPI_DOUBLE, comm,
                   &_fusedGatherRequest));
  }
}

bool exahype::reactive::PerformanceMonitor::isGloballyTerminated() {
  return _terminatedGlobally;
}

#endif
