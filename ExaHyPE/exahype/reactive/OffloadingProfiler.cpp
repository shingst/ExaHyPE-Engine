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


#include "../reactive/OffloadingProfiler.h"

#include "tarch/multicore/Core.h"
#include "tarch/parallel/Node.h"

#if defined(USE_TMPI)
#include "exahype/reactive/ReactiveContext.h"
#endif

tarch::logging::Log exahype::reactive::OffloadingProfiler::_log( "exahype::reactive::OffloadingProfiler" );

exahype::reactive::OffloadingProfiler::OffloadingProfiler():
  _executedTasksPhase(0),
  _executedTasks(0),
  _thresholdFailsPhase(0),
  _thresholdFails(0),
  _spawnedTasksPhase(0),
  _spawnedTasks(0),
  _performanceUpdatesPhase(0),
  _performanceUpdates(0),
  _latePerformanceUpdatesPhase(0),
  _latePerformanceUpdates(0),
  _offloadingDecisionsPhase(0),
  _offloadingDecisions(0),
  _accWaitTasksPhaseTime(0),
  _accUsefulCommunicationPhaseTime(0),
  _accIdleCommunicationPhaseTime(0),
  _accProgressPhaseTime(0),
  _accProgressRequestsPhaseTime(0),
  _accPollPhaseTime(0),
  _accComputationPhaseTime(0),
  _accHandlingPhaseTime(0),
  _accOffloadPhaseTime(0),
  _accWaitForWorkersPhaseTime(0),
  _accWaitForGlobalMasterPhaseTime(0),
  _accWaitTasksTime(0),
  _accUsefulCommunicationTime(0),
  _accIdleCommunicationTime(0),
  _accProgressTime(0),
  _accProgressRequestsTime(0),
  _accPollTime(0),
  _accComputationTime(0),
  _accHandlingTime(0),
  _accOffloadTime(0),
  _accWaitForWorkersTime(0),
  _accWaitForGlobalMasterTime(0)
{
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _offloadedTasksPerRankPhase = new std::atomic<int>[nnodes];
  _offloadedTasksPerRank = new std::atomic<int>[nnodes];
  _targetOffloadedTasksPerRankPhase = new std::atomic<int>[nnodes];
  _targetOffloadedTasksPerRank = new std::atomic<int>[nnodes];
  _receivedTasksPerRankPhase = new std::atomic<int>[nnodes];
  _receivedTasksPerRank = new std::atomic<int>[nnodes];
  std::fill(_offloadedTasksPerRankPhase, _offloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_offloadedTasksPerRank, _offloadedTasksPerRank+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRankPhase, _targetOffloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRank, _targetOffloadedTasksPerRank+nnodes, 0);
  std::fill(_receivedTasksPerRankPhase, _receivedTasksPerRankPhase+nnodes, 0);
  std::fill(_receivedTasksPerRank, _receivedTasksPerRank+nnodes, 0);
}

exahype::reactive::OffloadingProfiler::~OffloadingProfiler() {
  delete[] _offloadedTasksPerRank;
  delete[] _offloadedTasksPerRankPhase;
  delete[] _targetOffloadedTasksPerRank;
  delete[] _targetOffloadedTasksPerRankPhase;
  delete[] _receivedTasksPerRank;
  delete[] _receivedTasksPerRankPhase;
}

exahype::reactive::OffloadingProfiler& exahype::reactive::OffloadingProfiler::getInstance() {
  static exahype::reactive::OffloadingProfiler offloadingProfiler;
  return offloadingProfiler;
}

void exahype::reactive::OffloadingProfiler::beginProfilingPhase() {
#if defined(OffloadingUseProfiler)
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  std::fill(_offloadedTasksPerRankPhase, _offloadedTasksPerRankPhase+nnodes, 0);
  std::fill(_receivedTasksPerRankPhase, _receivedTasksPerRankPhase+nnodes, 0);
  std::fill(_targetOffloadedTasksPerRankPhase, _targetOffloadedTasksPerRankPhase+nnodes, 0);
  _executedTasksPhase              = 0;
  _spawnedTasksPhase               = 0;
  _thresholdFailsPhase             = 0;
  _performanceUpdatesPhase         = 0;
  _latePerformanceUpdatesPhase     = 0;
  _offloadingDecisionsPhase        = 0;
  _accWaitTasksPhaseTime           = 0;
  _accUsefulCommunicationPhaseTime = 0;
  _accIdleCommunicationPhaseTime   = 0;
  _accProgressPhaseTime            = 0;
  _accComputationPhaseTime         = 0;
  _accOffloadPhaseTime             = 0;
  _accWaitForWorkersPhaseTime      = 0;
  _accWaitForGlobalMasterPhaseTime = 0;
#endif
}

void exahype::reactive::OffloadingProfiler::endProfilingPhase() {
#if defined(OffloadingUseProfiler)
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

  for (int i=0;i<nnodes;i++) {
    _offloadedTasksPerRank[i]+=_offloadedTasksPerRankPhase[i];
    _targetOffloadedTasksPerRank[i]+=_targetOffloadedTasksPerRankPhase[i];
    _receivedTasksPerRank[i]+=_receivedTasksPerRankPhase[i];
  }
  _executedTasks+=_executedTasksPhase;
  _spawnedTasks+=_spawnedTasksPhase;
  _thresholdFails+=_thresholdFailsPhase;
  _performanceUpdates+=_performanceUpdatesPhase;
  _latePerformanceUpdates+=_latePerformanceUpdatesPhase;
  _offloadingDecisions+=_offloadingDecisionsPhase;
  _accWaitTasksTime+=_accWaitTasksPhaseTime;
  _accUsefulCommunicationTime+=_accUsefulCommunicationPhaseTime;
  _accIdleCommunicationTime+=_accIdleCommunicationPhaseTime;
  _accProgressTime+=_accProgressPhaseTime;
  _accProgressRequestsTime+=_accProgressRequestsPhaseTime;
  _accPollTime+=_accPollPhaseTime;
  _accComputationTime+=_accComputationPhaseTime;
  _accHandlingTime+=_accHandlingPhaseTime;
  _accOffloadTime+=_accOffloadPhaseTime;
  _accWaitForWorkersTime+=_accWaitForWorkersPhaseTime;
  _accWaitForGlobalMasterTime+=_accWaitForGlobalMasterPhaseTime;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyOffloadedTask(int rank) {
#if defined(OffloadingUseProfiler)
  _offloadedTasksPerRankPhase[rank]++;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyTargetOffloadedTask(int ntasks, int rank) {
#if defined(OffloadingUseProfiler)
  _targetOffloadedTasksPerRankPhase[rank]+=ntasks;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyReceivedTask(int rank) {
#if defined(OffloadingUseProfiler)
  _receivedTasksPerRankPhase[rank]++;
#endif
}

void exahype::reactive::OffloadingProfiler::notifySpawnedTask() {
#if defined(OffloadingUseProfiler)
  _spawnedTasks++;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyPerformanceUpdate() {
#if defined(OffloadingUseProfiler)
  _performanceUpdatesPhase++;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyLatePerformanceUpdate() {
#if defined(OffloadingUseProfiler)
  _latePerformanceUpdatesPhase++;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyOffloadingDecision() {
#if defined(OffloadingUseProfiler)
  _offloadingDecisionsPhase++;
#endif
}

void exahype::reactive::OffloadingProfiler::notifyThresholdFail() {
#if defined(OffloadingUseProfiler)
  _thresholdFailsPhase++;
#endif
}

//Todo: begin could start a threadsafe timer here (which is an implementation difficulty as we don't really have unique
//thread ids with TBB)
void exahype::reactive::OffloadingProfiler::beginComputation() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endComputation(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accComputationPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginProgress() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endProgress(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accProgressPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginProgressRequests() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endProgressRequests(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accProgressRequestsPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginPolling() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endPolling(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accPollPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginCommunication() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endCommunication(bool successful, double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed *1E6;
  if(successful)
    _accUsefulCommunicationPhaseTime+=elapsedTime;
  else
    _accIdleCommunicationPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginWaitForTasks() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endWaitForTasks(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accWaitTasksPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginWaitForWorker() {
#if defined(OffloadingUseProfiler)
#endif
}
    
void exahype::reactive::OffloadingProfiler::endWaitForWorker(double elapsed){
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accWaitForWorkersPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginWaitForGlobalMaster() {
#if defined(OffloadingUseProfiler)
#endif
}
    
void exahype::reactive::OffloadingProfiler::endWaitForGlobalMaster(double elapsed){
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accWaitForGlobalMasterPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::beginOffload() {
#if defined(OffloadingUseProfiler)
#endif
}

void exahype::reactive::OffloadingProfiler::endOffload(double elapsed) {
#if defined(OffloadingUseProfiler)
  const unsigned long long elapsedTime = elapsed*1E6;
  _accOffloadPhaseTime+=elapsedTime;
#endif
}

void exahype::reactive::OffloadingProfiler::printCumulativeStatistics() {
#if defined(OffloadingUseProfiler)
  int rank = tarch::parallel::Node::getInstance().getRank();
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();

#if defined (USE_TMPI)
  int team = exahype::reactive::ReactiveContext::getInstance().getTMPITeamNumber();
#endif

  std::string str="Offloading statistics for rank "+std::to_string(rank);
#if defined (USE_TMPI)
  str += " team = "+ std::to_string(team);
#endif
  str += ":\n";
  logInfo("printCumulativeStatistics", str);
  str="  spawned tasks: "+std::to_string(_spawnedTasks)+"\n";
  logInfo("printCumulativeStatistics", str);
  //str+="  executed tasks: "+std::to_string(_executedTasks)+"\n";
  str="  performance updates: "+std::to_string(_performanceUpdates)+"\n";
  logInfo("printCumulativeStatistics", str);
  //str+="  late performance updates: "+std::to_string(_latePerformanceUpdates)+"\n";
  str="  offloading decisions: "+std::to_string(_offloadingDecisions)+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  target offloaded tasks:\n";
  logInfo("printCumulativeStatistics", str);
  int totalTargetOffloaded=0;
  for(int i=0;i<nnodes;i++) {
    if(_targetOffloadedTasksPerRank[i]>0) {
      str="    to rank: "+std::to_string(i)+" : "+std::to_string(_targetOffloadedTasksPerRank[i])+"\n";
      logInfo("printCumulativeStatistics", str);
      totalTargetOffloaded+=_targetOffloadedTasksPerRank[i];
    }
  }
  str="  total target offloaded tasks: "+std::to_string(totalTargetOffloaded)+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  offloaded tasks:\n";
  logInfo("printCumulativeStatistics", str);
  int totalOffloaded=0;
  for(int i=0;i<nnodes;i++) {
    if(_offloadedTasksPerRank[i]>0){
      str="    to rank: "+std::to_string(i)+" : "+std::to_string(_offloadedTasksPerRank[i])+"\n";
      logInfo("printCumulativeStatistics", str);
    }
    totalOffloaded+=_offloadedTasksPerRank[i];
  }
  str="  total offloaded tasks: "+std::to_string(totalOffloaded)+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  offloadable tasks that failed threshold requirement: "+std::to_string(_thresholdFails)+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  received tasks:\n";
  logInfo("printCumulativeStatistics", str);
  int totalReceived=0;
  for(int i=0;i<nnodes;i++) {
    if(_receivedTasksPerRank[i]>0) {
      str="    from rank: "+std::to_string(i)+" : "+std::to_string(_receivedTasksPerRank[i])+"\n";
      logInfo("printCumulativeStatistics", str);
    }
    totalReceived+=_receivedTasksPerRank[i];
  }
  str="  total received tasks: "+std::to_string(totalReceived)+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total computation time: "+std::to_string(static_cast<double>(_accComputationTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total handling time: "+std::to_string(static_cast<double>(_accHandlingTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total useful communication time: "+std::to_string(static_cast<double>(_accUsefulCommunicationTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total idle communication time: "+std::to_string(static_cast<double>(_accIdleCommunicationTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total progress time: "+std::to_string(static_cast<double>(_accProgressTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total progressRequests time: "+std::to_string(static_cast<double>(_accProgressRequestsTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total poll time: "+std::to_string(static_cast<double>(_accPollTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total wait time for tasks: "+std::to_string(static_cast<double>(_accWaitTasksTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total wait time for workers: "+std::to_string(static_cast<double>(_accWaitForWorkersTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  str="  total wait time for global master: "+std::to_string(static_cast<double>(_accWaitForGlobalMasterTime/1E06))+"\n";
  logInfo("printCumulativeStatistics", str);
  

#endif
}

