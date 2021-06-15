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

#include "TimeStampAndTriggerTeamHistory.h"

#include "exahype/reactive/OffloadingContext.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"

#include "tarch/logging/Log.h"

#include <cassert>
#include <string>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <limits>

tarch::logging::Log exahype::reactive::TimeStampAndTriggerTeamHistory::_log("exahype::reactive::TimeStampAndLimiterHistory");


exahype::reactive::TimeStampAndTriggerTeamHistory::TimeStampAndTriggerTeamHistory() :
  _lastConsistentTimeStepPtr(-1)
{
  _timestamps = new std::vector<double> [exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()];
  _timestepSizes = new std::vector<double> [exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()];
  _triggerStatuses = new std::vector<bool> [exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams()];
}

exahype::reactive::TimeStampAndTriggerTeamHistory::~TimeStampAndTriggerTeamHistory() {
  delete[] _triggerStatuses;
  delete[] _timestepSizes;
  delete[] _timestamps;
}

exahype::reactive::TimeStampAndTriggerTeamHistory& exahype::reactive::TimeStampAndTriggerTeamHistory::getInstance() {
  static TimeStampAndTriggerTeamHistory history;
  return history;
}

void exahype::reactive::TimeStampAndTriggerTeamHistory::forwardLastConsistentTimeStepPtr() {
  tarch::multicore::Lock lock(_semaphore, true);
  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  size_t maxIdx = std::numeric_limits<int>::max();

  for(size_t i=0; i<exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams(); i++) {
    maxIdx = std::min(_timestamps[i].size(), maxIdx);
  }

  int i = std::max(0, _lastConsistentTimeStepPtr);
  while(i<maxIdx && _timestamps[otherTeam][i]==_timestamps[myTeam][i]
                 && _triggerStatuses[otherTeam][i]== _triggerStatuses[myTeam][i]  //we may only want to consider a time step consistent if limiter status is both zero (for linear applications)
                 && _timestepSizes[otherTeam][i]==_timestepSizes[myTeam][i]) {
    i++;
  }
  _lastConsistentTimeStepPtr = --i;
  assert(_lastConsistentTimeStepPtr>=-1);
  lock.free();
}

/*void exahype::reactive::TimeStampAndLimiterTeamHistory::trackEstimatedTimeStepSizeLocally(double timeStamp, double estTimeStepSize){
  tarch::multicore::Lock lock(_semaphore, true);
  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPIInterTeamRank();

  int size = _timestamps[myTeam].size();
  if(size>0 && _timestamps[myTeam][size-1]==timeStamp) {
    _estimatedTimestepSizes.push_back(estTimeStepSize);
  }
  else if (size>0 && _timestamps[myTeam][size-1]>timeStamp) {
    int idx = size-2;
    while(_timestamps[myTeam][idx]!=timeStamp) {
      idx--;
    }
    assert(idx>=0);
    _estimatedTimestepSizes[idx]=estTimeStepSize;
  }
  else if (size==0 || _timestamps[myTeam][size-1]<timeStamp){
    _timestamps[myTeam].push_back(timeStamp);
    _estimatedTimestepSizes.push_back(estTimeStepSize);
  }

  lock.free();
}*/

void exahype::reactive::TimeStampAndTriggerTeamHistory::trackTimeStepAndTriggerActive(int team, double timeStamp, double timeStepSize, bool triggerActive, double estimatedTimeStepSize) {
  tarch::multicore::Lock lock(_semaphore, true);
  logDebug("trackTimeStepAndLimiterActive", "Tracking time stamp="<<timeStamp<<" for team "<<team);

  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int size = _timestamps[team].size();
  if(size>0 && _timestamps[team][size-1]==timeStamp) {
    /**if(_timestepSizes[team][size-1]==timeStepSize) {
      //only update limiter
      _triggerStatuses[team][size-1] = _triggerStatuses[team][size-1] || limiterActive;

      if(team==myTeam) {
        _estimatedTimestepSizes[size-1]= std::min(_estimatedTimestepSizes[size-1], estimatedTimeStepSize);
      }
    }
    else {
      _timestamps[team].push_back(timeStamp);
      _timestepSizes[team].push_back(timeStepSize);
      _triggerStatuses[team].push_back(limiterActive);
      if(team==myTeam) {
        _estimatedTimestepSizes.push_back(estimatedTimeStepSize);
      }
    }*/
    ////only update limiter
    _triggerStatuses[team][size-1] = _triggerStatuses[team][size-1] || triggerActive;
    //assert(_timestepSizes[team][size-1]==std::min(_timestepSizes[team][size-1], timeStepSize));
    _timestepSizes[team][size-1] = std::min(_timestepSizes[team][size-1], timeStepSize);
    if(team==myTeam) {
      _estimatedTimestepSizes[size-1]= std::min(_estimatedTimestepSizes[size-1], estimatedTimeStepSize);
    }
  }
  else if (size>0 && _timestamps[team][size-1]>timeStamp) {
    int idx = size-2;
    while(idx>=0 && _timestamps[team][idx]!=timeStamp) {
      idx--;
    }
    if(idx<0) {
      logError("trackTimeStepAndLimiterActive", "Have received timestamp from another team which seems to be old but has not been observed. This is expected if the other team has done a rollback. Appending new time stamp..");
      _triggerStatuses[team].push_back(triggerActive);
      _timestamps[team].push_back(timeStamp);
      _timestepSizes[team].push_back(timeStepSize);
    }
    else {
      _triggerStatuses[team][idx] = _triggerStatuses[team][idx] || triggerActive;
      _timestepSizes[team][idx] = std::min(_timestepSizes[team][idx], timeStepSize);
      if(team==myTeam) {
        _estimatedTimestepSizes[idx]= std::min(_estimatedTimestepSizes[idx], estimatedTimeStepSize);
      }
    }
  }
  else if (size==0 || _timestamps[team][size-1]<timeStamp){
    _timestamps[team].push_back(timeStamp);
    _timestepSizes[team].push_back(timeStepSize);
    _triggerStatuses[team].push_back(triggerActive);
    if(team==myTeam) {
      _estimatedTimestepSizes.push_back(estimatedTimeStepSize);
    }
  }
  //printHistory();
  lock.free();
}

//todo: doesn't scale well with increasing numbers of timesteps -> can check more efficiently in principle
bool exahype::reactive::TimeStampAndTriggerTeamHistory::checkConsistency() {

  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  tarch::multicore::Lock lock(_semaphore, true);
  size_t maxIdx = std::numeric_limits<int>::max();

  for(int i=0; i<exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams(); i++) {
    maxIdx = std::min(_timestamps[i].size(), maxIdx);
  }

  logDebug("checkConsistency","checking arrays from " <<_lastConsistentTimeStepPtr<<" until "<<maxIdx);

  bool consistentTimeStamps = true;
  bool consistentTimeStepSizes = true;
  bool consistentLimiterStatuses = true;

  for(size_t i=_lastConsistentTimeStepPtr; i<maxIdx; i++) {
    std::vector<double> tmp_timestamps;
    std::vector<double> tmp_timestepsizes;
    std::vector<bool> tmp_statuses;
    for(int t=0; t<exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams(); t++) {
      tmp_timestamps.push_back(_timestamps[t][i]);
      tmp_statuses.push_back(_triggerStatuses[t][i]);
      tmp_timestepsizes.push_back(_timestepSizes[t][i]);
    }
    consistentTimeStamps = consistentTimeStamps && std::all_of(tmp_timestamps.begin(), tmp_timestamps.end(), [tmp_timestamps](double x){ return x==tmp_timestamps[0]; });
    consistentLimiterStatuses = consistentLimiterStatuses && std::all_of(tmp_statuses.begin(), tmp_statuses.end(), [tmp_statuses](double x){ return x==tmp_statuses[0]; });
    consistentTimeStepSizes = consistentTimeStepSizes && std::all_of(tmp_timestepsizes.begin(), tmp_timestepsizes.end(), [tmp_timestepsizes](double x){ return x==tmp_timestepsizes[0]; });

    /*if(!consistentTimeStamps) {
        logDebug("checkConsistency", std::setprecision(30)<<"time stamp "<<tmp_timestamps[0]<<", "<<tmp_timestamps[1]<<" equal = "<<(tmp_timestamps[0]==tmp_timestamps[1]));
    }

    if(!consistentTimeStepSizes) {
        logDebug("checkConsistency", "i = "<<i<<std::setprecision(30)<<" : time step sizes "<<tmp_timestepsizes[0]<<", "<<tmp_timestepsizes[1]<<" equal = "<<(tmp_timestepsizes[0]==tmp_timestepsizes[1]));
    }*/
  }

  if((!consistentTimeStamps || !consistentLimiterStatuses || !consistentTimeStepSizes)) {
    //logError("checkConsistency"," Time stamps or limiter statuses are diverged between teams! Consistent stamps = "<<consistentTimeStamps
    //    <<" consistent limiter statuses = "<<consistentLimiterStatuses<<" consistent time step sizes="<<consistentTimeStepSizes);
    //logError("checkConsistency","team="<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<": Time stamps or trigger statuses are inconsistent between teams! There must have been a soft error on at least one team.");
    if((long int)_triggerStatuses[myTeam].size()> _lastConsistentTimeStepPtr && (long int) _triggerStatuses[otherTeam].size()>_lastConsistentTimeStepPtr) {
      if(_triggerStatuses[myTeam][_lastConsistentTimeStepPtr+1]==1
         &&_triggerStatuses[otherTeam][_lastConsistentTimeStepPtr+1]==0) {
        logError("checkConsistency","team="<<myTeam<<" should be the faulty one, as the limiter was activated there.");
      }
      else if(_triggerStatuses[otherTeam][_lastConsistentTimeStepPtr+1]==1
          &&_triggerStatuses[myTeam][_lastConsistentTimeStepPtr+1]==0) {
         logError("checkConsistency","team="<<otherTeam<<" should be the faulty one, as the limiter was activated there.");
      }
    }

    //printHistory();
  }
  lock.free();

  return consistentTimeStamps && consistentLimiterStatuses && consistentTimeStepSizes;
}

bool exahype::reactive::TimeStampAndTriggerTeamHistory::otherTeamHasTimeStepData(double timeStamp, double timeStep) {
  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(size_t i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]==timeStamp && _timestepSizes[otherTeam][i]==timeStep)
      return true;
  }
  return false;
}

bool exahype::reactive::TimeStampAndTriggerTeamHistory::otherTeamHasTimeStamp(double timeStamp) {
  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(size_t i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]==timeStamp)
      return true;
  }
  return false;
}

bool exahype::reactive::TimeStampAndTriggerTeamHistory::otherTeamHasLargerTimeStamp(double timeStamp) {
  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  forwardLastConsistentTimeStepPtr();

  int startIndex = std::max(0, _lastConsistentTimeStepPtr);

  for(int i=startIndex; i<_timestamps[otherTeam].size(); i++) {
    if(_timestamps[otherTeam][i]>timeStamp)
      return true;
  }

  return false;
}

void exahype::reactive::TimeStampAndTriggerTeamHistory::getLastConsistentTimeStepData(double& timestamp, double& timestepSize, double& estimated) {
  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

//  size_t maxIdx = std::numeric_limits<int>::max();
//
//  for(int i=0; i<exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams(); i++) {
//    maxIdx = std::min(_timestamps[i].size(), maxIdx);
//  }
//
//  int i = 0;
//  while(i<maxIdx && _timestamps[otherTeam][i]==_timestamps[myTeam][i]
//                 && _limiterStatuses[otherTeam][i]==0
//                 && _limiterStatuses[myTeam][i]==0) {
//    i++;
//  }
//  i--;
//  if(i>0) {
//    timestamp = _timestamps[otherTeam][i];
//    timestepSize = _timestepSizes[otherTeam][i];
//  }
//  else {
//    timestamp = 0; // 0 is always assumed to be consistent
//    timestepSize = -1;
//  }
  forwardLastConsistentTimeStepPtr();

  tarch::multicore::Lock lock(_semaphore, true);
  if(_lastConsistentTimeStepPtr>=0) {
    timestamp = _timestamps[myTeam][_lastConsistentTimeStepPtr];
    timestepSize = _timestepSizes[myTeam][_lastConsistentTimeStepPtr];
    estimated = _estimatedTimestepSizes[_lastConsistentTimeStepPtr];
  }
  else {
    timestamp = 0;
    timestepSize = -1;
    estimated = -1;
  }
  lock.free();
}

void exahype::reactive::TimeStampAndTriggerTeamHistory::resetMyTeamHistoryToLastConsistentTimeStep() {
  tarch::multicore::Lock lock(_semaphore, true);

  int myTeam = exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber();
  int otherTeam = (myTeam + 1) % 2; //todo: support if more than 2 teams are used

  if(_lastConsistentTimeStepPtr>=0 && _lastConsistentTimeStepPtr<_timestamps[myTeam].size()) {
    _timestamps[myTeam].erase(_timestamps[myTeam].begin()+_lastConsistentTimeStepPtr+1, _timestamps[myTeam].end());
    _timestepSizes[myTeam].erase(_timestepSizes[myTeam].begin()+_lastConsistentTimeStepPtr+1, _timestepSizes[myTeam].end());
    _estimatedTimestepSizes.erase(_estimatedTimestepSizes.begin()+_lastConsistentTimeStepPtr+1, _estimatedTimestepSizes.end());
    _triggerStatuses[myTeam].erase(_triggerStatuses[myTeam].begin()+_lastConsistentTimeStepPtr+1, _triggerStatuses[myTeam].end());
  }

  logInfo("resetMyTeamToLastConsistentTimeStep", "reset (local) team history to last consistent time step. Printing history next...");
  printHistory();
  lock.free();
}


void exahype::reactive::TimeStampAndTriggerTeamHistory::printHistory() const {
  for(unsigned int i=0; i<exahype::reactive::OffloadingContext::getInstance().getTMPINumTeams(); i++) {
    /*if(tarch::parallel::Node::getInstance().getNumberOfNodes()>0
      && tarch::parallel::Node::getInstance().isGlobalMaster()) {
      logWarning("printHistory", "Caution: when using more than one rank, history on rank 0 only contains local information!");
    }*/

    std::ostringstream stream;
    stream <<" Timestamps: "<<std::setprecision(30);
    for(unsigned int j=0; j<_timestamps[i].size(); j++) {
      stream << _timestamps[i][j] << " , ";
    }
    logInfo("printHistory","TimeStampAndLimiterHistory for team "<<i<<" on team "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<":"<<stream.str());

    std::ostringstream stream2;
    stream2 << " TimeStepSizes: "<< std::setprecision(30);
    for(unsigned int j=0; j<_timestepSizes[i].size(); j++) {
      stream2 << _timestepSizes[i][j] << " , ";
    }
    logInfo("printHistory", "TimeStampAndLimiterHistory for team "<<i<<" on team "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<":"<<stream2.str());

    std::ostringstream stream3;
    stream3 << " Limiter statuses ";
    for(unsigned int j=0; j<_triggerStatuses[i].size(); j++) {
      stream3 << std::to_string(_triggerStatuses[i][j]) + " , ";
    }
    logInfo("printHistory","TimeStampAndLimiterHistory for team "<<i<<" on team "<<exahype::reactive::OffloadingContext::getInstance().getTMPITeamNumber()<<":"<<stream3.str());
  }
}
