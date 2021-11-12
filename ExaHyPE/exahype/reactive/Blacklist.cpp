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

#if defined(Parallel)

#include "Blacklist.h"

#include "exahype/reactive/ReactiveContext.h"
#include "exahype/reactive/AggressiveHybridDistributor.h"
#include "exahype/reactive/PerformanceMonitor.h"

#include "tarch/parallel/Node.h"

#include <algorithm>

tarch::logging::Log exahype::reactive::Blacklist::_log( "exahype::reactive::Blacklist" );

exahype::reactive::Blacklist::Blacklist(){
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  _localBlacklist = new double[nnodes];
  std::fill(&_localBlacklist[0], &_localBlacklist[nnodes], 0);
}

exahype::reactive::Blacklist::~Blacklist() {
  delete [] _localBlacklist;
}

exahype::reactive::Blacklist& exahype::reactive::Blacklist::getInstance() {
  static Blacklist blacklist;
  return blacklist;
}

void exahype::reactive::Blacklist::triggerEmergencyAndBlacklistRank(int rank) {
  switch(exahype::reactive::ReactiveContext::getInstance().getOffloadingStrategy()){
    case ReactiveContext::OffloadingStrategy::AggressiveHybrid:
      exahype::reactive::AggressiveHybridDistributor::getInstance().handleEmergencyOnRank(rank); break;
    default:
      //do nothing, no emergencies supported for other offloading strategies so far
      break;
  }

  _localBlacklist[rank]++;
  exahype::reactive::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_localBlacklist[rank], rank);
}

void exahype::reactive::Blacklist::recoverBlacklistedRanks() {
  logDebug("recoverBlacklistedRanks()","decrease heat of emergency heat map");
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  for(int i=0; i<nnodes;i++) {
    _localBlacklist[i]*= 0.9;
    if(_localBlacklist[i]>0)
      exahype::reactive::PerformanceMonitor::getInstance().submitBlacklistValueForRank(_localBlacklist[i], i);
  }
}

bool exahype::reactive::Blacklist::isBlacklisted(int rank) const {
  // global blacklist check might be moved out to reactive context
  const double* globalBlacklist = exahype::reactive::PerformanceMonitor::getInstance().getBlacklistGlobalSnapshot();
  return (globalBlacklist[rank]>0.5) || (_localBlacklist[rank]>0.5);
}

void exahype::reactive::Blacklist::printBlacklist() {
  int nnodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
  const double* globalBlacklist = exahype::reactive::PerformanceMonitor::getInstance().getBlacklistGlobalSnapshot();

  for(int i=0; i<nnodes; i++) {
    if(globalBlacklist[i]>0 || _localBlacklist[i]>0)
      logInfo("printBlacklist", "blacklist value for rank "<<i<<":"<<globalBlacklist[i]);
  }
}

#endif
