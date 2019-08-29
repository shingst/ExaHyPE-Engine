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

#include "ReplicationStatistics.h"
#include "exahype/offloading/OffloadingManager.h"

namespace exahype {
namespace offloading {


tarch::logging::Log ReplicationStatistics::_log( "exahype::offloading::ReplicationStatistics" );


ReplicationStatistics::ReplicationStatistics() :
  _spawnedTasks(0),
  _executedTasks(0),
  _savedTasks(0),
  _receivedTasks(0),
  _sentTasks(0)
{
	// TODO Auto-generated constructor stub

}

ReplicationStatistics::~ReplicationStatistics() {
	// TODO Auto-generated destructor stub
}

ReplicationStatistics& ReplicationStatistics::getInstance() {
	static ReplicationStatistics replicationStats;
	return replicationStats;
}

void ReplicationStatistics::notifyReceivedTask(){
	_receivedTasks++;
}

void ReplicationStatistics::notifySentTask(){
	_sentTasks++;
}

void ReplicationStatistics::notifySavedTask(){
	_savedTasks++;
}

void ReplicationStatistics::notifySpawnedTask(){
    _spawnedTasks++;
}

void ReplicationStatistics::notifyExecutedTask(){
	_executedTasks++;
}


void ReplicationStatistics::printStatistics() {
	int team = exahype::offloading::OffloadingManager::getInstance().getTMPIInterTeamRank();
	logInfo("printStatistics", " team "<<team
			                   <<" spawned tasks = "<<_spawnedTasks
							   <<" executed tasks = "<<_executedTasks
							   <<" saved tasks =  "<<_savedTasks
							   <<" sent tasks = "<<_sentTasks
							   <<" received tasks = "<<_receivedTasks);
}

} /* namespace offloading */
} /* namespace exahype */
