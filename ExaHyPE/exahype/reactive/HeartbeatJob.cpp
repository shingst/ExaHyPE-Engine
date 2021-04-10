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

#if  defined(SharedTBB)  && defined(Parallel)

#include "exahype/reactive/HeartbeatJob.h"

#include <mpi.h>
#include <unistd.h>
#include "peano/datatraversal/TaskSet.h"

namespace exahype {

namespace reactive {

HeartbeatJob* HeartbeatJob::_singleton;
tarch::logging::Log HeartbeatJob::_log( "exahype::reactive::HeartbeatJob");

HeartbeatJob::HeartbeatJob() :
 tarch::multicore::jobs::Job( tarch::multicore::jobs::JobType::BackgroundTask, 0 , tarch::multicore::DefaultPriority),
 _terminateTrigger(false),
 _hasTerminated(false),
 _lastTimeStampTriggered(0)
{
	// TODO Auto-generated constructor stub

}

HeartbeatJob::~HeartbeatJob() {
	// TODO Auto-generated destructor stub
}

void HeartbeatJob::startHeartbeatJob() {
  _singleton = new HeartbeatJob();
  _singleton->_lastTimeStampTriggered = MPI_Wtime();
  MPI_Sendrecv(MPI_IN_PLACE, 0, MPI_BYTE, MPI_PROC_NULL, 1, MPI_IN_PLACE, 0, MPI_BYTE, MPI_PROC_NULL, 0, MPI_COMM_SELF, MPI_STATUS_IGNORE);

  peano::datatraversal::TaskSet spawned(_singleton);

}

void HeartbeatJob::stopHeartbeatJob() {
  _singleton->_terminateTrigger = true;
  while(!_singleton->_hasTerminated) {
	  usleep(5);
  }
}

bool HeartbeatJob::run(bool runOnMasterThread) {
	if(!_terminateTrigger) {
	  double curTime = MPI_Wtime();

	  if(curTime-_lastTimeStampTriggered>1) {
		  MPI_Sendrecv(MPI_IN_PLACE, 0, MPI_BYTE, MPI_PROC_NULL, -1, MPI_IN_PLACE, 0, MPI_BYTE, MPI_PROC_NULL, 0, MPI_COMM_SELF, MPI_STATUS_IGNORE);
		  logInfo("run()", "triggering new heartbeat")
		  MPI_Sendrecv(MPI_IN_PLACE, 0, MPI_BYTE, MPI_PROC_NULL, 1, MPI_IN_PLACE, 0, MPI_BYTE, MPI_PROC_NULL, 0, MPI_COMM_SELF, MPI_STATUS_IGNORE);
          _lastTimeStampTriggered = curTime;
	  }

	  return true;
	}
	else {
	  _hasTerminated = true;
	  return false;
	}
}

}

}

#endif
