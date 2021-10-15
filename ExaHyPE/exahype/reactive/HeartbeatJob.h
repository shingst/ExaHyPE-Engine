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

#if !defined(EXAHYPE_OFFLOADING_HEARTBEATJOB_H_) && defined(SharedTBB)  && defined(Parallel)
#define EXAHYPE_OFFLOADING_HEARTBEATJOB_H_

#include "tarch/multicore/Jobs.h"
#include "tarch/logging/Log.h"
#include <atomic>

namespace exahype {

namespace reactive{

/**
 * A heartbeat job regularly posts teaMPI heartbeats after an interval of 1 s (currently hardcoded) has elapsed.
 * It reschedules itself until it has been stopped.
 * There can only be a single heartbeat job per rank.
 */
class HeartbeatJob : public tarch::multicore::jobs::Job {
  public:
	static tarch::logging::Log _log;

	/**
	 * Starts heartbeat job.
	 */
	static void startHeartbeatJob();

	/**
	 * Stops heartbeat job (busy polls until the task has finished).
	 */
	static void stopHeartbeatJob();
	bool run(bool runOnMasterThread);

  private:
	static HeartbeatJob* _singleton;

	/**
	 * Flag that indicates whether the heartbeat job should terminate.
	 */
	std::atomic <bool> _hasSetTerminateTrigger;

	/**
	 * Flag that indicates whether the heartbeat job has terminated.
	 */
	std::atomic <bool> _hasTerminated;

	/**
	 * Stores time stamp of last heartbeat.
	 */
	double _timestampOfLastHeartbeat;

	HeartbeatJob();
	virtual ~HeartbeatJob();
};

}

}
#endif /* EXAHYPE_OFFLOADING_HEARTBEATJOB_H_ */
