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

#ifndef EXAHYPE_EXAHYPE_OFFLOADING_REPLICATIONSTATISTICS_H_
#define EXAHYPE_EXAHYPE_OFFLOADING_REPLICATIONSTATISTICS_H_

#include <atomic>
#include "tarch/logging/Log.h"

namespace exahype {
namespace offloading {

class ReplicationStatistics {

private:

    static tarch::logging::Log _log;

	std::atomic<int> _spawnedTasks;
	std::atomic<int> _executedTasks;
	std::atomic<int> _savedTasks;
	std::atomic<int> _receivedTasks;
	std::atomic<int> _sentTasks;

	ReplicationStatistics();
	virtual ~ReplicationStatistics();
public:

    static ReplicationStatistics& getInstance();

    void printStatistics();

    void notifyReceivedTask();
    void notifySentTask();
    void notifySavedTask();
    void notifySpawnedTask();
    void notifyExecutedTask();
};

} /* namespace offloading */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_OFFLOADING_REPLICATIONSTATISTICS_H_ */
