#ifndef __HELPER_PRINTONCE__
#define __HELPER_PRINTONCE__

#include "tarch/logging/Log.h"
#include <iostream>

// variant 1
struct printonlyonce {
	bool printed;
	tarch::logging::Log _log;
	
	printonlyonce() : printed(false), _log("printonlyonce") {}
	printonlyonce(const std::string& className) : printed(false), _log(className) {}
	
	void log(std::stringstream& ostream) {
		if(!printed) {
			printed = true;
			_log.info("", ostream.str());
		}
	}
};

// variant 2, to be used as
//   static printonce msg("FooBar");
// exploting the static property -> only instanced once.

struct printonce {
	bool printed;
	tarch::logging::Log _log;
	
	printonce(const std::string& msg) : _log("printonce") {
		_log.info("printonce()", msg);
		//std::cout << "#PRINTONCE# " << msg << std::endl;
	}
	
	printonce(const std::string& identifier, const std::string& msg) : _log("printonce") {
		_log.info(identifier, msg);
		//std::cout << "#PRINTONCE# " << msg << std::endl;
	}
};

#endif /* __HELPER_PRINTONCE__ */
