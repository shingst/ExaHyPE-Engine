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

#if defined(SharedTBB)  && defined(Parallel) && defined(MemoryMonitoring)

#include "exahype/offloading/MemoryMonitor.h"
#include "tarch/services/ServiceFactory.h"
#include "tarch/multicore/Jobs.h"
#include "peano/utils/UserInterface.h"

registerService(exahype::offloading::MemoryMonitor);

tarch::logging::Log exahype::offloading::MemoryMonitor::_log("exahype::offloading::MemoryMonitor");

exahype::offloading::MemoryMonitor::MemoryMonitor()
 {}

exahype::offloading::MemoryMonitor::~MemoryMonitor() {}

std::size_t exahype::offloading::MemoryMonitor::getFreeMemMB() {

  char   work[256];
  FILE*  f;

  sprintf(work, "/proc/meminfo");
  f = fopen(work, "r");

  if (f == NULL) {
    std::ostringstream msg;
    msg << "can't open file " << work;
    _log.error("getFreeMem()", msg.str() );
    return(0);
  }

  char *line_buf = work;
  size_t size = 256;
  size_t line_size = 0;
  line_size = getline(&line_buf, &size, f);
  line_size = getline(&line_buf, &size, f);
  //logInfo("getFreeMemMB",  line_buf);

  fclose(f);

  size_t result;
  char *p;

  for(int i=0; i<size; i++) {
    if(isdigit(line_buf[i])) {
      result = strtol(&line_buf[i], &line_buf, 10);
      break;
    }
  }

  //logInfo("getFreeMemMB", "result "<<result);


  return result/1024;
}


void exahype::offloading::MemoryMonitor::receiveDanglingMessages() {

  size_t freeMem = getFreeMemMB();

  if(freeMem<1000)
    logInfo("receiveDanglingMessage(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB"
                                          <<" free memory = "<< getFreeMemMB()<< " MB ");

}

exahype::offloading::MemoryMonitor& exahype::offloading::MemoryMonitor::getInstance() {
  static MemoryMonitor service;
  return service;
}

#endif
