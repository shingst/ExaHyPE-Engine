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

#include "exahype/parser/Parser.h"

#include <fstream>
#include <iostream>

#ifdef TMPI
#include "teaMPI.h"
#endif

registerService(exahype::offloading::MemoryMonitor);

tarch::logging::Log exahype::offloading::MemoryMonitor::_log("exahype::offloading::MemoryMonitor");

exahype::offloading::MemoryMeasurement::MemoryMeasurement(unsigned long elapsed, std::size_t freeMem, std::size_t usedMem )
 : _elapsed(elapsed), _freeMem(freeMem), _usedMem(usedMem)
{}

std::string exahype::offloading::MemoryMeasurement::to_string() {
  std::string result = "";
  result = result +  std::to_string(_elapsed) + " " + std::to_string(_usedMem) + " " + std::to_string(_freeMem);
  return result;
}

exahype::offloading::MemoryMonitor::MemoryMonitor()
 : _start(std::chrono::system_clock::now()), _lastMeasurementTimestamp(std::chrono::system_clock::now())
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

  fclose(f);

  size_t result;
  char *p;

  for(int i=0; i<size; i++) {
    if(isdigit(line_buf[i])) {
      result = strtol(&line_buf[i], &line_buf, 10);
      break;
    }
  }

  return result/1024;
}

void exahype::offloading::MemoryMonitor::dumpMemoryUsage() {
  int rank = tarch::parallel::Node::getInstance().getRank();

  int team = 0;
#ifdef TMPI
  team = TMPI_GetTeamNumber();
#endif

  std::string outputFilename = _output_dir + "/memory_stats_r"+std::to_string(rank)+"_t"+std::to_string(team)+".csv";

  std::ofstream file;
  file.open(outputFilename);

  if (!file) {
    logInfo("dumpMemoryUsage()","Failed to open file for writing " << outputFilename);
  }

  file<<"timestamp usedMem freeMem\n";

  for(auto& measurement : _measurements) {
    file<<measurement.to_string()<<"\n";
  }

}

void exahype::offloading::MemoryMonitor::receiveDanglingMessages() {

  size_t freeMem = getFreeMemMB();
  size_t usedMem = peano::utils::UserInterface::getMemoryUsageMB();

  if(freeMem<1000)
    logInfo("receiveDanglingMessage(...)", "used memory =" << usedMem << " MB"
                                          <<" free memory ="<< freeMem << " MB ");

#if defined(MemoryMonitoringTrack)
  if(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-_lastMeasurementTimestamp).count()>=1) {
    _lastMeasurementTimestamp = std::chrono::system_clock::now();
    _measurements.push_back(MemoryMeasurement(std::chrono::duration_cast<std::chrono::seconds>((_lastMeasurementTimestamp-_start)).count(), freeMem, usedMem));
  }
#endif

}

void exahype::offloading::MemoryMonitor::setOutputDir(std::string output_dir) {
  _output_dir = output_dir;
}

exahype::offloading::MemoryMonitor& exahype::offloading::MemoryMonitor::getInstance() {
  static MemoryMonitor service;
  return service;
}

#endif
