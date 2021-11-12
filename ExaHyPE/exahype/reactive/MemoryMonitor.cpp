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

#if defined(SharedTBB)  && defined(Parallel)

#include "exahype/reactive/MemoryMonitor.h"
#include "tarch/services/ServiceFactory.h"
#include "tarch/multicore/Jobs.h"
#include "peano/utils/UserInterface.h"

#include "exahype/parser/Parser.h"

#include <fstream>
#include <iostream>

#ifdef TMPI
#include "teaMPI.h"
#endif

#if defined(MemoryMonitoring)
registerService(exahype::reactive::MemoryMonitor);
#endif

#define SAMPLING_INTERVAL 1 //samples every 1 seconds

tarch::logging::Log exahype::reactive::MemoryMonitor::_log("exahype::reactive::MemoryMonitor");

exahype::reactive::MemorySample::MemorySample(unsigned long elapsed, std::size_t freeMem, std::size_t usedMem )
 : _elapsedSeconds(elapsed), _freeMemMB(freeMem), _usedMemMB(usedMem)
{}

std::string exahype::reactive::MemorySample::to_string() {
  std::string result = "";
  result = result +  std::to_string(_elapsedSeconds) + " " + std::to_string(_usedMemMB) + " " + std::to_string(_freeMemMB);
  return result;
}

exahype::reactive::MemoryMonitor::MemoryMonitor()
 : _timestampLastSample(std::chrono::system_clock::now()),
   _timestampStart(std::chrono::system_clock::now())
{}

exahype::reactive::MemoryMonitor::~MemoryMonitor() {}

std::size_t exahype::reactive::MemoryMonitor::getFreeMemMB() {

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

  size_t read1, read2;
  read1 = getline(&line_buf, &size, f);
  read2 = getline(&line_buf, &size, f);

  if(read1<=0 || read2<=0) {
    logError("getFreeMemMB", "Error reading memory consumption info, exiting...")
    fclose(f);
    exit(-1);
  }

  fclose(f);

  size_t result = 0;

  for(unsigned int i=0; i<size; i++) {
    if(isdigit(line_buf[i])) {
      result = strtol(&line_buf[i], &line_buf, 10);
      break;
    }
  }
  return result/1024;
}

void exahype::reactive::MemoryMonitor::dumpMemoryUsage() {
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

  for(auto& sample : _samples) {
    file<<sample.to_string()<<"\n";
  }

}

void exahype::reactive::MemoryMonitor::receiveDanglingMessages() {

  size_t freeMem = getFreeMemMB();
  size_t usedMem = peano::utils::UserInterface::getMemoryUsageMB();

  if(freeMem<1000)
    logInfo("receiveDanglingMessage(...)", "used memory =" << usedMem << " MB"
                                          <<" free memory ="<< freeMem << " MB ");

#if defined(MemoryMonitoringTrack)
  if(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-_timestampLastSample).count()>=SAMPLING_INTERVAL) {
    _timestampLastSample = std::chrono::system_clock::now();
    _samples.push_back(MemorySample(std::chrono::duration_cast<std::chrono::seconds>((_timestampLastSample-_timestampStart)).count(), freeMem, usedMem));
  }
#endif

}

void exahype::reactive::MemoryMonitor::setOutputDir(std::string output_dir) {
  _output_dir = output_dir;
}

exahype::reactive::MemoryMonitor& exahype::reactive::MemoryMonitor::getInstance() {
  static MemoryMonitor service;
  return service;
}

#endif
