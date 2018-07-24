#include "SHMSharedMemoryBetweenTasks.h"
#include "SHMMacros.h"

#include <cstring>
#include <sys/types.h>
#include <signal.h>
#include <tbb/tbb.h>
#include <cassert>
#include <cstring>
#include <thread>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <linux/unistd.h>
#include <sys/syscall.h>
#include <sys/stat.h>
#include <signal.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sstream>


shminvade::SHMSharedMemoryBetweenTasks::SHMSharedMemoryBetweenTasks():
  _globalSHMData(nullptr) {

  int fd = shm_open(SHM_INVADE_SHM_FILE_NAME, (O_CREAT | O_RDWR), S_IRUSR | S_IWUSR);

  if (fd == -1)
    throw(std::string("Cannot create shared memory object ") + std::string(SHM_INVADE_SHM_FILE_NAME));

  if (ftruncate (fd, sizeof(SharedData)) == -1)
    throw(std::string("Cannot resize shared memory object"));

  _globalSHMData = (SharedData*) mmap (NULL, sizeof(SharedData), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (_globalSHMData == MAP_FAILED)
    throw(std::string("map object into memory"));

  _globalSHMData->isLocked = false;
}


shminvade::SHMSharedMemoryBetweenTasks& shminvade::SHMSharedMemoryBetweenTasks::getInstance() {
  static shminvade::SHMSharedMemoryBetweenTasks singleton;
  return singleton;
}


void shminvade::SHMSharedMemoryBetweenTasks::setSharedUserData(
  const char*   data,
  int           numberOfEntriesInData
) {
  int myIndex = getProcessIndexInSharedDataTable();

  std::memcpy((void*)&((_globalSHMData->userDataPerProcess[myIndex][0])), (data), numberOfEntriesInData);
}


const char* shminvade::SHMSharedMemoryBetweenTasks::getSharedUserDataPointer(
  int i
) const {
  if (i >= SHM_INVADE_MAX_PROGRAMS)
    throw(std::string("Too many processes! SHM_INVADE_MAX_PROGRAMS too small"));

  return _globalSHMData->userDataPerProcess[i];
}


void shminvade::SHMSharedMemoryBetweenTasks::cleanUp() {
  _globalSHMData->lock();

  _globalSHMData->noOfRegisteredProcesses = 0;

  for (int i=0; i<SHM_INVADE_MAX_CORES; i++) {
	_globalSHMData->owningProcessOfCore[i] = -1;
  }

  _globalSHMData->unlock();

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "cleaned up shared memory region (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
  #endif
}


void shminvade::SHMSharedMemoryBetweenTasks::SharedData::lock() {
  bool hasBeenLockedByThisThread = false;

  while (!hasBeenLockedByThisThread) {
    hasBeenLockedByThisThread = isLocked.fetch_and_store(true);
  }
}


void shminvade::SHMSharedMemoryBetweenTasks::SharedData::unlock() {
  isLocked.fetch_and_store(false);
}


int shminvade::SHMSharedMemoryBetweenTasks::getProcessIndexInSharedDataTable(int myId) {
  int numberOfProcesses = _globalSHMData->noOfRegisteredProcesses;

  for (int i=0; i<numberOfProcesses; i++) {
    if (_globalSHMData->registeredProcessPids[i] == myId) return i;
  }

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Process " << myId << " not known yet: insert entry into shared memory region (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
  #endif

  _globalSHMData->lock();
  numberOfProcesses = _globalSHMData->noOfRegisteredProcesses;
  _globalSHMData->registeredProcessPids[numberOfProcesses] = myId;
  _globalSHMData->noOfRegisteredProcesses++;

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Association: " << getCoreProcessAssociation() << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  _globalSHMData->unlock();

  return numberOfProcesses;
}


bool shminvade::SHMSharedMemoryBetweenTasks::tryToBookCoreForProcess(int coreNumber) {
  const int previousOwner = (_globalSHMData->owningProcessOfCore[coreNumber]).compare_and_swap(
    (pid_t) syscall (__NR_getpid), -1
  );

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Tried to book core. Association: " << getCoreProcessAssociation() << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  return previousOwner==-1;
}



void shminvade::SHMSharedMemoryBetweenTasks::freeCore(int coreNumber) {
  _globalSHMData->owningProcessOfCore[coreNumber] = -1;
}


bool shminvade::SHMSharedMemoryBetweenTasks::isBooked(int coreNumber) const {
  return _globalSHMData->owningProcessOfCore[coreNumber] != -1;
}


int shminvade::SHMSharedMemoryBetweenTasks::getNumberOfRegisteredProcesses() const {
  return _globalSHMData->noOfRegisteredProcesses;
}


std::string shminvade::SHMSharedMemoryBetweenTasks::getCoreProcessAssociation() const {
  std::ostringstream out;

  out << "{";
  for (int i=0; i<std::thread::hardware_concurrency(); i++) {
    if (i!=0) out << ",";
    out << "("
    	<< "core:" << i << ","
		<< "owning process:" << _globalSHMData->owningProcessOfCore[i]
        << ")";
  }
  out << ",local process:" << (pid_t) syscall (__NR_getpid);
  out << ",no of registered processes:" << getNumberOfRegisteredProcesses();
  out << ",processes:";
  for (int i=0; i<getNumberOfRegisteredProcesses(); i++) {
	if (i!=0) out << ",";
    out << _globalSHMData->registeredProcessPids[i];
  }
  out << "}";

  return out.str();
}

