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
  global_shm_data(nullptr) {

  int fd = shm_open(SHM_INVADE_SHM_FILE_NAME, (O_CREAT | O_RDWR), S_IRUSR | S_IWUSR);

  if (fd == -1)
    throw(std::string("Cannot create shared memory object ") + std::string(SHM_INVADE_SHM_FILE_NAME));

  if (ftruncate (fd, sizeof(SharedData)) == -1)
    throw(std::string("Cannot resize shared memory object"));

  global_shm_data = (SharedData*) mmap (NULL, sizeof(SharedData), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (global_shm_data == MAP_FAILED)
    throw(std::string("map object into memory"));

  global_shm_data->is_locked = false;
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

  std::memcpy((void*)&((global_shm_data->user_data_per_process[myIndex][0])), (data), numberOfEntriesInData);
}


const char* shminvade::SHMSharedMemoryBetweenTasks::getSharedUserDataPointer(
  int i
) const {
  if (i >= SHM_INVADE_MAX_PROGRAMS)
    throw(std::string("Too many processes! SHM_INVADE_MAX_PROGRAMS too small"));

  #if SHM_INVADE_DEBUG>=8
  std::cout << SHM_DEBUG_PREFIX <<  "computed pointer from shared memory environment for entry " << i << " (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
  #endif

  return global_shm_data->user_data_per_process[i];
}


void shminvade::SHMSharedMemoryBetweenTasks::cleanUp() {
  global_shm_data->lock();

  global_shm_data->num_registered_threads   = 0;
  global_shm_data->num_registered_processes = 0;

  for (int i=0; i<SHM_INVADE_MAX_THREADS; i++) {
    global_shm_data->owning_process_of_thread[i] = -1;
  }

  global_shm_data->unlock();

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "cleaned up shared memory region (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
  #endif

}


void shminvade::SHMSharedMemoryBetweenTasks::SharedData::lock() {
  bool hasBeenLockedByThisThread = false;

  while (!hasBeenLockedByThisThread) {
    hasBeenLockedByThisThread = is_locked.fetch_and_store(true);
  }
}


void shminvade::SHMSharedMemoryBetweenTasks::SharedData::unlock() {
  is_locked.fetch_and_store(false);
}


int shminvade::SHMSharedMemoryBetweenTasks::getProcessIndexInSharedDataTable(int myId) {
  int numberOfProcesses = global_shm_data->num_registered_processes;

  for (int i=0; i<numberOfProcesses; i++) {
    if (global_shm_data->registered_process_pids[i] == myId) return i;
  }

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "process " << myId << " not known yet: insert entry into shared memory region (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
  #endif

  global_shm_data->lock();
  numberOfProcesses = global_shm_data->num_registered_processes;
  global_shm_data->registered_process_pids[numberOfProcesses] = myId;
  global_shm_data->num_registered_processes++;
  global_shm_data->unlock();

  // register master thread
  int threadsLineInTable = getThreadIndexInSharedDataTable();
  global_shm_data->owning_process_of_thread[threadsLineInTable] = myId;

  return numberOfProcesses;
}


bool shminvade::SHMSharedMemoryBetweenTasks::tryToBookThreadForProcess(int threadId) {
  const int indexOfBookedThread = getThreadIndexInSharedDataTable(threadId);

  const int previousOwner = (global_shm_data->owning_process_of_thread[indexOfBookedThread]).compare_and_swap(
    (pid_t) syscall (__NR_getpid), -1
  );

  return previousOwner==-1;
}



void shminvade::SHMSharedMemoryBetweenTasks::freeThread(int threadId) {
  const int indexOfBookedThread = getThreadIndexInSharedDataTable(threadId);

  global_shm_data->owning_process_of_thread[indexOfBookedThread] = -1;
}


int shminvade::SHMSharedMemoryBetweenTasks::getThreadIndexInSharedDataTable(int myId) {
  int numberOfThreads = global_shm_data->num_registered_threads;

  for (int i=0; i<numberOfThreads; i++) {
    if (global_shm_data->registered_thread_pids[i] == myId) return i;
  }

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "thread " << myId << " not known yet: insert entry into shared memory region (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
  #endif

  global_shm_data->lock();
  numberOfThreads = global_shm_data->num_registered_threads;
  global_shm_data->registered_thread_pids[numberOfThreads] = myId;
  global_shm_data->num_registered_threads++;
  global_shm_data->unlock();

  return numberOfThreads;
}


int shminvade::SHMSharedMemoryBetweenTasks::getNumberOfRegisteredProcesses() const {
  return global_shm_data->num_registered_processes;
}


int shminvade::SHMSharedMemoryBetweenTasks::getNumberOfRegisteredThreads() const {
  return global_shm_data->num_registered_threads;
}


std::string shminvade::SHMSharedMemoryBetweenTasks::getThreadProcessAssociation() const {
  std::ostringstream out;

  out << "{";
  for (int i=0; i<getNumberOfRegisteredThreads(); i++) {
    out << "("
        << global_shm_data->registered_thread_pids[i]
        << ":"
        << global_shm_data->owning_process_of_thread[i]
        << ")";
  }

  out << "}";

  return out.str();
}
