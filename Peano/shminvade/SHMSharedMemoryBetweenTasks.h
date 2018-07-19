/**
Copyright (C) 2018, Martin Schreiber and Tobias Weinzierl

All rights reserved.

Open Source License
===================

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

4. Scientific and commercial work using the software should cite the authors' corresponding papers.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _SHMINVADE_SHM_SHARED_MEMORY_BETWEEN_TASKS_H_
#define _SHMINVADE_SHM_SHARED_MEMORY_BETWEEN_TASKS_H_


#include <sys/types.h>
#include <set>
#include <tbb/spin_mutex.h>


// Maximum number of processs
#if !defined(SHM_INVADE_MAX_PROGRAMS)
#define SHM_INVADE_MAX_PROGRAMS 256
#endif



#ifndef SHMINVADE_USER_DATA_SIZE
/**
 * The threads can all exchange data through their shared memory environment.
 * By default, this shared area is four doubles big. You can however also make
 * it bigger by resetting the SHMINVADE_USER_DATA_SIZE at compile time
 * manually.
 */
#define SHMINVADE_USER_DATA_SIZE  (sizeof(double)*4)
#endif



// Maximum number of cores for large NUMA systems (Ultraviolet, e.g.)
#if !defined(SHM_INVADE_MAX_CORES)
#define SHM_INVADE_MAX_CORES 1024
#endif


// Name of shm filename to use for interprocess shared memory area
#if !defined(SHM_INVADE_SHM_FILE_NAME)
#define SHM_INVADE_SHM_FILE_NAME  "/shminvade"
#endif


namespace shminvade {
  class SHMSharedMemoryBetweenTasks;
}


class shminvade::SHMSharedMemoryBetweenTasks {
  private:
    struct SharedData {
      void lock();
      void unlock();

      tbb::atomic<bool> isLocked;

      tbb::atomic<int> noOfRegisteredProcesses;

      /**
       * Array with pids. Unused pids are -1 - the array is fixed size, so we
       * need dummies.
       */
      volatile pid_t    registeredProcessPids[SHM_INVADE_MAX_PROGRAMS];

      tbb::atomic<int>  owningProcessOfCore[SHM_INVADE_MAX_CORES];

      /**
       * User-specific data per process
       */
      char userDataPerProcess[SHM_INVADE_MAX_PROGRAMS][SHMINVADE_USER_DATA_SIZE];
    };

    SharedData* volatile  _globalSHMData;

    SHMSharedMemoryBetweenTasks();
  public:
    static SHMSharedMemoryBetweenTasks& getInstance();

    int getNumberOfRegisteredProcesses() const;

    template < typename T >
    void setSharedUserData(
      const T *i_data,
      int numberOfEntriesInData
    ) {
      setSharedUserData((const char*)i_data,sizeof(T)*numberOfEntriesInData);
    }

    bool tryToBookCoreForProcess(int coreNumber);
    void freeCore(int coreNumber);

//    template <>
    void setSharedUserData(
      const char*   data,
      int           numberOfEntriesInData
    );

    template <typename T>
    T getSharedUserData(
      int i_process_idx,  ///< index of one particular process, this is *NOT* the process ID!!!
      int i     ///< ID of entry
    ) const {
      T* pointer = (T*)( getSharedUserDataPointer(i_process_idx) );
      return pointer[i];
      //((T*)(getSharedUserDataPointer(i_process_idx)[0]))[i];
    }

    const char* getSharedUserDataPointer(
      int i_process_idx  ///< index of one particular process, this is *NOT* the process ID!!!
    ) const;

    void cleanUp();

    /**
     * Returns the entry in the shared table that corresponds to this very
     * process.The routine also is a lazy initialisation. If the process is
     * not yet known in the shared region, then it adds it to the list.
     */
    int getProcessIndexInSharedDataTable( int myId = (pid_t) syscall (__NR_getpid) );

    std::string getCoreProcessAssociation() const;
};

#endif
