// @todo Martin

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
#if !defined(SHM_INVADE_MAX_THREADS)
#define SHM_INVADE_MAX_THREADS 1024
#endif


// Name of shm filename to use for interprocess shared memory area
#if !defined(SHM_INVADE_SHM_FILE_NAME)
#define SHM_INVADE_SHM_FILE_NAME  "/shminvade"
#endif




// @todo Remove Namespace
namespace shminvade {
  class SHMSharedMemoryBetweenTasks;
}


class shminvade::SHMSharedMemoryBetweenTasks {
  private:
    struct SharedData {
      void lock();
      void unlock();

      // this value is true (1), if the shm was initialized
      tbb::atomic<bool> is_locked;

      // current number of registered processs
      tbb::atomic<int> num_registered_threads;
      tbb::atomic<int> num_registered_processes;

      // array with pids, unused pids are specified with -1. All the
      // others carry the process number of their owner.
      // (pid_t is of type int)
      volatile pid_t registered_process_pids[SHM_INVADE_MAX_PROGRAMS];
      volatile pid_t registered_thread_pids[SHM_INVADE_MAX_THREADS];

      tbb::atomic<int>  owning_process_of_thread[SHM_INVADE_MAX_THREADS];

      // User-specific data per process
      char user_data_per_process[SHM_INVADE_MAX_PROGRAMS][SHMINVADE_USER_DATA_SIZE];
    };

    SharedData* volatile  global_shm_data;

    SHMSharedMemoryBetweenTasks();

    /**
     * Is not const as it creates entries in the table on-the-fly.
     */
    int getThreadIndexInSharedDataTable(int myId = (pid_t) syscall (__NR_gettid));
  public:
    static SHMSharedMemoryBetweenTasks& getInstance();

    int getNumberOfRegisteredProcesses() const;
    int getNumberOfRegisteredThreads() const;

    template < typename T >
    void setSharedUserData(
      const T *i_data,
      int numberOfEntriesInData
    ) {
      setSharedUserData((const char*)i_data,sizeof(T)*numberOfEntriesInData);
    }

    bool tryToBookThreadForProcess(int threadId);
    void freeThread(int threadId);

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

    int getProcessIndexInSharedDataTable( int myId = (pid_t) syscall (__NR_getpid) );

    std::string getThreadProcessAssociation() const;
};

#endif
