// @todo Martin

#ifndef _SHMINVADE_SHMCONTROLLER_H_
#define _SHMINVADE_SHMCONTROLLER_H_


#include <tbb/atomic.h>
#include <tbb/spin_mutex.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/spin_mutex.h>
#include <tbb/task.h>



// @todo Remove Namespace
namespace shminvade {
  class SHMController;
  class SHMLockTask;
  class SHMStrategy;
  class SHMOccupyAllCoresStrategy;
  class SHMMultipleRanksPerNodeStrategy;
}



class shminvade::SHMController {
  public:
    enum class ThreadType {
      Master,
      ExclusivelyOwned,
      NotOwned,
      Shutdown
    };

    struct ThreadState {
      typedef tbb::spin_mutex  Mutex;

      ThreadState(ThreadType type_):
        type(type_),
        numberOfExistingLockTasks(0) {}

      Mutex      mutex;

      ThreadType type;

      /**
       * It is really important to ensure that not too many lock threads are
       * flying around in the system. So we count them per thread and do issue
       * lock threads if and only if there is a need to do so. Could obviously
       * be a map over bools, but I just prefer to count up and down today.
       */
      int        numberOfExistingLockTasks;
    };

    static std::string toString( ThreadState state );

  private:
    /**
     * TBB otherwise might destroy the context once it thinks that all tasks
     * have terminated. See my own post to Intel at
     *
     * https://software.intel.com/en-us/forums/intel-threading-building-blocks/topic/700057
     *
     * So we create a special task group context for SHMInvade which notably
     * allows us to enqueue all lock tasks into this one.
     */
    static tbb::task_group_context  InvasiveTaskGroupContext;

    tbb::atomic<bool> _switchedOn;

    /**
     * I'd prefer to use the thread states directly here. However, I have to
     * work with pointers as each entry contains a mutex and TBB does not
     * allow us to copy mutexes.
     */
    typedef tbb::concurrent_hash_map<pid_t, ThreadState*> ThreadTable;
    ThreadTable  _threads;

    /**
     * Read-only operation mainly required by lock tasks
     */
    ThreadType getThreadType( pid_t pid ) const;

    /**
     * Just register the master thread through registerNewThread(). This
     * implies that a lock task is launched for the master thread which
     * eventually will go down - the latest when we shut down the system.
     * Before it does so, it will identify lots of other threads that would
     * in theory be available.
     */
    SHMController();

    /**
     * Insert a new thread entry and then launch the lock task for it
     * immediately.
     *
     * May not be called if thread already is known
     */
    void registerNewThread(int threadId, ThreadType initialstate=ThreadType::NotOwned );

    /**
     * Is typically used by the strategies when they go down.
     */
    void retreatFromAllCores();

    friend class SHMLockTask;
    friend class SHMStrategy;
    friend class SHMOccupyAllCoresStrategy;
    friend class SHMMultipleRanksPerNodeStrategy;
  public:
    ~SHMController();

    static SHMController&  getInstance();

    void switchOn();
    void switchOff();

    int getMaxAvailableCores() const;
    int getFreeCores() const;
    int getBookedCores() const;

    /**
     * The shutdown sets all internal thread states to Shutdown such that all
     * lock threads hanging around in the system shut down, too. It then goes
     * to sleep to allow the lock threads to see this new state being a signal
     * to them. Before, we switch the invasion overall off and thus can then
     * terminate.
     */
    void shutdown();


    void retreatFromThread( pid_t pid );

    bool tryToBookThread( pid_t pid );
};

#endif
