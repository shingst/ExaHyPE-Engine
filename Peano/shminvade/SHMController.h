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
#ifndef _SHMINVADE_SHMCONTROLLER_H_
#define _SHMINVADE_SHMCONTROLLER_H_


#include <tbb/atomic.h>
#include <tbb/spin_mutex.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/spin_mutex.h>
#include <tbb/task.h>
// This seems to be an Intel requirement as this feature isnt' released yet officially.
#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>


#include "SHMPinningObserver.h"


namespace shminvade {
  class SHMController;
  class SHMLockTask;
  class SHMStrategy;
  class SHMOccupyAllCoresStrategy;
  class SHMMultipleRanksPerNodeStrategy;
}



/**
 *
 *
 * <h2> Usage </h2>
 *
 * - Init your TBB environment with shminvade::SHMController::getInstance().getMaxAvailableCores()
 *   threads.
 * - Use shminvade::SHMStrategy::setStrategy to set a strategy if you want
 *   another one than let all ranks invade all cores simultaneously.
 * - Initialise shared memory regions through shminvade::SHMSharedMemoryBetweenTasks
 *   if you want ranks to communicate via shared memory.
 */
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

      std::string toString() const;
   };


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

    SHMPinningObserver   _pinningObserver;
    tbb::global_control  _globalThreadCountControl;

    tbb::atomic<bool>    _switchedOn;

    /**
     * I'd prefer to use the thread states directly here. However, I have to
     * work with pointers as each entry contains a mutex and TBB does not
     * allow us to copy mutexes. The table maps the core numbers to the state
     * of the thread pinned to this very core.
     */
    typedef tbb::concurrent_hash_map<int, ThreadState*> ThreadTable;
    ThreadTable  _cores;

    /**
     * Read-only operation mainly required by lock tasks
     */
    ThreadType getThreadType( int core ) const;

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
     * May not be called if thread already is known.
     */
    void registerNewCore(int coreNumber, ThreadType initialstate );

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

    void retreat( int core );

    bool tryToBookCore( int core );
};

#endif
