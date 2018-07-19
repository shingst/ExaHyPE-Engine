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
#ifndef _SHMINVADE_SHMPINNINGOBSERVER_H_
#define _SHMINVADE_SHMPINNINGOBSERVER_H_


#include <sched.h>
#include <stdlib.h>
#include <unistd.h>
#include <bitset>


#include <tbb/task_arena.h>
#include <tbb/task_scheduler_observer.h>
#include <tbb/atomic.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/tick_count.h>


#include <map>


namespace shminvade {
  class SHMPinningObserver;
}


/**
 * This implementations follows very closely
 *
 * https://software.intel.com/en-us/blogs/2013/10/31/applying-intel-threading-building-blocks-observers-for-thread-affinity-on-intel
 *
 * It is basically also the same as we use it in Peano's tarch though we make
 * it slightly more specific, i.e. invasion-focused.
 *
 * The file is an observer which is called every time TBB launches a thread.
 * The observer routine then pins this very thread to the next free core. As it
 * is called for every new thread, this observer always knows how many threads
 * do physically exist.
 *
 * Furthermore, it can maintain a map that tells us exactly which thread id
 * from Linux is mapped/pinned onto which core.
 *
 * @author Leonhard Rannabauer
 * @author Tobias Weinzierl
 */
class shminvade::SHMPinningObserver: public tbb::task_scheduler_observer {
  private:
    /**
     * Masking being available to process. This is basically a bitfield which
     * holds an entry for each core (hardware thread) the present application
     * is allowed to run on. If you run multiple MPI ranks for example, this
     * is a subset of the actual cores available on a node. The fild is
     * initialised in the constructor.
     */
    cpu_set_t*    _mask;
    //AffinityMask*    _mask;

    int              _ncpus;

    /**
     * How many threads have been registered through callback
     */
    tbb::atomic<int> _numThreads;

    std::map<pid_t,int> _threadIdToCoreMap;

    /**
     * If the observer is switched on, it automatically pins all threads.
     */
    void pinCurrentThread();

    typedef std::bitset<sizeof(long int)*8> CPUSetBitfield;
    static CPUSetBitfield cpuSetMaskToBitfield( cpu_set_t   _mask);
//    static std::string toString( CPUSetBitfield bitfield );
  public:
    SHMPinningObserver();
    virtual ~SHMPinningObserver();

    void on_scheduler_entry( bool ) override;
    void on_scheduler_exit( bool ) override;

    int getNumberOfRegisteredThreads() const;

    int getCoreOfThread(pid_t threadId) const;
};

#endif

