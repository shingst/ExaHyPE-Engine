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
#include <tbb/spin_mutex.h>


#include <map>


namespace shminvade {
  class SHMPinningObserver;
}


/**
 * This implementations follows very closely the PinningObserver of the Peano
 * framework. We originally intended to follow
 *
 * https://software.intel.com/en-us/blogs/2013/10/31/applying-intel-threading-building-blocks-observers-for-thread-affinity-on-intel
 *
 * but this code snippet seems to be buggy as it does not take into account
 * that TBB might terminate and re-issue threads throughout the lifetime.
 *
 * @author Tobias Weinzierl
 */
class shminvade::SHMPinningObserver: public tbb::task_scheduler_observer {
  private:
    static constexpr int MaxCores = sizeof(long int)*8;

    /**
     * Masks are bitfields which hold an entry for each core (hardware thread)
     * the present application is allowed to run on. If you run multiple MPI
     * ranks for example, this is a subset of the actual cores available on a
     * node. Pinning is realised by disabling all bits but one for a particular
     * thread. Obviously, the pin mask has to be a subset of the process mask.
     */
    typedef std::bitset<MaxCores> CPUSetBitfield;

    /**
     * We may not initialise this field before we
     * Is initialised in the constructor PinningObserver().
     */
    CPUSetBitfield           _availableCores;

    tbb::spin_mutex          _mutex;

    int                      _numberOfCPUs;
  public:
    /**
     * Initialise the field _availableCores.
     */
    SHMPinningObserver();
    virtual ~SHMPinningObserver();

    /**
      * Search for a free core in the bitset.
     */
    void on_scheduler_entry( bool ) override;
    void on_scheduler_exit( bool ) override;

    /**
     * Is an override though TBB developers didn't want it to be overriden for
     * whatever reason
     */
    void observe(bool toggle);
};

#endif

