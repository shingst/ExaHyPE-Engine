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
#ifndef _SHMINVADE_SHMLOCKTASK_H_
#define _SHMINVADE_SHMLOCKTASK_H_


#include "SHMMacros.h"

#include <tbb/tbb.h>


namespace shminvade {
  class SHMLockTask;
}


/**
 * A lock task is an interesting beast: We launch it
 * once at startup and we do assign it a thread id. A lock task
 * terminates if when it runs into its own thread and this thread
 * is shutdown or actively used by the process.
 */
class shminvade::SHMLockTask: public tbb::task {
  private:
    const int _core;

    int    _sleepTime;

    void reenqueue();

    void terminate();

  public:
    SHMLockTask( int core, int sleepTime = SHM_MIN_SLEEP );

    /**
     * Please note that we use enqueue for task requeuing as TBB does
     * depth-first task enqueuing for efficiency reasons. This means, if we
     * re-enqueues stuff again, it would most likely be ran on the same thread
     * which is exactly not what we want.
     */
    tbb::task* execute();
};


#endif
