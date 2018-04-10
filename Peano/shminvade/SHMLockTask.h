// @todo Martin

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
    const pid_t  _pid_t;

    int    _sleepTime;

    void reenqueue();

    void terminate();

  public:
    SHMLockTask( pid_t pid_t, int sleepTime = SHM_MIN_SLEEP );

    /**
     * Please note that we use enqueue for task requeuing as TBB does
     * depth-first task enqueuing for efficiency reasons. This means, if we
     * re-enqueues stuff again, it would most likely be ran on the same thread
     * which is exactly not what we want.
     */
    tbb::task* execute();
};


#endif
