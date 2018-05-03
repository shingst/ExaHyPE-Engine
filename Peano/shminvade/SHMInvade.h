// @todo Martin

#ifndef _SHMINVADE_SHMINVADE_H_
#define _SHMINVADE_SHMINVADE_H_


#include <sys/types.h>
#include <set>


// @todo Remove Namespace
namespace shminvade {
  class SHMInvade;
}


/**
 * Fundamental invasion object. If you create it, you have tell the object
 * how many cores you'd like to invade. You may also try to invade all cores.
 * Then, the object asks the strategy to identify whether there are cores
 * available and, if successful, memorizes those.
 *
 * If you call retreat() or if the object is destroyed, it tells the
 * SHMController that all booked cores are not required anymore.
 */
class shminvade::SHMInvade {
  private:
    std::set<pid_t> _occupiedThreads;
  public:
    static constexpr int MaxThreads = -1;

    SHMInvade(int threads);
    ~SHMInvade();
    void retreat();
};

#endif
