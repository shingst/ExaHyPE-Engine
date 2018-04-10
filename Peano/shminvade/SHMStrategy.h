// @todo Martin

#ifndef _SHMINVADE_SHMSTRATEGY_H_
#define _SHMINVADE_SHMSTRATEGY_H_


#include <sys/types.h>
#include <set>


// @todo Remove Namespace
namespace shminvade {
  class SHMStrategy;
}


class shminvade::SHMStrategy {
  private:
    static SHMStrategy*  _activeStrategy;
  public:
    virtual ~SHMStrategy();

    static SHMStrategy& getInstance();

    /**
     * Ownership is transferred to strategy, i.e. you don't have to destroy
     * this instance.
     */
    static void setStrategy(SHMStrategy* strategy);

    virtual std::set<pid_t> invadeThreads(int wantedNumberOfThreads) = 0;

    /**
     * Whatever you do in the retreat, ensure that you forward the call to the
     * SHMController afterwards.
     */
    virtual void retreat(const std::set<pid_t>& threadIds) = 0;

    virtual void cleanUp() = 0;
};


#endif
