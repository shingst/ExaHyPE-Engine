// @todo Martin

#ifndef _SHMINVADE_SHM_MULTIPLE_RANKS_PER_NODE_STRATEGY_H_
#define _SHMINVADE_SHM_OCCUPY_ALL_CORES_STRATEGY_H_


#include "SHMStrategy.h"


// @todo Remove Namespace
namespace shminvade {
  class SHMMultipleRanksPerNodeStrategy;
}


class shminvade::SHMMultipleRanksPerNodeStrategy: public shminvade::SHMStrategy {
  public:
    SHMMultipleRanksPerNodeStrategy();
    virtual ~SHMMultipleRanksPerNodeStrategy();

    std::set<pid_t> invadeThreads(int wantedNumberOfThreads) override;
    void cleanUp() override;

    void retreat(const std::set<pid_t>& threadIds) override;
};

#endif
