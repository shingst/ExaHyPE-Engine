// @todo Martin

#ifndef _SHMINVADE_SHM_OCCUPY_ALL_CORES_STRATEGY_H_
#define _SHMINVADE_SHM_OCCUPY_ALL_CORES_STRATEGY_H_


#include "SHMStrategy.h"


// @todo Remove Namespace
namespace shminvade {
  class SHMOccupyAllCoresStrategy;
}


class shminvade::SHMOccupyAllCoresStrategy: public shminvade::SHMStrategy {
  public:
    SHMOccupyAllCoresStrategy();
    virtual ~SHMOccupyAllCoresStrategy();

    std::set<pid_t> invadeThreads(int wantedNumberOfThreads) override;
    void cleanUp() override;

    void retreat(const std::set<pid_t>& threadIds) override;
};

#endif
