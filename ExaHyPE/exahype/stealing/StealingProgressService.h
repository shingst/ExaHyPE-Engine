#if !defined(_EXAHYPE_STEALINGPROGRESSSERVICE_H_) && defined(SharedTBB)  && defined(Parallel) && defined(DistributedStealing)
#define _EXAHYPE_STEALINGPROGRESSSERVICE_H_

#include "tarch/logging/Log.h"
#include "tarch/services/Service.h"

#include "exahype/solvers/ADERDGSolver.h"

namespace exahype {
  namespace stealing {
    class StealingProgressService;
  }
}


class exahype::stealing::StealingProgressService: public tarch::services::Service {
  private:
    exahype::solvers::ADERDGSolver* _solver;
    bool _isSet;
    static tarch::logging::Log _log;
  public:
    StealingProgressService();
    static StealingProgressService& getInstance();
    virtual ~StealingProgressService();
    virtual void receiveDanglingMessages();
    void setSolver(exahype::solvers::ADERDGSolver *solver);
};

#endif

