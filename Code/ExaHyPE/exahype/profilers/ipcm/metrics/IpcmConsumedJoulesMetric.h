#ifndef _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CONSUMED_JOULES_METRIC_H_
#define _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CONSUMED_JOULES_METRIC_H_

#ifdef IPCM_AVAILABLE

#include "GenericIpcmMetric.h"

namespace exahype {
namespace profilers {
namespace ipcm {

struct __IpcmConsumedJoulesMetric {
  using method_return_t = double;

  static method_return_t method(const SystemCounterState& before,
                                const SystemCounterState& after) {
    return getConsumedJoules(before, after);
  }

  static const char* method_tag() { return "ConsumedEnergy_Joule"; }
};

using IpcmConsumedJoulesMetric = GenericIpcmMetric<__IpcmConsumedJoulesMetric>;

}  // namespace ipcm
}  // namespace profilers
}  // namespace exahype

#endif  // IPCM_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_IPCM_METRICS_IPCM_CONSUMED_JOULES_METRIC_H_
