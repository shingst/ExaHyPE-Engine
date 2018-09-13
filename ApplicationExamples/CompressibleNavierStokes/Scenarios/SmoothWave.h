//
// Created by lukas on 27/08/18.
//

#ifndef NAVIERSTOKES_SMOOTHWAVE_H
#define NAVIERSTOKES_SMOOTHWAVE_H

#include "Scenario.h"

namespace NavierStokes {
class SmoothWave : public Scenario {
  void initialValues(const double *const x, const PDE &ns,
                     Variables &vars) final override;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_SMOOTHWAVE_H
