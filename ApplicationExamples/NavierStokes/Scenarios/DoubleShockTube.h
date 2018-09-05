#ifndef NAVIERSTOKES_DOUBLESHOCKTUBE_H
#define NAVIERSTOKES_DOUBLESHOCKTUBE_H

#include "Scenario.h"

namespace NavierStokes {
class DoubleShockTube : public Scenario {
  void initialValues(const double *const x, const NavierStokes &ns,
                     Variables &vars) final override;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_DOUBLESHOCKTUBE_H
