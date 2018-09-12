//
// Created by lukas on 27/08/18.
//

#ifndef NAVIERSTOKES_TAYLORGREEN_H
#define NAVIERSTOKES_TAYLORGREEN_H

#include "Scenario.h"

namespace NavierStokes {

class TaylorGreen : public Scenario {
  void initialValues(const double* const x, const NavierStokes& ns,
                     Variables& vars) final override;
  void analyticalSolution(const double* const x, const double t,
                          const NavierStokes& ns, Variables& vars,
                          double* gradState) final override;
  BoundaryType getBoundaryType(int faceId) final override;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_TAYLORGREEN_H
