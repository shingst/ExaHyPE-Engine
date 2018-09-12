#ifndef NAVIERSTOKES_CONVERGENCETEST_H
#define NAVIERSTOKES_CONVERGENCETEST_H

#include "../Scenario.h"

namespace NavierStokes {
class ConvergenceTest : public Scenario {
  void analyticalSolution(const double* const x, double t,
                          const NavierStokes& ns, Variables& vars,
                          double* gradState) override;
  void source(const tarch::la::Vector<DIMENSIONS, double>& x, double t,
              const NavierStokes& ns, const double* const Q,
              double* S) override;
  BoundaryType getBoundaryType(int faceId) override;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_CONVERGENCETEST_H
