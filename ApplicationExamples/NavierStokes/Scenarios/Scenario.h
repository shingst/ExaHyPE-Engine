#ifndef NAVIERSTOKES_SCENARIO_H
#define NAVIERSTOKES_SCENARIO_H

#include "../NavierStokes.h"
#include "tarch/la/Vector.h"
#include <array>

namespace NavierStokes {
enum class BoundaryType { analytical, wall };

class Scenario {
 public:
  Scenario();
  virtual void initialValues(const double* const x, const NavierStokes& ns,
                             Variables& vars);
  virtual void analyticalSolution(const double* const x, double t,
                                  const NavierStokes& ns, Variables& vars,
                                  double* gradState);
  virtual void source(const tarch::la::Vector<DIMENSIONS, double>& x, double t,
                      const NavierStokes& ns, const double* const Q, double* S);
  virtual BoundaryType getBoundaryType(int faceId);
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_SCENARIO_H
