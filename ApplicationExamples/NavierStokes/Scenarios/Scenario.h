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

  virtual const double getGamma() const;
  virtual const double getPr() const;
  virtual const double getC_v() const;
  virtual const double getC_p() const;
  virtual const double getGasConstant() const;
  virtual const double getReferencePressure() const;

  const double gamma = 1.4;
  const double Pr = 0.7;
  const double c_v = 1;
  const double c_p = c_v * gamma;
  const double gasConstant = c_p - c_v;
  const double referencePressure = 10000;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_SCENARIO_H
