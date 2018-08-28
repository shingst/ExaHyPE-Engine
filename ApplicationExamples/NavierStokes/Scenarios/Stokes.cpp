#include "Stokes.h"

void NavierStokes::Stokes::initialValues(const double* const x, const NavierStokes& ns, Variables& vars) {
    vars.rho() = 1.0;
#if DIMENSIONS == 2
    vars.j(1.0, 0.0);
#elif DIMENSIONS == 3
    vars.j(1.0, 0.0, 0.0);
#endif
    const double pressure = 100/ns.GAMMA;
    vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

}
void NavierStokes::Stokes::analyticalSolution(const double* const x, const double t,
                        const NavierStokes& ns, Variables& vars,double* gradState) {
      kernels::idx2 idxGradQ(DIMENSIONS,vars.SizeVariables);
      const double rho = 1;
      const double v = ns.referenceViscosity / rho;
      vars.rho() = 1.0;
      vars.j(0) = std::erf(x[1] / std::sqrt(2 * v * t));
      vars.j(1) = 0.0;
#if DIMENSIONS == 3
      vars.j(2) = 0.0;
#endif
      const double pressure = 100/ns.GAMMA;
      vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

      // Assuming rho is constant.
      gradState[idxGradQ(1,1)] = 1.0;
#if DIMENSIONS == 2
      gradState[idxGradQ(1,3)] = 1/vars.rho() * vars.j(0);
#elif DIMENSIONS == 3
      gradState[idxGradQ(1,4)] = 1/vars.rho() * vars.j(0);
#endif

}

NavierStokes::BoundaryType NavierStokes::Stokes::getBoundaryType(int faceId) {
    if (faceId == 2) {
        return BoundaryType::analytical;
    }
    return BoundaryType::wall;
}
