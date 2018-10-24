#ifndef COMPRESSIBLENAVIERSTOKES_ATMOSPHERE_H
#define COMPRESSIBLENAVIERSTOKES_ATMOSPHERE_H

#include "../PDE.h"

namespace NavierStokes {
double computeHydrostaticPressure(const PDE &ns, double g, double posZ,
                                  double backgroundTemperature);
double potentialTToT(const PDE &ns, double pressure, double potentialT);
}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_ATMOSPHERE_H
