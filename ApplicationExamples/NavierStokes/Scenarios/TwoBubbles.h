#ifndef NAVIERSTOKES_TWOBUBBLES_H
#define NAVIERSTOKES_TWOBUBBLES_H

#include "Scenario.h"

namespace NavierStokes {

class TwoBubbles : public Scenario {
    struct Bubble {
        const double tempDifference;

        Bubble(const double tempDifference, const double size, const double decay, const double centerX,
               const double centerZ);

        // temp. difference
        const double size; // [m]
        const double decay; // [m]
        const double centerX; // [m]
        const double centerZ; // [m]
    };

public:
    void initialValues(const double *const x, const NavierStokes &ns, Variables &vars) override;

    void source(const double *const Q, double *S) override;
};
}

#endif //NAVIERSTOKES_TWOBUBBLES_H
