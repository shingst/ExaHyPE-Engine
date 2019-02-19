#ifndef __EXAHYPE_TOVSOLVER_ADAPTER__
#define __EXAHYPE_TOVSOLVER_ADAPTER__

#include "InitialData.h"

namespace TOV {
	class TOVSolver; // Forward declaration
}

class TovSolverAdapter : public idobj {
	TOV::TOVSolver *tov;
public:
	TovSolverAdapter();
	void Interpolate(const double* const x, double t, double* const Q) override;
};


#endif /* __EXAHYPE_TOVSOLVER_ADAPTER__ */
