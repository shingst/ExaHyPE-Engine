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
	void Interpolate(const double* x, double t, double* Q) override;
};


#endif /* __EXAHYPE_TOVSOLVER_ADAPTER__ */
