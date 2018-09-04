#ifndef __EXAHYPE_TOVSOLVER_ADAPTER__
#define __EXAHYPE_TOVSOLVER_ADAPTER__

#include "InitialData.h"

namespace TOV {
	class TOVSolver; // Forward declaration
}

class TovSolverAdapter : public InitialDataCode {
	TOV::TOVSolver *tov;
public:
	bool hasBeenPrepared; ///< a guard to ensure the preparation took place
	
	TovSolverAdapter();
	
	void readParameters(const mexa::mexafile& parameters) override;
	void prepare() override;
	void Interpolate(const double* x, double t, double* Q) override;
};


#endif /* __EXAHYPE_TOVSOLVER_ADAPTER__ */
