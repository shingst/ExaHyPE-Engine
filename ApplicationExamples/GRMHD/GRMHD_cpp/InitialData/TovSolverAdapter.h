#ifndef __EXAHYPE_TOVSOLVER_ADAPTER__
#define __EXAHYPE_TOVSOLVER_ADAPTER__

#include "InitialData.h"

#include "GRMHDSolver_ADERDG_Variables.h"

// Trivial inlined interface to the C++ TovSolver which is included in the code

#include "../tovsolver_lib/tov.h"

class TovSolverAdapter : public InitialDataCode {
	TOV::TOVSolver tov;
public:
	bool hasBeenPrepared; ///< a guard to ensure the preparation took place
	
	TovSolverAdapter() : hasBeenPrepared(false) {}
	
	void readParameters(const mexa::mexafile& parameters) override {
		// well, that could be then passed by the parameters, right...
		tov.TOV_Rho_Central[0]     = 1.28e-3 ;
		tov.TOV_Combine_Method = "maximum"   ;
		tov.TOV_Num_Radial     = 40000000    ;
		tov.TOV_dr[0]          = 0.00001     ;
		tov.Perturb[0]         = false       ;
		tov.Perturb_Pressure[0]   = true     ;
		tov.Pert_Press_Amplitude[0] = 0.01   ;
	}
	
	void prepare() override {
		tov.Setup();
		hasBeenPrepared = true;
	}
	
	void Interpolate(const double* x, double t, double* Q) override {
		TOV::idvars id;
		tov.Interpolate(x, id);
		
		using namespace GRMHD::GRMHDSolver_ADERDG_Variables::shortcuts;
		
		Q[rho] = id.rho;
		Q[E] = id.eps;
		for(int i=0; i<3; i++) {
			Q[vel + i] = id.vel[i];
			Q[B + i] = 0;
			Q[shift + i] = id.beta[i];
		}
		Q[psi] = 0;
		Q[lapse] = id.alp;
		
		// Caveats with the ordering, here it is Tensish (C)
		Q[gij + 0] = id.gam[0][0];
		Q[gij + 1] = id.gam[0][1];
		Q[gij + 2] = id.gam[1][1];
		Q[gij + 3] = id.gam[0][2];
		Q[gij + 4] = id.gam[1][2];
		Q[gij + 5] = id.gam[2][2];
	}
};


#endif /* __EXAHYPE_TOVSOLVER_ADAPTER__ */
