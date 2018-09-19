#include "TovSolverAdapter.h"
#include "Fortran/PDE.h" // P2C

#include "GRMHDSolver_ADERDG_Variables.h"
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
using namespace GRMHD::GRMHDSolver_ADERDG_Variables::shortcuts;

// Trivial inlined interface to the C++ TovSolver which is included in the code

#include "../tovsolver_lib/tov.h"

TovSolverAdapter::TovSolverAdapter() {
	tov = new TOV::TOVSolver();
	
	// well, that could be then passed by the parameters one day, right...
	
	tov->TOV_Rho_Central[0]     = 1.28e-3 ;
	tov->TOV_Combine_Method = "maximum"   ;
	tov->TOV_Num_Radial     = 40000000    ;
	tov->TOV_dr[0]          = 0.00001     ;
	tov->Perturb[0]         = false       ;
	tov->Perturb_Pressure[0]   = true     ;
	tov->Pert_Press_Amplitude[0] = 0.01   ;

	tov->Setup();
}
	
void TovSolverAdapter::Interpolate(const double* x, double t, double* Q) {
	double V[nVar];
	
	TOV::idvars id;
	tov->Interpolate(x, id);
	
	V[rho] = id.rho;
	V[E] = id.press;
	for(int i=0; i<3; i++) {
		V[vel + i] = id.vel[i];
		V[B + i] = 0;
		V[shift + i] = id.beta[i];
	}
	V[psi] = 0;
	V[lapse] = id.alp;
	
	// Floor atmosphere
	if(V[rho] < Parameters::atmo_rho) V[rho] = Parameters::atmo_rho;
	if(V[E] < Parameters::atmo_press) V[E] = Parameters::atmo_press;
	
	// Caveats with the ordering, here it is Tensish (C)
	V[gij + 0] = id.gam[0][0];
	V[gij + 1] = id.gam[0][1];
	V[gij + 2] = id.gam[1][1];
	V[gij + 3] = id.gam[0][2];
	V[gij + 4] = id.gam[1][2];
	V[gij + 5] = id.gam[2][2];
	
	pdeprim2cons_(Q, V);

	//printf("At x=(%f,%f,%f), rho=%f\n", x[0],x[1],x[2],Q[rho]);
}
