#include "TovSolverAdapter.h"

#include "PDE/PDE.h" // P2C
using namespace SVEC::GRMHD;

#include "GRMHDSolver_ADERDG_Variables.h"
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
using namespace GRMHD::GRMHDSolver_ADERDG_Variables::shortcuts;

// Trivial inlined interface to the C++ TovSolver which is included in the code

#include "../tovsolver_lib/tov.h"

TovSolverAdapter::TovSolverAdapter() : hasBeenPrepared(false) {
	tov = new TOV::TOVSolver();
}
	
void TovSolverAdapter::readParameters(const mexa::mexafile& parameters) {
	// well, that could be then passed by the parameters, right...
	tov->TOV_Rho_Central[0]     = 1.28e-3 ;
	tov->TOV_Combine_Method = "maximum"   ;
	tov->TOV_Num_Radial     = 40000000    ;
	tov->TOV_dr[0]          = 0.00001     ;
	tov->Perturb[0]         = false       ;
	tov->Perturb_Pressure[0]   = true     ;
	tov->Pert_Press_Amplitude[0] = 0.01   ;
}
	
void TovSolverAdapter::prepare() {
	tov->Setup();
	hasBeenPrepared = true;
}
	
void TovSolverAdapter::Interpolate(const double* x, double t, double* Q) {
	if(!hasBeenPrepared) std::abort();

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
	
	Prim2Cons(Q, V).copyFullStateVector();
	//printf("At x=(%f,%f,%f), rho=%f\n", x[0],x[1],x[2],Q[rho]);
}
