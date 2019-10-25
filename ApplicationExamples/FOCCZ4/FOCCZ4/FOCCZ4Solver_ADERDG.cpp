// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "FOCCZ4Solver_ADERDG.h"
#include "FOCCZ4Solver_FV.h"

#include <algorithm>

#include "FOCCZ4Solver_ADERDG_Variables.h"
#include "kernels/GaussLegendreBasis.h"

#include "kernels/KernelUtils.h"
#include "peano/utils/Loop.h"

#include "PDE.h"
#include "InitialData.h"
#include "Tools.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"
tarch::logging::Log FOCCZ4::FOCCZ4Solver_ADERDG::_log( "FOCCZ4::FOCCZ4Solver_ADERDG" );


void FOCCZ4::FOCCZ4Solver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {

    const int order = FOCCZ4::FOCCZ4Solver_ADERDG::Order;
  	constexpr int basisSize = AbstractFOCCZ4Solver_FV::PatchSize;
	constexpr int Ghostlayers = AbstractFOCCZ4Solver_FV::GhostLayerWidth;
	
	static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
  
  
    tarch::multicore::Lock lock(initializationSemaphoreDG);	

  if (constants.isValueValidString("reference")) {
    std::string reference = constants.getValueAsString("reference");
	const int length=reference.length();
	logInfo("init(...)","Reference setup:"<<reference);
  	printf("\n******************************************************************");
	printf("\n**************<<<  INIT TECPLOT    >>>****************************");
	printf("\n******************************************************************");
    inittecplot_(&order,&order,&basisSize,&Ghostlayers);
	//inittecplot_(&order,&order);
	printf("\n******************************************************************");
	printf("\n**************<<<  INIT PDE SETUP  >>>****************************");
	printf("\n******************************************************************\n");
    initparameters_(&length,&reference[0]);
	printf("\n******************************************************************");
	printf("\n**************<<<       DONE       >>>****************************");
	printf("\n******************************************************************");
  } else {
    logInfo("init(...)","Not recognized setup.");
	std::abort();
  }	
	
	lock.free();
    
}

void FOCCZ4::FOCCZ4Solver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = FOCCZ4::FOCCZ4Solver_ADERDG::Order;
    std::fill_n(Q,96,0.0);

    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Q);
  }
  for(int i = 0; i< 96 ; i++){
    assert(std::isfinite(Q[i]));
  }
}

void FOCCZ4::FOCCZ4Solver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  const int nVar = FOCCZ4::FOCCZ4Solver_ADERDG::NumberOfVariables;
  const int order = FOCCZ4::FOCCZ4Solver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;
  double Qgp[nVar],*F[nDim], Fs[nDim][nVar];

  double x_3[3];
  x_3[2]=0;
  std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
  
  int md=0;
  double cms=0;
	
  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut , 0, nVar * sizeof(double));
	
  //std::copy_n(stateIn,nVar,stateOut);
  //std::copy_n(fluxIn,nVar,fluxOut);
  for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
  
  for(int i=0; i < basisSize; i++)  { // i == time
    const double weight = kernels::legendre::weights[order][i];
    const double xi = kernels::legendre::nodes[order][i];
    double ti = t + xi * dt;
  
    //    initialdata_(x, &ti, Qgp,&md,&cms,&order);
    initialdata_(x_3, &ti, Qgp);
    flux(Qgp, F);
    for(int m=0; m < nVar; m++) {
      stateOut[m] += weight * Qgp[m];
      fluxOut[m] += weight * Fs[direction][m];
    }
  }
  
}

exahype::solvers::Solver::RefinementControl FOCCZ4::FOCCZ4Solver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,double t,const int level) {
  // @todo Please implement/augment if required
  //return exahype::solvers::Solver::RefinementControl::Keep;
  /*if(DIMENSIONS == 2){
    if(std::abs(cellCentre[0]) < 10){
      if(std::abs(cellCentre[1]) < 10){
	return exahype::solvers::Solver::RefinementControl::Refine;
      }
    }
  }else{
    if(std::abs(cellCentre[0]) < 5){
      if(std::abs(cellCentre[1]) < 5){
		  if(std::abs(cellCentre[2]) < 5){
			return exahype::solvers::Solver::RefinementControl::Refine;
		  }
      }
    }	  
	  
  };	 
  
  //return exahype::solvers::Solver::RefinementControl::Keep;
  if ( level > getCoarsestMeshLevel() ) {
    return exahype::solvers::Solver::RefinementControl::Erase;
  }
return exahype::solvers::Solver::RefinementControl::Keep;*/


  const int nVar = FOCCZ4::AbstractFOCCZ4Solver_ADERDG::NumberOfVariables;
  const int order = FOCCZ4::AbstractFOCCZ4Solver_ADERDG::Order;
  const int basisSize = order + 1;
  int refine_flag;
  double max_luh[nVar];
  double min_luh[nVar];
  
  //return exahype::solvers::Solver::RefinementControl::Keep;
 
 // 	 
 
  for(int m = 0; m <nVar ; m++){
	max_luh[m]=-1.e+14;
	min_luh[m]=1.e+14;
  }
  
#if DIMENSIONS==3
	kernels::idx4 id_xyz_dof(basisSize,basisSize,basisSize,nVar);
#else
	kernels::idx3 id_xy_dof(basisSize,basisSize,nVar);
#endif
  
  for(int i = 0; i < basisSize; i++){
		for(int j = 0; j <basisSize ; j++){
#if DIMENSIONS==3	
			for(int k = 0; k <basisSize ; k++){
#endif
#if DIMENSIONS==3
				for(int m = 0; m <nVar ; m++){
					if(luh[id_xyz_dof(i,j,k,m)]<min_luh[m]){
						min_luh[m]=luh[id_xyz_dof(i,j,k,m)];
					}
					if(luh[id_xyz_dof(i,j,k,m)]>max_luh[m]){
						max_luh[m]=luh[id_xyz_dof(i,j,k,m)];
					}
				}
#else
				for(int m = 0; m <nVar ; m++){
					if(luh[id_xy_dof(i,j,m)]<min_luh[m]){
						min_luh[m]=luh[id_xy_dof(i,j,m)];
					}
					if(luh[id_xy_dof(i,j,m)]>max_luh[m]){
						max_luh[m]=luh[id_xy_dof(i,j,m)];
					}
				}	
#endif
		

#if DIMENSIONS==3				
			}
#endif			
		}
	}		  
  pderefinecriteria_(&refine_flag,&max_luh[0],&min_luh[0],&cellCentre[0]);
  if(refine_flag>1){
	  return exahype::solvers::Solver::RefinementControl::Refine;
  }else{
		if(refine_flag>0){
			return exahype::solvers::Solver::RefinementControl::Keep;
		}else{
			//return exahype::solvers::Solver::RefinementControl::Recoarse;
			return exahype::solvers::Solver::RefinementControl::Keep;
		};
  }


}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void FOCCZ4::FOCCZ4Solver_ADERDG::eigenvalues(const double* const Q,const int direction,double* const lambda) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_ADERDG.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[direction] = 1;
  pdeeigenvalues_(lambda, Q, nv);

  for(int i = 0; i< 96 ; i++){
    assert(std::isfinite(lambda[i]));
  }
}





void FOCCZ4::FOCCZ4Solver_ADERDG::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_ADERDG.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  constexpr int nVar = FOCCZ4::FOCCZ4Solver_ADERDG::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }

  for(int d = 0; d< DIMENSIONS ; d++){
    for(int i = 0; i< 96 ; i++){
      assert(std::isfinite(F[d][i]));
    }
  }
}


//You can either implement this method or modify fusedSource
void FOCCZ4::FOCCZ4Solver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_ADERDG.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  // @todo Please implement/augment if required
   pdesource_(S, Q);

  for(int i = 0; i< 96 ; i++){
    assert(std::isfinite(S[i]));
  }
}

void  FOCCZ4::FOCCZ4Solver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_ADERDG.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_ADERDG.h" for toolkit generated compile-time 
  //      constants such as Order, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
   pdencp_(BgradQ, Q, gradQ);
  for(int i = 0; i< 96 ; i++){
    assert(std::isfinite(BgradQ[i]));
  }
}



bool FOCCZ4::FOCCZ4Solver_ADERDG::isPhysicallyAdmissible(
      const double* const                         solution,
      const double* const                         localDMPObservablesMin,
      const double* const                         localDMPObservablesMax,
      const bool                                  wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const double                               timeStamp) const
	  {
  		  
	  int limvalue;
	   
	  pdelimitervalue_(&limvalue,&cellCentre[0],&NumberOfDMPObservables, localDMPObservablesMin, localDMPObservablesMax);
	  bool ret_value;
	  limvalue > 0 ? ret_value=false : ret_value=true;
	  return ret_value;
}
void FOCCZ4::FOCCZ4Solver_ADERDG::fusedSource(const double* const restrict Q, const double* const restrict gradQ, double* const restrict S){
	//static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
  
  
    //tarch::multicore::Lock lock(initializationSemaphoreDG);	
	pdefusedsrcncp_(S,Q,gradQ);
	//fusedSource(Q, gradQ, S);
	//lock.free();
}
