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
  else
  {
	  enforceccz4constraints_(Q); 
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

bool isInRefinementZone(const tarch::la::Vector<DIMENSIONS,double>& center, double dx){
	double radius = 8.12514;
	// lower left, upper right radius of cell
	double cen = tarch::la::norm2(center);
	double dr = std::max(0.5,dx);
  bool shouldRefine = (cen > (radius -dr) ) && ( cen  <= (radius+dr) ); 
  return shouldRefine;
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

/*
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
			if(tarch::la::equals(t,0.0) && tarch::la::norm2(cellCentre)<100)
				return exahype::solvers::Solver::RefinementControl::Refine;
			//return exahype::solvers::Solver::RefinementControl::Recoarse;
			return exahype::solvers::Solver::RefinementControl::Keep;
		};
  }
*/

    if(isInRefinementZone(cellCentre,std::sqrt(tarch::la::norm2(cellSize))))
        return exahype::solvers::Solver::RefinementControl::Refine;
    return exahype::solvers::Solver::RefinementControl::Erase;


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
	   const int num = 3;
	  pdelimitervalue_(&limvalue,&cellCentre[0],&num, localDMPObservablesMin, localDMPObservablesMax);
	  return !limvalue;
}
void FOCCZ4::FOCCZ4Solver_ADERDG::fusedSource(const double* const restrict Q, const double* const restrict gradQ, double* const restrict S){
	//static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
  
  
    //tarch::multicore::Lock lock(initializationSemaphoreDG);	
	pdefusedsrcncp_(S,Q,gradQ);
	//fusedSource(Q, gradQ, S);
	//lock.free();
}


#ifdef CCZ4GRHD

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "kernels/KernelUtils.h"
#include "kernels/RiemannSolverUtils.h"

// from kernels:aderdg:generic:c:3d
void FOCCZ4::FOCCZ4Solver_ADERDG::MyNewRusanovSolver(
    double* FL, double* FR, const double* const QL, const double* const QR,
    const double t, const double dt,
    const tarch::la::Vector<DIMENSIONS, double>& dx, const int direction) {
  constexpr int numberOfVariables =
      FOCCZ4::AbstractFOCCZ4Solver_ADERDG::NumberOfVariables;
  constexpr int numberOfParameters =
      FOCCZ4::AbstractFOCCZ4Solver_ADERDG::NumberOfParameters;
  constexpr int numberOfData = numberOfVariables + numberOfParameters;
  constexpr int order = FOCCZ4::AbstractFOCCZ4Solver_ADERDG::Order;
  constexpr int basisSize = order + 1;

  // Compute the average variables and parameters from the left and the right
  double QavL[numberOfData] = {0.0};
  double QavR[numberOfData] = {0.0};
  kernels::riemannsolvers::util::averageRiemannInputs<basisSize, numberOfData>(
      QL, FOCCZ4::AbstractFOCCZ4Solver_ADERDG::weights, QavL);
  kernels::riemannsolvers::util::averageRiemannInputs<basisSize, numberOfData>(
      QR, FOCCZ4::AbstractFOCCZ4Solver_ADERDG::weights, QavR);

  double LL[numberOfVariables] = {
      0.0};  // do not need to store material parameters
  double LR[numberOfVariables] = {0.0};

  // Hyperbolic eigenvalues
  eigenvalues(QavL, direction, LL);
  eigenvalues(QavR, direction, LR);
  // skip parameters
  std::transform(LL, LL + numberOfVariables, LL, std::abs<double>);
  std::transform(LR, LR + numberOfVariables, LR, std::abs<double>);
  const double maxHyperbolicEigenvalueL =
      *std::max_element(LL, LL + numberOfVariables);
  const double maxHyperbolicEigenvalueR =
      *std::max_element(LR, LR + numberOfVariables);
  const double maxHyperbolicEigenvalue =
      std::max(maxHyperbolicEigenvalueL, maxHyperbolicEigenvalueR);

  double smax = maxHyperbolicEigenvalue;

  double smaxtwo;
  double smaxtmpone;
  double smaxtmptwo;
  smaxtwo = 0.;

  double nv[3] = {0.};
  nv[direction] = 1;

  {
    kernels::idx3 idx_QLR(basisSize, basisSize, numberOfData);
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        smaxtmpone = smaxtwo;
        // pdeeigenvalues_(lambda, Q, nv);
        maxhyperboliceigenvaluegrhd_(&smaxtmptwo, &QR[idx_QLR(i, j, 0)],
                                     &QL[idx_QLR(i, j, 0)], nv);
        smaxtwo = std::max(smaxtmpone, smaxtmptwo);
      }
    }
  }

  // compute fluxes (and fluctuations for non-conservative PDEs)
  double Qavg[numberOfData];
  kernels::idx2 idx_gradQ(DIMENSIONS, numberOfVariables);
  double gradQ[DIMENSIONS][numberOfVariables] = {0.0};
  double ncp[numberOfVariables] = {0.0};
  {
    kernels::idx3 idx_FLR(basisSize, basisSize, numberOfVariables);
    kernels::idx3 idx_QLR(basisSize, basisSize, numberOfData);
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        // if (useNCP) {  // we don't use matrixB but the NCP call here.
        for (int l = 0; l < numberOfVariables; l++) {
          gradQ[direction][l] = QR[idx_QLR(i, j, l)] - QL[idx_QLR(i, j, l)];
        }
        for (int l = 0; l < numberOfData; l++) {
          Qavg[l] = 0.5 * (QR[idx_QLR(i, j, l)] + QL[idx_QLR(i, j, l)]);
        }

        nonConservativeProduct(Qavg, gradQ[0], ncp);
        //}

        // skip parameters
        for (int k = 0; k < 59; k++) {
          FL[idx_FLR(i, j, k)] =
              0.5 * (FR[idx_FLR(i, j, k)] + FL[idx_FLR(i, j, k)]) -
              0.5 * smax * (QR[idx_QLR(i, j, k)] - QL[idx_QLR(i, j, k)]);

          // if (useNCP) {
          FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] - 0.5 * ncp[k];
          FL[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] + 0.5 * ncp[k];
          //} else {
          //  FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)];
          //}
        }
        // this is the hydro
        for (int k = 59; k < 64; k++) {
          FL[idx_FLR(i, j, k)] =
              0.5 * (FR[idx_FLR(i, j, k)] + FL[idx_FLR(i, j, k)]) -
              0.5 * smaxtwo * (QR[idx_QLR(i, j, k)] - QL[idx_QLR(i, j, k)]);

          // if (useNCP) {
          FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] - 0.5 * ncp[k];
          FL[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] + 0.5 * ncp[k];
          //} else {
          //  FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)];
          //}
        }
        // this is the GLM cleaning
        for (int k = 64; k < numberOfVariables; k++) {
          FL[idx_FLR(i, j, k)] =
              0.5 * (FR[idx_FLR(i, j, k)] + FL[idx_FLR(i, j, k)]) -
              0.5 * smax * (QR[idx_QLR(i, j, k)] - QL[idx_QLR(i, j, k)]);

          // if (useNCP) {
          FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] - 0.5 * ncp[k];
          FL[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)] + 0.5 * ncp[k];
          //} else {
          //  FR[idx_FLR(i, j, k)] = FL[idx_FLR(i, j, k)];
          //}
        }
      }
    }
  }
}

#endif

 