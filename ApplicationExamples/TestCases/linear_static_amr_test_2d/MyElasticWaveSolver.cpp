// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
// 
// ========================
//   www.exahype.eu
// ========================

#include "MyElasticWaveSolver.h"

#include "MyElasticWaveSolver_Variables.h"

#include "../../../ExaHyPE/kernels/KernelUtils.h"
#include "../../../ExaHyPE/kernels/GaussLegendreQuadrature.h"

#ifdef OPT_KERNELS
#include "kernels/MyElasticWaveSolver/converter.h"
using namespace Elastic::MyElasticWaveSolver_kernels::aderdg;
#endif



tarch::logging::Log Elastic::MyElasticWaveSolver::_log( "Elastic::MyElasticWaveSolver" );


void Elastic::MyElasticWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required

  initPointSourceLocations(cmdlineargs,constants);
}

void Elastic::MyElasticWaveSolver::adjustSolution(double *luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  // @todo Please implement/augment if required

  int level=std::round(log(_domainSize[0]/dx[0])/log(3.)) + 1;
  if (tarch::la::equals(t,0.0)) {
    
    constexpr int basisSize = MyElasticWaveSolver::Order+1;
    constexpr int numberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

    kernels::idx3 id_xyzf(basisSize,basisSize,numberOfData);
    kernels::idx2 id_xyz(basisSize,basisSize);

    int num_nodes = basisSize;

    double offset_x=center[0]-0.5*dx[0];
    double offset_y=center[1]-0.5*dx[1];

    double width_x=dx[0];
    double width_y=dx[1];

    for (int j=0; j< num_nodes; j++){
      for (int i=0; i< num_nodes; i++){
	double x  =  (offset_x+width_x*kernels::gaussLegendreNodes[basisSize-1][i]);
	double y  =  (offset_y+width_y*kernels::gaussLegendreNodes[basisSize-1][j]);


	// velocity
	luh[id_xyzf(j,i,0)]  = 0;
	luh[id_xyzf(j,i,1)]  = 0;
	
	// stress field
	luh[id_xyzf(j,i,2)]  = 0;
	luh[id_xyzf(j,i,3)]  = 0;
	luh[id_xyzf(j,i,4)]  = 0;

	  if( level <= getCoarsestMeshLevel()){	  
	    // material parameters for loh.1
	    luh[id_xyzf(j,i,5)]  = x*x+2; //rho
	    luh[id_xyzf(j,i,6)] = y*y+2; //cp
	    luh[id_xyzf(j,i,7)] = 1; //cs
	  }
	}
      }
    }

}

void Elastic::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  const int NumberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

  for (int i = 0; i<NumberOfData; i++){
    stateOut[i] = stateIn[i];
  }
 
  for (int i = 0; i< NumberOfVariables; i++){
    fluxOut[i] =  fluxIn[i];
  }
}

exahype::solvers::Solver::RefinementControl Elastic::MyElasticWaveSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  
  if (level <= getCoarsestMeshLevel()){
    //std::cout << "erase" << std::endl;
    return exahype::solvers::Solver::RefinementControl::Refine;
  }

  //std::cout << "keep" << std::endl;
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Elastic::MyElasticWaveSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  
  double cp = Q[6];
  double cs = Q[7];
   

  lambda[0] = cp;
  lambda[1] = cs;
  lambda[2] = -cp;
  lambda[3] = -cs;
  lambda[4] = 0.0;
}


void Elastic::MyElasticWaveSolver::flux(const double* const Q,double** F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3

  double sxx = Q[2];
  double syy = Q[3];
  double sxy = Q[4];
  
  
  F[0][0] = -sxx;
  F[0][1] = -sxy;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  
  F[1][0] = -sxy;
  F[1][1] = -syy;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
}


void  Elastic::MyElasticWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  double vx_x = gradQ[0];
  double vy_x = gradQ[1];

  double vx_y = gradQ[5];
  double vy_y = gradQ[6];

  
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = -vx_x;
  BgradQ[3] = 0.0;
  BgradQ[4] = -vy_x;

  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = -vy_y;
  BgradQ[9] = -vx_y;

}


void  Elastic::MyElasticWaveSolver::initPointSourceLocations(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants){
  pointSourceLocation[0][0]=0.0;
  pointSourceLocation[0][1]=2.0;
  pointSourceLocation[0][2]=0.0;

}

void  Elastic::MyElasticWaveSolver::pointSource(const double* const Q,const double* const x,const double t,const double dt, double* forceVector, int n) {
  constexpr double t0 = 0.1;
  constexpr double M0 = 1000.0;

  if(n==0){
    double f = M0*t/(t0*t0)*std::exp(-t/t0);
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
    forceVector[4] = 0.0;
    forceVector[5] = 0.0;
    forceVector[6] = 0.0;
    forceVector[7] = f;
  }
}

/**
 * @TODO LR : document
 */
void Elastic::MyElasticWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* rhs) {



  const double rho  = Q[5];   // km/s
  const double cp   = Q[6];   // km/s
  const double cs   = Q[7];   // km/s
  //double jacobian = Q[8];
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp-2.0*mu;


  rhs[0] = 1.0/rho*rhs[0];
  rhs[1] = 1.0/rho*rhs[1];
  
  double rhs_2= (2*mu+lam)*rhs[2]+lam*rhs[3];
  double rhs_3= (2*mu+lam)*rhs[3]+lam*rhs[2];
  
  rhs[2]=rhs_2;
  rhs[3]=rhs_3;  
  rhs[4]=mu*rhs[4];


  rhs[5] = 1.0/rho*rhs[5];
  rhs[6] = 1.0/rho*rhs[6];
  
  double rhs_7= (2*mu+lam)*rhs[7]+lam*rhs[8];
  double rhs_8= (2*mu+lam)*rhs[8]+lam*rhs[7];
  
  rhs[7]=rhs_7;
  rhs[8]=rhs_8;  
  rhs[9]=mu*rhs[9];
}
