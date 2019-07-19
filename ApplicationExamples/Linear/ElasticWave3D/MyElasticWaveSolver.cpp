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
#include "../../../ExaHyPE/kernels/GaussLegendreBasis.h"

#ifdef OPT_KERNELS
#include "kernels/Elastic_MyElasticWaveSolver/converter.h"
using namespace Elastic::MyElasticWaveSolver_kernels::aderdg;
#endif

#include <algorithm>



tarch::logging::Log Elastic::MyElasticWaveSolver::_log( "Elastic::MyElasticWaveSolver" );


void Elastic::MyElasticWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
  initPointSourceLocations(cmdlineargs,constants);
  amr_regularization = false; //constants.getValueAsBool ("amr_regularization");
}

void Elastic::MyElasticWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    
    constexpr int basisSize = MyElasticWaveSolver::Order+1;
    constexpr int numberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

    kernels::idx4 id_xyzf(basisSize,basisSize,basisSize,numberOfData);
    kernels::idx3 id_xyz(basisSize,basisSize,basisSize);

    int num_nodes = basisSize;

    double offset_x=center[0]-0.5*dx[0];
    double offset_y=center[1]-0.5*dx[1];
    double offset_z=center[2]-0.5*dx[2];  

    double width_x=dx[0];
    double width_y=dx[1];
    double width_z=dx[2];

    for (int k=0; k< num_nodes; k++){
      for (int j=0; j< num_nodes; j++){
	for (int i=0; i< num_nodes; i++){

	  double x  =  (offset_x+width_x*kernels::legendre::nodes[basisSize-1][i]);
	  double y  =  (offset_y+width_y*kernels::legendre::nodes[basisSize-1][j]);
	  double z  =  (offset_z+width_z*kernels::legendre::nodes[basisSize-1][k]);

	  // velocity
	  luh[id_xyzf(k,j,i,0)]  = 0;
	  luh[id_xyzf(k,j,i,1)]  = 0;
	  luh[id_xyzf(k,j,i,2)]  = 0;
	
	  // stress field
	  luh[id_xyzf(k,j,i,3)]  = 0;
	  luh[id_xyzf(k,j,i,4)]  = 0;
	  luh[id_xyzf(k,j,i,5)]  = 0;
	  luh[id_xyzf(k,j,i,6)]  = 0;
	  luh[id_xyzf(k,j,i,7)]  = 0;
	  luh[id_xyzf(k,j,i,8)]  = 0;


	  // material parameters for loh.1
	  luh[id_xyzf(k,j,i,9)]  = 2.7; //rho
	  luh[id_xyzf(k,j,i,10)] = 6.0; //cp
	  luh[id_xyzf(k,j,i,11)] = 3.343; //cs
	  
	  if( y < 1.0) {
	    luh[id_xyzf(k,j,i,9)]  = 2.6; //rho
	    luh[id_xyzf(k,j,i,10)] = 4.0; //cp
	    luh[id_xyzf(k,j,i,11)] = 2.0; //cs
	  }

	}
      }
    }
    
  }
}

void Elastic::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
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

exahype::solvers::Solver::RefinementControl Elastic::MyElasticWaveSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  
  double left_vertex[3];
  double right_vertex[3];

  for(int i = 0 ; i<DIMENSIONS; i++){
    left_vertex[i]  = center[i] - dx[i]*0.5;
    right_vertex[i] = center[i] + dx[i]*0.5;
  }

  bool pointSourceInElement= true;
  for(int i = 0 ; i<DIMENSIONS; i++){
    pointSourceInElement &= ((left_vertex[i] <= pointSourceLocation[0][i]) && (right_vertex[i] >= pointSourceLocation[0][i]));
  }

  double refined_domain;
  if(amr_regularization){
    refined_domain = getCoarsestMeshSize();
  }else{
    refined_domain = std::ceil(pointSourceLocation[0][1]/getCoarsestMeshSize())*getCoarsestMeshSize();
  }

  bool elementOnSurface = left_vertex[1] < refined_domain;
  if (tarch::la::equals(t,0.0)) {
    if(pointSourceInElement || elementOnSurface){
      return exahype::solvers::Solver::RefinementControl::Refine;
    }else{
      return exahype::solvers::Solver::RefinementControl::Keep;
    }
  }

 if(pointSourceInElement){
   if(t <2.25){
     return exahype::solvers::Solver::RefinementControl::Keep;
   }
 }



 //compute velocity vector length
//  kernels::idx4 id((Order+1),(Order+1),(Order+1),(NumberOfParameters+NumberOfVariables));
//  double max_velocity=0;
//  double min_velocity=std::numeric_limits<double>::max();
//  for(int k = 0 ; k<(Order+1); k++){
//    for(int j = 0 ; j<(Order+1); j++){
//      for(int i = 0 ; i<(Order+1); i++){
//        double u=luh[id(k,j,i,0)];
//        double v=luh[id(k,j,i,1)];
//        double w=luh[id(k,j,i,2)];
//        double velocity=sqrt(u*u+v*v+w*w);
//        max_velocity=std::max(max_velocity,velocity);
//        min_velocity=std::min(min_velocity,velocity);
//      }
//    }
//  }
 
 /* constexpr double wave_threshold=0.5;
 // double wave_threshold= std::sqrt(center[0]*center[0]+center[1]*center[1]) * 10.0;
     
 if(max_velocity > wave_threshold){
   return exahype::solvers::Solver::RefinementControl::Refine;
 }

 if(max_velocity < wave_threshold*0.1){
   return exahype::solvers::Solver::RefinementControl::Erase;
   }*/
 
 return exahype::solvers::Solver::RefinementControl::Keep;
}


//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Elastic::MyElasticWaveSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  
  double cp = Q[10];
  double cs = Q[11];
   
  lambda[0] = cp;
  lambda[1] = cs;
  lambda[2] = cs;
  
  lambda[3] = -cp;
  lambda[4] = -cs;
  lambda[5] = -cs;
  
  lambda[6] = 0.0;
  lambda[7] = 0.0;
  lambda[8] = 0.0;
}


void Elastic::MyElasticWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 3
  double sigma_xx=Q[3];
  double sigma_yy=Q[4];
  double sigma_zz=Q[5];  
  double sigma_xy=Q[6];
  double sigma_xz=Q[7];
  double sigma_yz=Q[8];

  F[0][ 0] = -sigma_xx;
  F[0][ 1] = -sigma_xy;
  F[0][ 2] = -sigma_xz;
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;
  F[0][ 6] = 0.0;
  F[0][ 7] = 0.0;
  F[0][ 8] = 0.0;
  
  F[1][ 0] = -sigma_xy;
  F[1][ 1] = -sigma_yy;
  F[1][ 2] = -sigma_yz;
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
  F[1][ 5] = 0.0;
  F[1][ 6] = 0.0;
  F[1][ 7] = 0.0;
  F[1][ 8] = 0.0;

  F[2][ 0] = -sigma_xz;
  F[2][ 1] = -sigma_yz;
  F[2][ 2] = -sigma_zz;
  F[2][ 3] = 0.0;
  F[2][ 4] = 0.0;
  F[2][ 5] = 0.0;
  F[2][ 6] = 0.0;
  F[2][ 7] = 0.0;
  F[2][ 8] = 0.0;

}



void  Elastic::MyElasticWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  double u_x = gradQ[0];
  double v_x = gradQ[1];
  double w_x = gradQ[2];

  double u_y = gradQ[9];
  double v_y = gradQ[10];
  double w_y = gradQ[11];

  double u_z = gradQ[18];
  double v_z = gradQ[19];
  double w_z = gradQ[20];
 
  BgradQ[0] = 0;
  BgradQ[1] = 0;
  BgradQ[2] = 0;  
  BgradQ[3] =-u_x;
  BgradQ[4] =-0;
  BgradQ[5] =-0;
  BgradQ[6] =-v_x; //sigma_xy
  BgradQ[7] =-w_x; //sigma_xz
  BgradQ[8] =-0.0; //sigma_yz

  BgradQ[9] = 0;
  BgradQ[10]= 0;
  BgradQ[11]= 0;  
  BgradQ[12]= 0;
  BgradQ[13]= -v_y;
  BgradQ[14]= 0;
  BgradQ[15]= -u_y; //sigma_xy
  BgradQ[16]=- 0; //sigma_xz
  BgradQ[17]=-w_y; //sigma_yz

  BgradQ[18]= 0;
  BgradQ[19]= 0;
  BgradQ[20]= 0;  
  BgradQ[21]= 0;
  BgradQ[22]= 0;
  BgradQ[23]=- w_z;
  BgradQ[24]=-0; //sigma_xy
  BgradQ[25]=-u_z; //sigma_xz
  BgradQ[26]=-v_z; //sigma_yz    
}


void  Elastic::MyElasticWaveSolver::initPointSourceLocations(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants){
  pointSourceLocation[0][0]=0.0;
  pointSourceLocation[0][1]=2.0;
  pointSourceLocation[0][2]=0.0;

}


//Point source taken from: http://www.sismowine.org/model/WP2_LOH1.pdf
void  Elastic::MyElasticWaveSolver::pointSource(const double* const Q,const double* const x,const double t,const double dt, double* const forceVector, int n) {
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
    forceVector[8] = 0.0;
  }
}

/**
 * @TODO LR : document
 */
void Elastic::MyElasticWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
   double rho = Q[9];  
  double c_p = Q[10];
  double c_s = Q[11];  
  double mu     = rho*c_s*c_s;
  double lambda = rho*c_p*c_p-2*mu;
  double rho_inv=1.0/rho;

  rhs[0]=rho_inv * rhs[0];
  rhs[1]=rho_inv * rhs[1];
  rhs[2]=rho_inv * rhs[2];
  
  double lam_temp = lambda * (rhs[3] + rhs[4] + rhs[5]);
  rhs[3]=(2*mu) * rhs[3] +lam_temp;
  rhs[4]=(2*mu) * rhs[4] +lam_temp;
  rhs[5]=(2*mu) * rhs[5] +lam_temp;

  rhs[6]= mu*rhs[6];
  rhs[7]= mu*rhs[7];
  rhs[8]= mu*rhs[8];


  rhs[9] =rho_inv * rhs[9];
  rhs[10]=rho_inv * rhs[10];
  rhs[11]=rho_inv * rhs[11];

  
  lam_temp = lambda * (rhs[12] + rhs[13] + rhs[14]);

  rhs[12]=(2*mu) * rhs[12] +lam_temp;
  rhs[13]=(2*mu) * rhs[13] +lam_temp;
  rhs[14]=(2*mu) * rhs[14] +lam_temp;

  rhs[15]= mu*rhs[15];
  rhs[16]= mu*rhs[16];
  rhs[17]= mu*rhs[17];

  rhs[18] = rho_inv * rhs[18];
  rhs[19] = rho_inv * rhs[19];
  rhs[20] = rho_inv * rhs[20];
  
  lam_temp = lambda * (rhs[21] + rhs[22] + rhs[23]);

  rhs[21]=(2*mu) * rhs[21] +lam_temp;
  rhs[22]=(2*mu) * rhs[22] +lam_temp;
  rhs[23]=(2*mu) * rhs[23] +lam_temp;

  rhs[24]= mu*rhs[24];
  rhs[25]= mu*rhs[25];
  rhs[26]= mu*rhs[26];
}

void Elastic::MyElasticWaveSolver::riemannSolver(double* const FL_,double* const FR_,const double* const QL_,const double* const QR_,const double t,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex){

#ifdef OPT_KERNELS
  double FL[converter::getFFaceGenArraySize()];
  double FR[converter::getFFaceGenArraySize()];
  double QL[converter::getQFaceGenArraySize()];
  double QR[converter::getQFaceGenArraySize()];
  
  converter::FFace_optimised2generic(FL_,FL);
  converter::FFace_optimised2generic(FR_,FR);
  converter::QFace_optimised2generic(QL_,QL);
  converter::QFace_optimised2generic(QR_,QR);
#else
  double* FL=FL_;
  double* FR=FR_;
  const double* QL=QL_;
  const double* QR=QR_;
#endif
  
  constexpr int numberOfVariables  = MyElasticWaveSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyElasticWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyElasticWaveSolver::Order+1;
  constexpr int order              = basisSize - 1;

  kernels::idx3 idx_QLR(basisSize,basisSize,numberOfData);
  kernels::idx3 idx_FLR(basisSize,basisSize,numberOfVariables);

  double n[3]={0,0,0};
  n[normalNonZeroIndex]=1;

  double n_p[3]={0,0,0};
  double n_m[3]={0,0,0};

  double m_p[3]={0,0,0};
  double m_m[3]={0,0,0};

  double l_p[3]={0,0,0};  
  double l_m[3]={0,0,0};

  double norm_p_qr = 1.0;
  double norm_m_qr = 1.0;
  
  double FLn, FLm, FLl, FRn,FRm,FRl;
  double FL_n,FL_m,FL_l,FR_n,FR_m,FR_l;
  double FLx,FLy,FLz,FRx,FRy,FRz;
  double FL_x,FL_y,FL_z,FR_x,FR_y,FR_z;

   for (int k = 0; k < 3; k++){

     n_m[k] = n[k];
     n_p[k] = n[k];
   }

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      double rho_p=QR[idx_QLR(i,j,9)];
      double c_p_p=QR[idx_QLR(i,j,10)];
      double c_s_p=QR[idx_QLR(i,j,11)];
      double mu_p=c_s_p*c_s_p*rho_p;
      double lam_p = rho_p*c_p_p*c_p_p-2*mu_p;      

      double rho_m=QL[idx_QLR(i,j,9)];
      double c_p_m=QL[idx_QLR(i,j,10)];
      double c_s_m=QL[idx_QLR(i,j,11)];
      double mu_m=c_s_m*c_s_m*rho_m;
      double lam_m = rho_m*c_p_m*c_p_m-2*mu_m;      
  

      double Tx_m,Ty_m,Tz_m,Tx_p,Ty_p,Tz_p;
      double vx_m,vy_m,vz_m,vx_p,vy_p,vz_p;
      
      extract_tractions_and_particle_velocity(n_p,QR+idx_QLR(i,j,0),Tx_p,Ty_p,Tz_p,vx_p,vy_p,vz_p );
      extract_tractions_and_particle_velocity(n_m,QL+idx_QLR(i,j,0),Tx_m,Ty_m,Tz_m,vx_m,vy_m,vz_m ); 
      
      localBasis(n_p, m_p, l_p, 3);
      localBasis(n_m, m_m, l_m, 3);

      double Tn_m,Tm_m,Tl_m,vn_m,vm_m,vl_m;
      double Tn_p,Tm_p,Tl_p,vn_p,vm_p,vl_p;

      // rotate fields into l, m, n basis
      rotate_into_orthogonal_basis(n_m,m_m,l_m,Tx_m,Ty_m,Tz_m,Tn_m,Tm_m,Tl_m);
      rotate_into_orthogonal_basis(n_m,m_m,l_m,vx_m,vy_m,vz_m,vn_m,vm_m,vl_m);
      rotate_into_orthogonal_basis(n_p,m_p,l_p,Tx_p,Ty_p,Tz_p,Tn_p,Tm_p,Tl_p);
      rotate_into_orthogonal_basis(n_p,m_p,l_p,vx_p,vy_p,vz_p,vn_p,vm_p,vl_p);      
      
  
      // extract local s-wave and p-wave impedances
      double zs_p=rho_p*c_s_p;
      double zp_p=rho_p*c_p_p;      
      double zs_m=rho_m*c_s_m;
      double zp_m=rho_m*c_p_m;
      
      // impedance must be greater than zero !
      if (zp_p <= 0.0 || zp_m <= 0.0){
	std::cout<<zs_p<<" "<<zs_m<<" "<<zp_p<<" "<<zp_m<<"\n";
	std::cout<<" Impedance must be greater than zero ! "<< std::endl;
	std::exit(-1);
      }

      // generate interface data preserving the amplitude of the outgoing charactertritics
      // and satisfying interface conditions exactly.
      double vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p;        
      double vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m;    

      if (isBoundaryFace) {
	// 1 absorbing 0 free surface
	double r= faceIndex==2 ? 1 : 0;
	//	double r=1.;
	riemannSolver_boundary(faceIndex,r,vn_m,vm_m,vl_m,Tn_m,Tm_m,Tl_m,zp_m,zs_m,vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m);
	riemannSolver_boundary(faceIndex,r,vn_p,vm_p,vl_p,Tn_p,Tm_p,Tl_p,zp_p,zs_p,vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p);      
      }else {
	riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
	riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
	riemannSolver_Nodal(vl_p,vl_m, Tl_p, Tl_m, zs_p , zs_m, vl_hat_p , vl_hat_m, Tl_hat_p, Tl_hat_m);
      }

      //generate fluctuations in the local basis coordinates: n, m, l
      generate_fluctuations_left(zp_m,Tn_m,Tn_hat_m,vn_m,vn_hat_m,FLn);
      generate_fluctuations_left(zs_m,Tm_m,Tm_hat_m,vm_m,vm_hat_m,FLm);
      generate_fluctuations_left(zs_m,Tl_m,Tl_hat_m,vl_m,vl_hat_m,FLl);

      generate_fluctuations_right(zp_p,Tn_p,Tn_hat_p,vn_p,vn_hat_p,FRn);
      generate_fluctuations_right(zs_p,Tm_p,Tm_hat_p,vm_p,vm_hat_p,FRm);
      generate_fluctuations_right(zs_p,Tl_p,Tl_hat_p,vl_p,vl_hat_p,FRl);

      FL_n = FLn/zp_m;
      if(zs_m > 0){
	FL_m = FLm/zs_m;
	FL_l = FLl/zs_m;
      }else{
	FL_m=0;
	FL_l=0;
      }
    
      FR_n = FRn/zp_p;
      if(zs_p > 0){    
	FR_m = FRm/zs_p;
	FR_l = FRl/zs_p;
      }else{
	FR_m=0;
	FR_l=0;
      }
    
      // rotate back to the physical coordinates x, y, z
      rotate_into_physical_basis(n_m,m_m,l_m,FLn,FLm,FLl,FLx,FLy,FLz);
      rotate_into_physical_basis(n_p,m_p,l_p,FRn,FRm,FRl,FRx,FRy,FRz);
      rotate_into_physical_basis(n_m,m_m,l_m,FL_n,FL_m,FL_l,FL_x,FL_y,FL_z);
      rotate_into_physical_basis(n_p,m_p,l_p,FR_n,FR_m,FR_l,FR_x,FR_y,FR_z);
     
      // construct flux fluctuation vectors obeying the eigen structure of the PDE
      // and choose physically motivated penalties such that we can prove
      // numerical stability.

      FR[idx_FLR(i,j, 0)] = norm_p_qr/rho_p*FRx;
      FL[idx_FLR(i,j, 0)] = norm_m_qr/rho_m*FLx;
    
      FR[idx_FLR(i,j, 1)] = norm_p_qr/rho_p*FRy;
      FL[idx_FLR(i,j, 1)] = norm_m_qr/rho_m*FLy;

      FR[idx_FLR(i,j, 2)] = norm_p_qr/rho_p*FRz;
      FL[idx_FLR(i,j, 2)] = norm_m_qr/rho_m*FLz;

      FL[idx_FLR(i,j, 3)] = norm_m_qr*((2*mu_m+lam_m)*n_m[0]*FL_x+lam_m*n_m[1]*FL_y+lam_m*n_m[2]*FL_z);
      FL[idx_FLR(i,j, 4)] = norm_m_qr*((2*mu_m+lam_m)*n_m[1]*FL_y+lam_m*n_m[0]*FL_x+lam_m*n_m[2]*FL_z);
      FL[idx_FLR(i,j, 5)] = norm_m_qr*((2*mu_m+lam_m)*n_m[2]*FL_z+lam_m*n_m[0]*FL_x+lam_m*n_m[1]*FL_y);

      FR[idx_FLR(i,j, 3)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[0]*FR_x+lam_p*n_p[1]*FR_y+lam_p*n_p[2]*FR_z);
      FR[idx_FLR(i,j, 4)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[1]*FR_y+lam_p*n_p[0]*FR_x+lam_p*n_p[2]*FR_z);
      FR[idx_FLR(i,j, 5)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[2]*FR_z+lam_p*n_p[0]*FR_x+lam_p*n_p[1]*FR_y);
    
      FL[idx_FLR(i,j, 6)] =  norm_m_qr*mu_m*(n_m[1]*FL_x + n_m[0]*FL_y);
      FL[idx_FLR(i,j, 7)] =  norm_m_qr*mu_m*(n_m[2]*FL_x + n_m[0]*FL_z);
      FL[idx_FLR(i,j, 8)] =  norm_m_qr*mu_m*(n_m[2]*FL_y + n_m[1]*FL_z);

      FR[idx_FLR(i,j, 6)] = -norm_p_qr*mu_p*(n_p[1]*FR_x + n_p[0]*FR_y);
      FR[idx_FLR(i,j, 7)] = -norm_p_qr*mu_p*(n_p[2]*FR_x + n_p[0]*FR_z);
      FR[idx_FLR(i,j, 8)] = -norm_p_qr*mu_p*(n_p[2]*FR_y + n_p[1]*FR_z);
    }    
  }

#ifdef OPT_KERNELS
  converter::FFace_generic2optimised(FL,FL_);
  converter::FFace_generic2optimised(FR,FR_);
#endif 

}

//Gram Schmidt orthonormalization
void Elastic::MyElasticWaveSolver::Gram_Schmidt(double* const y, double* const z){
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
  }
  
  double norm_z = std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  
  for (int i = 0; i< 3; i++){
    z[i] =  z[i]/norm_z;
  }
}

void Elastic::MyElasticWaveSolver::localBasis(double* const n, double* const m, double* const l, int d){
  if (d == 2){
      l[0] = 0.;
      l[1] = 0.;
      l[2] = 1.0;
      
      m[0] = n[1]*l[2]-n[2]*l[1];
      m[1] = -(n[0]*l[2]-n[2]*l[0]);
      m[2] = n[0]*l[1]-n[1]*l[0];
  }else if (d == 3){
      double tol, diff_norm1, diff_norm2;
      tol = 1e-12;
      m[0] = 0.;
      m[1] = 1.;
      m[2] = 0.;
      
      diff_norm1 = pow2(n[0]-m[0]) + pow2(n[1]-m[1]) + pow2(n[2]-m[2]);
      diff_norm2 = pow2(n[0]+m[0]) + pow2(n[1]+m[1]) + pow2(n[2]+m[2]);
      
      if (diff_norm1 >= pow2(tol) && diff_norm2 >= pow2(tol)){
      	Gram_Schmidt(n, m);
      }else{
      	  m[0] = 0.;
      	  m[1] = 0.;
      	  m[2] = 1.;
      	  Gram_Schmidt(n, m);
      }
      l[0] = n[1]*m[2]-n[2]*m[1];
      l[1] = -(n[0]*m[2]-n[2]*m[0]);
      l[2] = n[0]*m[1]-n[1]*m[0];
  }
} 



void Elastic::MyElasticWaveSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
  double p=0;
  double q=0;
  double phi=0;
  double v_hat=0;
  double eta=0;

  p=z_m*v_p + sigma_p;
  q=z_p*v_m - sigma_m;

  if(z_p > 0 && z_m > 0){
    eta=(z_p*z_m)/(z_p+z_m);

    phi= eta*(p/z_p - q/z_m);
     
    sigma_hat_p=phi;
    sigma_hat_m=phi;

    v_hat_p=(q+phi)/z_m;     
    v_hat_m=(p-phi)/z_p;
  }else if(z_p > 0){
    sigma_hat_p=0;
    sigma_hat_m=sigma_m;

    v_hat_p=v_p;     
    v_hat_m=v_m;
  }else if(z_m > 0){
    sigma_hat_p=sigma_p;
    sigma_hat_m=0;

    v_hat_p=v_p;     
    v_hat_m=v_m;
  }else{
    sigma_hat_p=sigma_p;
    sigma_hat_m=sigma_m;
     
    v_hat_p=v_p;
    v_hat_m=v_m;     
  }
 }

void Elastic::MyElasticWaveSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
   double p = 0.5*(z*v + sigma);
   if(z > 0){
     v_hat = (1+r)/z*p;
     sigma_hat = (1-r)*p;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }
}

void Elastic::MyElasticWaveSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
   double q = 0.5*(z*v - sigma);
   if(z > 0){
     v_hat = (1+r)/z*q;
     sigma_hat = -(1-r)*q;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }
}



void Elastic::MyElasticWaveSolver::extract_tractions_and_particle_velocity(double* const n,const double* const Q, double& Tx,double& Ty,double& Tz,double& vx,double& vy,double& vz ){
  double sigma_xx = Q[3];
  double sigma_yy = Q[4];
  double sigma_zz = Q[5];
  double sigma_xy = Q[6];
  double sigma_xz = Q[7];
  double sigma_yz = Q[8];
  
  Tx = n[0]*sigma_xx + n[1]*sigma_xy + n[2]*sigma_xz;
  Ty = n[0]*sigma_xy + n[1]*sigma_yy + n[2]*sigma_yz;
  Tz = n[0]*sigma_xz + n[1]*sigma_yz + n[2]*sigma_zz;    
  
  vx = Q[0];
  vy = Q[1];
  vz = Q[2];    
}

void Elastic::MyElasticWaveSolver::rotate_into_orthogonal_basis(double* const n,double* const m,double* const l, double Tx,double Ty,double Tz, double& Tn, double& Tm, double& Tl){
    Tn= Tx*n[0] + Ty*n[1] + Tz*n[2];
    Tm= Tx*m[0] + Ty*m[1] + Tz*m[2];
    Tl= Tx*l[0] + Ty*l[1] + Tz*l[2];
}

void Elastic::MyElasticWaveSolver::rotate_into_physical_basis(double* const n,double* const m,double* const l, double Fn,double Fm,double Fl, double& Fx, double& Fy, double& Fz){
  Fx = n[0]*Fn + m[0]*Fm + l[0]*Fl;
  Fy = n[1]*Fn + m[1]*Fm + l[1]*Fl;
  Fz = n[2]*Fn + m[2]*Fm + l[2]*Fl;
}

void Elastic::MyElasticWaveSolver::generate_fluctuations_left(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) + (T-T_hat));
}

void Elastic::MyElasticWaveSolver::generate_fluctuations_right(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) - (T-T_hat));
}

void Elastic::MyElasticWaveSolver::riemannSolver_boundary(int faceIndex,double r, double vn , double vm , double vl, double Tn , double Tm ,double Tl , double zp, double zs , double& vn_hat , double& vm_hat ,double& vl_hat , double& Tn_hat , double& Tm_hat ,double& Tl_hat)
{
  if (faceIndex % 2  == 0) {
    riemannSolver_BC0(vn, Tn, zp, r, vn_hat, Tn_hat);
    riemannSolver_BC0(vm, Tm, zs, r, vm_hat, Tm_hat);
    riemannSolver_BC0(vl, Tl, zs, r, vl_hat, Tl_hat);	
  }
      
  if (faceIndex % 2 == 1) {
    riemannSolver_BCn(vn, Tn, zp, r, vn_hat, Tn_hat);
    riemannSolver_BCn(vm, Tm, zs, r, vm_hat, Tm_hat);
    riemannSolver_BCn(vl, Tl, zs, r, vl_hat, Tl_hat);	
  }
}


