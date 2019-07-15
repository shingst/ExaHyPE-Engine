#include "ElasticWaveEquation.h"

#include "ElasticWaveEquation_Variables.h"

#include "../../../ExaHyPE/kernels/KernelUtils.h"

#ifdef OPT_KERNELS
#include kernels/ElasticWaveEquation/converter.h
#endif  

tarch::logging::Log ElasticWaveEquation3D::ElasticWaveEquation::_log( "ElasticWaveEquation3D::ElasticWaveEquation" );

ElasticWaveEquation3D::ElasticWaveEquation::ElasticWaveEquation(double maximumMeshSize,int maximumAdaptiveMeshDepth,int DMPObservables,exahype::solvers::Solver::TimeStepping timeStepping,const std::vector<std::string>& cmdlineargs):
 AbstractElasticWaveEquation(maximumMeshSize,maximumAdaptiveMeshDepth,DMPObservables,timeStepping)
,crt(ElasticWaveEquation::Order+1,1.0/9.0)
{
  init(cmdlineargs);
}


void ElasticWaveEquation3D::ElasticWaveEquation::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

//   return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PatchWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
// }

void ElasticWaveEquation3D::ElasticWaveEquation::adjustSolution(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt){

  if(t != 0.0) return;
  
  constexpr int num_nodes = ElasticWaveEquation::Order+1;
  constexpr int numberOfData = ElasticWaveEquation::NumberOfVariables + ElasticWaveEquation::NumberOfParameters;
  kernels::idx4 id_4(num_nodes,num_nodes,num_nodes,numberOfData);
  kernels::idx3 id_3(num_nodes,num_nodes,num_nodes);

  double gl_vals_x[num_nodes*num_nodes*num_nodes];
  double gl_vals_y[num_nodes*num_nodes*num_nodes];
  double gl_vals_z[num_nodes*num_nodes*num_nodes];

  double jacobian[num_nodes*num_nodes*num_nodes];
  double q_x[num_nodes*num_nodes*num_nodes];
  double q_y[num_nodes*num_nodes*num_nodes];
  double q_z[num_nodes*num_nodes*num_nodes];
  
  double r_x[num_nodes*num_nodes*num_nodes];
  double r_y[num_nodes*num_nodes*num_nodes];
  double r_z[num_nodes*num_nodes*num_nodes];

  double s_x[num_nodes*num_nodes*num_nodes];
  double s_y[num_nodes*num_nodes*num_nodes];
  double s_z[num_nodes*num_nodes*num_nodes];  

  crt.genCoordinates(center,dx,
		     gl_vals_x,gl_vals_y,gl_vals_z,
		     jacobian,
		     q_x,q_y,q_z,
		     r_x,r_y,r_z,
		     s_x,s_y,s_z);

  int n = crt.getBlock(center,dx);
  
  for (int k=0; k< num_nodes; k++){
    for (int j=0; j< num_nodes; j++){
      for (int i=0; i< num_nodes; i++){
	double x= gl_vals_x[id_3(k,j,i)];
	double y= gl_vals_y[id_3(k,j,i)];
	double z= gl_vals_z[id_3(k,j,i)];

	if(n == 0){  
	  luh[id_4(k,j,i,0)]  = std::exp(-10*((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	  luh[id_4(k,j,i,1)]  = std::exp(-10*((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	  luh[id_4(k,j,i,2)]  = std::exp(-10*((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	}else{
	  luh[id_4(k,j,i,0)]  = 0;
	  luh[id_4(k,j,i,1)]  = 0;
	  luh[id_4(k,j,i,2)]  = 0;
	}

	// //Velocity
	// luh[id_4(k,j,i,0)]  = 0;
	// luh[id_4(k,j,i,1)]  = 0;
	// luh[id_4(k,j,i,2)]  = 0;

	// stress
	luh[id_4(k,j,i,3)]  = 0;
	luh[id_4(k,j,i,4)]  = 0;
	luh[id_4(k,j,i,5)]  = 0;  
	// luh[id_4(k,j,i,3)]  = std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	// luh[id_4(k,j,i,4)]  = luh[id_4(k,j,i,3)];
	// luh[id_4(k,j,i,5)]  = luh[id_4(k,j,i,3)];  

	luh[id_4(k,j,i,6)]  = 0;
	luh[id_4(k,j,i,7)]  = 0;
	luh[id_4(k,j,i,8)]  = 0;
  
	
	if(n == 0){
	  // elastic solid:
	  luh[id_4(k,j,i,9)]   = 2.7;   //rho
	  luh[id_4(k,j,i,10)]  = 3.343; //c(1)
	  luh[id_4(k,j,i,11)]  = 6.0; //c(2)
	}else{
	  // water column
	  // luh[id_4(k,j,i,9)]   = 1.0;   //rho
	  // luh[id_4(k,j,i,10)]  = 0.0; //c(1)
	  // luh[id_4(k,j,i,11)]  = 1.484; //c(2)
	  
	  // luh[id_4(k,j,i,9)]   = 2.7;   //rho
	  // luh[id_4(k,j,i,10)]  = 3.343; //c(1)
	  // luh[id_4(k,j,i,11)]  = 6.0; //c(2)
	  
	  luh[id_4(k,j,i,9)]   = 1.0;   //rho
	  luh[id_4(k,j,i,10)]  = 0.0; //c(1)
	  luh[id_4(k,j,i,11)]  = 1.484; //c(2)
	}
	
	
	
	luh[id_4(k,j,i,12)]  = jacobian[id_3(k,j,i)];
	
	luh[id_4(k,j,i,13)]  = q_x[id_3(k,j,i)];
	luh[id_4(k,j,i,14)]  = q_y[id_3(k,j,i)];
	luh[id_4(k,j,i,15)]  = q_z[id_3(k,j,i)];
	
	luh[id_4(k,j,i,16)] = r_x[id_3(k,j,i)];
	luh[id_4(k,j,i,17)] = r_y[id_3(k,j,i)];
	luh[id_4(k,j,i,18)] = r_z[id_3(k,j,i)];
	
	luh[id_4(k,j,i,19)] = s_x[id_3(k,j,i)];
	luh[id_4(k,j,i,20)] = s_y[id_3(k,j,i)];
	luh[id_4(k,j,i,21)] = s_z[id_3(k,j,i)];
	
	luh[id_4(k,j,i,22)] = gl_vals_x[id_3(k,j,i)];
	luh[id_4(k,j,i,23)] = gl_vals_y[id_3(k,j,i)];
	luh[id_4(k,j,i,24)] = gl_vals_z[id_3(k,j,i)];
	
      }
    }
  }
}

void ElasticWaveEquation3D::ElasticWaveEquation::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required

  double c_s=Q[10];  
  double c_p=Q[11];

  double q_x=Q[13];
  double q_y=Q[14];
  double q_z=Q[15];

  double r_x=Q[16];
  double r_y=Q[17];
  double r_z=Q[18];

  double s_x=Q[19];
  double s_y=Q[20];
  double s_z=Q[21];
  
  lambda[ 0] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*c_p;
  lambda[ 1] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*c_p;
  lambda[ 2] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*c_p;
  lambda[ 3] = 0;
  lambda[ 4] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*c_s;
  lambda[ 5] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*c_s;
  lambda[ 6] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*c_s;
  lambda[ 7] = 0;
  lambda[ 8] = 0;  
  
}


void ElasticWaveEquation3D::ElasticWaveEquation::flux(const double* const Q,double** F) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required
  
  double u=Q[0];
  double v=Q[1];
  double w=Q[2];
  
  double sigma_xx=Q[3];
  double sigma_yy=Q[4];
  double sigma_zz=Q[5];  
  double sigma_xy=Q[6];
  double sigma_xz=Q[7];
  double sigma_yz=Q[8];    

  
  double jacobian=Q[12];
  
  double q_x=Q[13];
  double q_y=Q[14];
  double q_z=Q[15];
  
  double r_x=Q[16];    
  double r_y=Q[17];
  double r_z=Q[18];

  double s_x=Q[19];    
  double s_y=Q[20];
  double s_z=Q[21];          

  
  F[0][ 0] = -jacobian*(q_x*sigma_xx+q_y*sigma_xy+q_z*sigma_xz);
  F[0][ 1] = -jacobian*(q_x*sigma_xy+q_y*sigma_yy+q_z*sigma_yz);
  F[0][ 2] = -jacobian*(q_x*sigma_xz+q_y*sigma_yz+q_z*sigma_zz);
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;
  F[0][ 6] = 0.0;
  F[0][ 7] = 0.0;
  F[0][ 8] = 0.0;
  
  F[1][ 0] = -jacobian*(r_x*sigma_xx+r_y*sigma_xy+r_z*sigma_xz);
  F[1][ 1] = -jacobian*(r_x*sigma_xy+r_y*sigma_yy+r_z*sigma_yz);
  F[1][ 2] = -jacobian*(r_x*sigma_xz+r_y*sigma_yz+r_z*sigma_zz);
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
  F[1][ 5] = 0.0;
  F[1][ 6] = 0.0;
  F[1][ 7] = 0.0;
  F[1][ 8] = 0.0;

  F[2][ 0] = -jacobian*(s_x*sigma_xx+s_y*sigma_xy+s_z*sigma_xz);
  F[2][ 1] = -jacobian*(s_x*sigma_xy+s_y*sigma_yy+s_z*sigma_yz);
  F[2][ 2] = -jacobian*(s_x*sigma_xz+s_y*sigma_yz+s_z*sigma_zz);
  F[2][ 3] = 0.0;
  F[2][ 4] = 0.0;
  F[2][ 5] = 0.0;
  F[2][ 6] = 0.0;
  F[2][ 7] = 0.0;
  F[2][ 8] = 0.0;
}


void ElasticWaveEquation3D::ElasticWaveEquation::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {


  double u_q = gradQ[0];
  double v_q = gradQ[1];
  double w_q = gradQ[2];

  double u_r = gradQ[9];
  double v_r = gradQ[10];
  double w_r = gradQ[11];

  double u_s = gradQ[18];
  double v_s = gradQ[19];
  double w_s = gradQ[20];


  double jacobian=Q[12];
  
  double q_x=Q[13];
  double q_y=Q[14];
  double q_z=Q[15];
  
  double r_x=Q[16];    
  double r_y=Q[17];
  double r_z=Q[18];

  double s_x=Q[19];    
  double s_y=Q[20];
  double s_z=Q[21];          

  

  BgradQ[0] = 0;
  BgradQ[1] = 0;
  BgradQ[2] = 0;  
  BgradQ[3] = -q_x*u_q;
  BgradQ[4] = -q_y*v_q;
  BgradQ[5] = -q_z*w_q;
  BgradQ[6] = -(q_y*u_q+q_x*v_q); //sigma_xy
  BgradQ[7] = -(q_z*u_q+q_x*w_q); //sigma_xz
  BgradQ[8] = -(q_z*v_q+q_y*w_q); //sigma_yz

  BgradQ[9] = 0;
  BgradQ[10] = 0;
  BgradQ[11] = 0;  
  BgradQ[12] = -r_x*u_r;
  BgradQ[13] = -r_y*v_r;
  BgradQ[14] = -r_z*w_r;
  BgradQ[15] = -(r_y*u_r+r_x*v_r); //sigma_xy
  BgradQ[16] = -(r_z*u_r+r_x*w_r); //sigma_xz
  BgradQ[17] = -(r_z*v_r+r_y*w_r); //sigma_yz

  BgradQ[18] = 0;
  BgradQ[19] = 0;
  BgradQ[20] = 0;  
  BgradQ[21] = -s_x*u_s;
  BgradQ[22] = -s_y*v_s;
  BgradQ[23] = -s_z*w_s;
  BgradQ[24] = -(s_y*u_s+s_x*v_s); //sigma_xy
  BgradQ[25] = -(s_z*u_s+s_x*w_s); //sigma_xz
  BgradQ[26] = -(s_z*v_s+s_y*w_s); //sigma_yz    

}



void ElasticWaveEquation3D::ElasticWaveEquation::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters

  // @todo Please implement/augment if required
  
  for(int i = 0 ; i< 25; i++){
    stateOut[ i] = stateIn[i];
  }

  for(int i = 0 ; i< 9; i++){
    fluxOut[ i] = fluxIn[ i];
  }

  

}


void ElasticWaveEquation3D::ElasticWaveEquation::coefficientMatrix(const double* const Q,const int d,double* Bn){
  std::cout << "Error: Coefficient Matrix should not be used" <<std::endl;
  std::exit(-1);
  double n[3]={0.0,0.0,0.0};
  n[d]=1.0;

  double B1[4][4];
  double B2[4][4];
  double B3[4][4];

  kernels::idx2 idx_Bn(4,4);

  for (int i =0; i< 4 ; i++){
    for (int j =0; j< 4 ; j++){
      B1[i][j]=0;
      B2[i][j]=0;
      B3[i][j]=0;  
    }
  }

  B1[0][1] = -1.0;
  B1[1][0] = -1.0;

  B2[0][2] = -1.0;
  B2[2][0] = -1.0;

  B3[0][3] = -1.0;
  B3[3][0] = -1.0;  

  for (int i =0; i< 4 ; i++){
    for (int j =0; j< 4 ; j++){
      Bn[idx_Bn(i,j)] = n[0]*B1[i][j]+n[1]*B2[i][j]+n[2]*B3[i][j];
    }
  }
}



void ElasticWaveEquation3D::ElasticWaveEquation::algebraicSource(const double* const Q,double* S) {
  // S[0] = 10*Q[0];
  // S[1] = 10*Q[1];
  // S[2] = 10*Q[2];
  // S[3] = 10*Q[3];
  // S[4] = 10*Q[0];
  // S[5] = 10*Q[1];
  // S[6] = 10*Q[2];
  // S[7] = 10*Q[3];
  // S[8] = 10*Q[2];
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = 0.0;
  S[6] = 0.0;
  S[7] = 0.0;
  S[8] = 0.0;
}


void ElasticWaveEquation3D::ElasticWaveEquation::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0, int n){

  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.1;
  double f = 0.0;
  double M0 = 1.0;

  if(n == 0){
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    
    x0[0] = 0.25;
    x0[1] = 0.5;
    x0[2] = 0.5;
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 1.*f;
    forceVector[4] = 1.*f;
    forceVector[5] = 1.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
    
    
  }else if(n == 1){
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    
    x0[0] = 0.75;
    x0[1] = 0.5;
    x0[2] = 0.5;
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 1.*f;
    forceVector[4] = 1.*f;
    forceVector[5] = 1.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
  }else if(n == 2){
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    
    x0[0] = 0.5;
    x0[1] = 0.25;
    x0[2] = 0.5;

    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 1.*f;
    forceVector[4] = 1.*f;
    forceVector[5] = 1.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;

  }else if(n == 3){
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    
    x0[0] = 0.5;
    x0[1] = 0.75;
    x0[2] = 0.5;
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 1.*f;
    forceVector[4] = 1.*f;
    forceVector[5] = 1.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
  }
    

}
void ElasticWaveEquation3D::ElasticWaveEquation::multiplyMaterialParameterMatrix(const double* const Q, double* rhs){

  double rho = Q[9];  
  double c_s   = Q[10];
  double c_p   = Q[11];

   double mu     = rho*c_s*c_s;
   double lambda = rho*c_p*c_p-2*mu;

  double jacobian = Q[12];  
 

  double rho_jacobian_inv=1.0/(rho*jacobian);

  rhs[0]=rho_jacobian_inv * rhs[0];
  rhs[1]=rho_jacobian_inv * rhs[1];
  rhs[2]=rho_jacobian_inv * rhs[2];
  double lam_temp = lambda * (rhs[3] + rhs[4] + rhs[5]);

  rhs[3]=(2*mu) * rhs[3] +lam_temp;
  rhs[4]=(2*mu) * rhs[4] +lam_temp;
  rhs[5]=(2*mu) * rhs[5] +lam_temp;

  rhs[6]= mu*rhs[6];
  rhs[7]= mu*rhs[7];
  rhs[8]= mu*rhs[8];


  rhs[9] =rho_jacobian_inv * rhs[9];
  rhs[10]=rho_jacobian_inv * rhs[10];
  rhs[11]=rho_jacobian_inv * rhs[11];

  
  lam_temp = lambda * (rhs[12] + rhs[13] + rhs[14]);

  rhs[12]=(2*mu) * rhs[12] +lam_temp;
  rhs[13]=(2*mu) * rhs[13] +lam_temp;
  rhs[14]=(2*mu) * rhs[14] +lam_temp;

  rhs[15]= mu*rhs[15];
  rhs[16]= mu*rhs[16];
  rhs[17]= mu*rhs[17];

  rhs[18] = rho_jacobian_inv * rhs[18];
  rhs[19] = rho_jacobian_inv * rhs[19];
  rhs[20] = rho_jacobian_inv * rhs[20];
  
  lam_temp = lambda * (rhs[21] + rhs[22] + rhs[23]);

  rhs[21]=(2*mu) * rhs[21] +lam_temp;
  rhs[22]=(2*mu) * rhs[22] +lam_temp;
  rhs[23]=(2*mu) * rhs[23] +lam_temp;

  rhs[24]= mu*rhs[24];
  rhs[25]= mu*rhs[25];
  rhs[26]= mu*rhs[26];  

}


// void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex,bool isBoundaryFace, int faceIndex){

//   constexpr int numberOfVariables  = ElasticWaveEquation::NumberOfVariables;
//   constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
//   constexpr int numberOfParameters = ElasticWaveEquation::NumberOfParameters;
//   constexpr int numberOfData       = numberOfVariables+numberOfParameters;
//   constexpr int basisSize          = ElasticWaveEquation::Order+1;
//   constexpr int order              = basisSize - 1; 

//   //changed for dynamic rupture
//   kernels::idx3 idx_QLR(basisSize,basisSize,numberOfData);
//   kernels::idx3 idx_FLR(basisSize,basisSize,numberOfVariables);


//   double n[3]={0,0,0};

//   n[normalNonZeroIndex]=1;

//   double n_p[3]={0,0,0};
//   double n_m[3]={0,0,0};

//   double m_p[3]={0,0,0};
//   double m_m[3]={0,0,0};

//   double l_p[3]={0,0,0};  
//   double l_m[3]={0,0,0};

//   //std::cout<<isBoundaryFace<<std::endl;


//   double norm_p_qr;
//   double norm_m_qr;
  
//   double FLn, FLm, FLl; 
//   double FRn,FRm,FRl;
//   double FL_n,FL_m,FL_l;
//   double FR_n,FR_m,FR_l;
//   double FLx,FLy,FLz ;
//   double FRx,FRy,FRz ;
//   double FL_x,FL_y,FL_z;
//   double FR_x,FR_y,FR_z;
  


//   //std::cout<<isBoundaryFace<<std::endl;
//   for (int i = 0; i < basisSize; i++) {
//     for (int j = 0; j < basisSize; j++) {
//       //implemeneted for dynamic rupture


//       double rho_p=QR[idx_QLR(i,j,9)];
//       double c_s_p=QR[idx_QLR(i,j,10)];
//       double c_p_p=QR[idx_QLR(i,j,11)];
      

//       double mu_p=c_s_p*c_s_p*rho_p;
//       double lam_p = rho_p*c_p_p*c_p_p-2*mu_p;      


//       double rho_m=QL[idx_QLR(i,j,9)];
      
//       double c_s_m=QL[idx_QLR(i,j,10)];
//       double c_p_m=QL[idx_QLR(i,j,11)];

//       double mu_m=c_s_m*c_s_m*rho_m;
//       double lam_m = rho_m*c_p_m*c_p_m-2*mu_m;      

//       get_normals(normalNonZeroIndex, norm_p_qr, n_p, QR + idx_QLR(i,j,0));
//       get_normals(normalNonZeroIndex, norm_m_qr, n_m, QL + idx_QLR(i,j,0));    
      

//       double Tx_m,Ty_m,Tz_m;
//       double Tx_p,Ty_p,Tz_p;
      
//       double vx_m,vy_m,vz_m;
//       double vx_p,vy_p,vz_p;
      
//       extract_tractions_and_particle_velocity(n_p,QR+idx_QLR(i,j,0),Tx_p,Ty_p,Tz_p,vx_p,vy_p,vz_p );
//       extract_tractions_and_particle_velocity(n_m,QL+idx_QLR(i,j,0),Tx_m,Ty_m,Tz_m,vx_m,vy_m,vz_m ); 
      
//       localBasis(n_p, m_p, l_p, 3);
//       localBasis(n_m, m_m, l_m, 3);

//       double Tn_m,Tm_m,Tl_m;
//       double vn_m,vm_m,vl_m;
//       double Tn_p,Tm_p,Tl_p;
//       double vn_p,vm_p,vl_p;

//       // rotate fields into l, m, n basis
//       rotate_into_orthogonal_basis(n_m,m_m,l_m,Tx_m,Ty_m,Tz_m,Tn_m,Tm_m,Tl_m);
//       rotate_into_orthogonal_basis(n_m,m_m,l_m,vx_m,vy_m,vz_m,vn_m,vm_m,vl_m);
//       rotate_into_orthogonal_basis(n_p,m_p,l_p,Tx_p,Ty_p,Tz_p,Tn_p,Tm_p,Tl_p);
//       rotate_into_orthogonal_basis(n_p,m_p,l_p,vx_p,vy_p,vz_p,vn_p,vm_p,vl_p);      
      
  
//     // extract local s-wave and p-wave impedances
//       double zs_p=rho_p*c_s_p;
//       double zs_m=rho_m*c_s_m;
      
//       double zp_p=rho_p*c_p_p;
//       double zp_m=rho_m*c_p_m;
      
//     // impedance must be greater than zero !
//     if (zp_p <= 0.0 || zp_m <= 0.0){
//       std::cout<<zs_p<<" "<<zs_m<<" "<<zp_p<<" "<<zp_m<<"\n";
//       std::cout<<" Impedance must be greater than zero ! "<< std::endl;
//       std::exit(-1);
//     }

//     // generate interface data preserving the amplitude of the outgoing charactertritics
//     // and satisfying interface conditions exactly.
    
//     double vn_hat_p,vm_hat_p,vl_hat_p;    
//     double vn_hat_m,vm_hat_m,vl_hat_m;    
//     double Tn_hat_p,Tm_hat_p,Tl_hat_p;    
//     double Tn_hat_m,Tm_hat_m,Tl_hat_m;    


//     if (isBoundaryFace) {
//       double r= faceIndex==1 ? 1 : 0;
//       riemannSolver_boundary(faceIndex,r,vn_m,vm_m,vl_m,Tn_m,Tm_m,Tl_m,zp_m,zs_m,vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m);
//       riemannSolver_boundary(faceIndex,r,vn_p,vm_p,vl_p,Tn_p,Tm_p,Tl_p,zp_p,zs_p,vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p);      

//     }else {
//       riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
//       riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
//       riemannSolver_Nodal(vl_p,vl_m, Tl_p, Tl_m, zs_p , zs_m, vl_hat_p , vl_hat_m, Tl_hat_p, Tl_hat_m);

//     }

//     // generate fluctuations in the local basis coordinates: n, m
//     generate_fluctuations_left(zp_m,Tn_m,Tn_hat_m,vn_m,vn_hat_m,FLn);
//     generate_fluctuations_left(zs_m,Tm_m,Tm_hat_m,vm_m,vm_hat_m,FLm);
//     generate_fluctuations_left(zs_m,Tl_m,Tl_hat_m,vl_m,vl_hat_m,FLl);

//     generate_fluctuations_right(zp_p,Tn_p,Tn_hat_p,vn_p,vn_hat_p,FRn);
//     generate_fluctuations_right(zs_p,Tm_p,Tm_hat_p,vm_p,vm_hat_p,FRm);
//     generate_fluctuations_right(zs_p,Tl_p,Tl_hat_p,vl_p,vl_hat_p,FRl);


    
//     FL_n = FLn/zp_m;
//     FL_m=0;
//     FL_l=0;
    
//     if(zs_m > 0){
//     FL_m = FLm/zs_m;
//     FL_l = FLl/zs_m;
//      }

    
//     FR_n = FRn/zp_p;
//     FR_m=0;
//     FR_l=0;
//     if(zs_p > 0){    
//     FR_m = FRm/zs_p;
//     FR_l = FRl/zs_p;
//     }

    
//     // rotate back to the physical coordinates x, y, z
//     rotate_into_physical_basis(n_m,m_m,l_m,FLn,FLm,FLl,FLx,FLy,FLz);
//     rotate_into_physical_basis(n_p,m_p,l_p,FRn,FRm,FRl,FRx,FRy,FRz);
//     rotate_into_physical_basis(n_m,m_m,l_m,FL_n,FL_m,FL_l,FL_x,FL_y,FL_z);
//     rotate_into_physical_basis(n_p,m_p,l_p,FR_n,FR_m,FR_l,FR_x,FR_y,FR_z);
     
//     // construct flux fluctuation vectors obeying the eigen structure of the PDE
//     // and choose physically motivated penalties such that we can prove
//     // numerical stability.

//     FR[idx_FLR(i,j, 0)] = norm_p_qr/rho_p*FRx;
//     FL[idx_FLR(i,j, 0)] = norm_m_qr/rho_m*FLx;
    
//     FR[idx_FLR(i,j, 1)] = norm_p_qr/rho_p*FRy;
//     FL[idx_FLR(i,j, 1)] = norm_m_qr/rho_m*FLy;

//     FR[idx_FLR(i,j, 2)] = norm_p_qr/rho_p*FRz;
//     FL[idx_FLR(i,j, 2)] = norm_m_qr/rho_m*FLz;
    

//     FL[idx_FLR(i,j, 3)] = norm_m_qr*((2*mu_m+lam_m)*n_m[0]*FL_x+lam_m*n_m[1]*FL_y+lam_m*n_m[2]*FL_z);
//     FL[idx_FLR(i,j, 4)] = norm_m_qr*((2*mu_m+lam_m)*n_m[1]*FL_y+lam_m*n_m[0]*FL_x+lam_m*n_m[2]*FL_z);
//     FL[idx_FLR(i,j, 5)] = norm_m_qr*((2*mu_m+lam_m)*n_m[2]*FL_z+lam_m*n_m[0]*FL_x+lam_m*n_m[1]*FL_y);

//     FR[idx_FLR(i,j, 3)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[0]*FR_x+lam_p*n_p[1]*FR_y+lam_p*n_p[2]*FR_z);
//     FR[idx_FLR(i,j, 4)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[1]*FR_y+lam_p*n_p[0]*FR_x+lam_p*n_p[2]*FR_z);
//     FR[idx_FLR(i,j, 5)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[2]*FR_z+lam_p*n_p[0]*FR_x+lam_p*n_p[1]*FR_y);
    
//     FL[idx_FLR(i,j, 6)] =  norm_m_qr*mu_m*(n_m[1]*FL_x + n_m[0]*FL_y);
//     FL[idx_FLR(i,j, 7)] =  norm_m_qr*mu_m*(n_m[2]*FL_x + n_m[0]*FL_z);
//     FL[idx_FLR(i,j, 8)] =  norm_m_qr*mu_m*(n_m[2]*FL_y + n_m[1]*FL_z);

//     FR[idx_FLR(i,j, 6)] = -norm_p_qr*mu_p*(n_p[1]*FR_x + n_p[0]*FR_y);
//     FR[idx_FLR(i,j, 7)] = -norm_p_qr*mu_p*(n_p[2]*FR_x + n_p[0]*FR_z);
//     FR[idx_FLR(i,j, 8)] = -norm_p_qr*mu_p*(n_p[2]*FR_y + n_p[1]*FR_z);
    
//     // double x = QR[idx_QLR(i,j,22)];
//     // // std::cout<< pow(x-0.5,2)<< std::endl;
    
//     // if (pow(x-0.5,2) < 1e-5) {
//     //   std::cout<<"x: "<< x<< n_m[0]-n_p[0] << " " << n_m[1]-n_p[1] << " " << n_m[2]-n_p[2] << std::endl;
//     //   }
//     }    
//   }

// }

void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver(double* FL_,double* FR_,const double* const QL_,const double* const QR_,const double dt,const int normalNonZeroIndex,bool isBoundaryFace, int faceIndex){


#ifdef OPT_KERNELS
  double* FL = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getFFaceGenArraySize()];
  double* FR = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getFFaceGenArraySize()];
  const double* QL = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getQFaceGenArraySize()];
  const double* QR = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getQFaceGenArraySize()];

  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_optimised2generic(FL_,FL);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_optimised2generic(FR_,FR);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::QFace_optimised2generic(QL_,QL);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::QFace_optimised2generic(QR_,QR);
#else
  double* FL=FL_;
  double* FR=FR_;
  const double* QL=QL_;
  const double* QR=QR_;
#endif  

  constexpr int numberOfVariables  = ElasticWaveEquation::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = ElasticWaveEquation::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = ElasticWaveEquation::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx3 idx_QLR(basisSize, basisSize,numberOfData);

  kernels::idx3 idx_FLR(basisSize, basisSize,NumberOfVariables);

  double n[3]={0,0,0};
  n[normalNonZeroIndex]=1;
  
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {    
      double qm_x=QL[idx_QLR(i,j,13)];
      double qm_y=QL[idx_QLR(i,j,14)];
      double qm_z=QL[idx_QLR(i,j,15)];    
      double rm_x=QL[idx_QLR(i,j,16)];
      double rm_y=QL[idx_QLR(i,j,17)];
      double rm_z=QL[idx_QLR(i,j,18)];
      double sm_x=QL[idx_QLR(i,j,19)];
      double sm_y=QL[idx_QLR(i,j,20)];
      double sm_z=QL[idx_QLR(i,j,21)];
      
      double qp_x=QR[idx_QLR(i,j,13)];
      double qp_y=QR[idx_QLR(i,j,14)];
      double qp_z=QR[idx_QLR(i,j,15)];    
      double rp_x=QR[idx_QLR(i,j,16)];
      double rp_y=QR[idx_QLR(i,j,17)];
      double rp_z=QR[idx_QLR(i,j,18)];
      double sp_x=QR[idx_QLR(i,j,19)];
      double sp_y=QR[idx_QLR(i,j,20)];
      double sp_z=QR[idx_QLR(i,j,21)];
      
      
      double n_p[3]={0,0,0};
      double n_m[3]={0,0,0};
      
      double m_p[3]={0,0,0};
      double m_m[3]={0,0,0};
      
      double l_p[3]={0,0,0};
      double l_m[3]={0,0,0};
      
      
      double norm_p_qr;
      double norm_m_qr;
      
      if (normalNonZeroIndex == 0){
  	norm_m_qr = std::sqrt(qm_x*qm_x + qm_y*qm_y + qm_z*qm_z);
	n_m[0] = qm_x/norm_m_qr;
	n_m[1] = qm_y/norm_m_qr;
	n_m[2] = qm_z/norm_m_qr;  
	
	norm_p_qr = std::sqrt(qp_x*qp_x + qp_y*qp_y  + qp_z*qp_z );
	n_p[0] = qp_x/norm_p_qr;
	n_p[1] = qp_y/norm_p_qr;
	n_p[2] = qp_z/norm_p_qr;  
      }
      
      if (normalNonZeroIndex == 1){
  
  norm_m_qr = std::sqrt(rm_x*rm_x + rm_y*rm_y + rm_z*rm_z);
  n_m[0] = rm_x/norm_m_qr;
  n_m[1] = rm_y/norm_m_qr;
  n_m[2] = rm_z/norm_m_qr;  
  
  norm_p_qr = std::sqrt(rp_x*rp_x + rp_y*rp_y  + rp_z*rp_z);
  n_p[0] = rp_x/norm_p_qr;
  n_p[1] = rp_y/norm_p_qr;
  n_p[2] = rp_z/norm_p_qr;  
      }
      
      if (normalNonZeroIndex == 2){
  
  norm_m_qr = std::sqrt(sm_x*sm_x + sm_y*sm_y + sm_z*sm_z);
  n_m[0] = sm_x/norm_m_qr;
  n_m[1] = sm_y/norm_m_qr;
  n_m[2] = sm_z/norm_m_qr;  
  
  norm_p_qr = std::sqrt(sp_x*sp_x + sp_y*sp_y  + sp_z*sp_z);
  n_p[0] = sp_x/norm_p_qr;
  n_p[1] = sp_y/norm_p_qr;
  n_p[2] = sp_z/norm_p_qr;  
      }
      

      
      double rho_m  = QL[idx_QLR(i,j,9)];   // km/s
      double cs_m   = QL[idx_QLR(i,j,10)];   // km/s
      double cp_m   = QL[idx_QLR(i,j,11)];   // km/s
      
      double rho_p  = QR[idx_QLR(i,j,9)];   // km/s
      double cs_p   = QR[idx_QLR(i,j,10)];   // km/s
      double cp_p   = QR[idx_QLR(i,j,11)];   // km/s
      
      
      
      double mu_m = rho_m*cs_m*cs_m;
      double lam_m = rho_m*cp_m*cp_m - 2*mu_m;
      
      double mu_p = rho_p*cs_p*cs_p;
      double lam_p = rho_p*cp_p*cp_p - 2*mu_p;
      
      
      // extract tractions and particle velocities
      double sigma_m_xx =  QL[idx_QLR(i,j,3)];
      double sigma_m_yy =  QL[idx_QLR(i,j,4)];
      double sigma_m_zz =  QL[idx_QLR(i,j,5)];
      double sigma_m_xy =  QL[idx_QLR(i,j,6)];
      double sigma_m_xz =  QL[idx_QLR(i,j,7)];
      double sigma_m_yz =  QL[idx_QLR(i,j,8)];
      
      double sigma_p_xx =  QR[idx_QLR(i,j,3)];
      double sigma_p_yy =  QR[idx_QLR(i,j,4)];
      double sigma_p_zz =  QR[idx_QLR(i,j,5)];
      double sigma_p_xy =  QR[idx_QLR(i,j,6)];
      double sigma_p_xz =  QR[idx_QLR(i,j,7)];
      double sigma_p_yz =  QR[idx_QLR(i,j,8)];    
      
      
      double Tx_m = n_m[0]*sigma_m_xx + n_m[1]*sigma_m_xy + n_m[2]*sigma_m_xz;
      double Ty_m = n_m[0]*sigma_m_xy + n_m[1]*sigma_m_yy + n_m[2]*sigma_m_yz;
      double Tz_m = n_m[0]*sigma_m_xz + n_m[1]*sigma_m_yz + n_m[2]*sigma_m_zz;
      
      double Tx_p = n_p[0]*sigma_p_xx + n_p[1]*sigma_p_xy + n_p[2]*sigma_p_xz;
      double Ty_p = n_p[0]*sigma_p_xy + n_p[1]*sigma_p_yy + n_p[2]*sigma_p_yz;
      double Tz_p = n_p[0]*sigma_p_xz + n_p[1]*sigma_p_yz + n_p[2]*sigma_p_zz;
      
      
      double vx_m = QL[idx_QLR(i,j,0)];
      double vy_m = QL[idx_QLR(i,j,1)];
      double vz_m = QL[idx_QLR(i,j,2)];    

      double vx_p = QR[idx_QLR(i,j,0)];
      double vy_p = QR[idx_QLR(i,j,1)];
      double vz_p = QR[idx_QLR(i,j,2)];    
      
      localBasis(n_m, m_m, l_m,3);
      localBasis(n_p, m_p, l_p,3);
      
      // rotate tractions and particle velocities into orthogonal coordinates: n, m
      double Tn_m= Tx_m*n_m[0] + Ty_m*n_m[1] + Tz_m*n_m[2];
      double Tm_m= Tx_m*m_m[0] + Ty_m*m_m[1] + Tz_m*m_m[2];
      double Tl_m= Tx_m*l_m[0] + Ty_m*l_m[1] + Tz_m*l_m[2];
      
      double Tn_p= Tx_p*n_p[0] + Ty_p*n_p[1] + Tz_p*n_p[2];
      double Tm_p= Tx_p*m_p[0] + Ty_p*m_p[1] + Tz_p*m_p[2];
      double Tl_p= Tx_p*l_p[0] + Ty_p*l_p[1] + Tz_p*l_p[2];
      
      double vn_m= vx_m*n_m[0] + vy_m*n_m[1] + vz_m*n_m[2];
      double vm_m= vx_m*m_m[0] + vy_m*m_m[1] + vz_m*m_m[2];
      double vl_m= vx_m*l_m[0] + vy_m*l_m[1] + vz_m*l_m[2];
      
      double vn_p= vx_p*n_p[0] + vy_p*n_p[1] + vz_p*n_p[2];
      double vm_p= vx_p*m_p[0] + vy_p*m_p[1] + vz_p*m_p[2];
      double vl_p= vx_p*l_p[0] + vy_p*l_p[1] + vz_p*l_p[2];    
      
      // extract local s-wave and p-wave impedances
      double zs_p=rho_p*cs_p;
      double zs_m=rho_m*cs_m;
      
      double zp_p=rho_p*cp_p;
      double zp_m=rho_m*cp_m;
      
      // // impedance must be greater than zero !
      // if (zs_p <= 0.0 || zs_m <= 0.0 || zp_p <= 0.0 || zp_m <= 0.0){
      //   std::cout<<zs_p<<' '<<zs_m<<' '<<zp_p<<' '<<zp_m<<'\n';
      //   std::cout<<' Impedance must be greater than zero ! '<<'\n';
      //   std::exit(-1);
      // }
      
      
      // generate interface data preserving the amplitude of the outgoing charactertritics
      // and satisfying interface conditions exactly.
      double vn_hat_p=0;
      double vm_hat_p=0;
      double vl_hat_p=0;    
      
      double vn_hat_m=0;
      double vm_hat_m=0;
      double vl_hat_m=0;    
      
      double Tn_hat_p=0;
      double Tm_hat_p=0;
      double Tl_hat_p=0;    
      
      double Tn_hat_m=0;
      double Tm_hat_m=0;
      double Tl_hat_m=0;
      
      // data is generated by solving a local Riemann problem and contraining the solutions against
      // physical interface conditions
      
      
      if (isBoundaryFace) {
  if (faceIndex == 0) {
    
    double r = 0.;
    
    riemannSolver_BC0(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
    riemannSolver_BC0(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
    riemannSolver_BC0(vl_p, Tl_p, zs_p, r, vl_hat_p, Tl_hat_p);
    
    
    riemannSolver_BC0(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
    riemannSolver_BC0(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
    riemannSolver_BC0(vl_m, Tl_m, zs_m, r, vl_hat_m, Tl_hat_m);  
  }
  
  
  if (faceIndex == 1) {
    double r = 0.;
    
    riemannSolver_BCn(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
    riemannSolver_BCn(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
    riemannSolver_BCn(vl_p, Tl_p, zs_p, r, vl_hat_p, Tl_hat_p);
    
    riemannSolver_BCn(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
    riemannSolver_BCn(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
    riemannSolver_BCn(vl_m, Tl_m, zs_m, r, vl_hat_m, Tl_hat_m);  
  }
  
  
  if (faceIndex == 2) {
    double r = 1.;
    
    riemannSolver_BC0(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
    riemannSolver_BC0(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
    riemannSolver_BC0(vl_p, Tl_p, zs_p, r, vl_hat_p, Tl_hat_p);  
    
    riemannSolver_BC0(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
    riemannSolver_BC0(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
    riemannSolver_BC0(vl_m, Tl_m, zs_m, r, vl_hat_m, Tl_hat_m);  
    
    
  }
  
  if (faceIndex == 3) {
    
    double r = 0.;
    
    riemannSolver_BCn(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
    riemannSolver_BCn(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
    riemannSolver_BCn(vl_p, Tl_p, zs_p, r, vl_hat_p, Tl_hat_p);  
    
    riemannSolver_BCn(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
    riemannSolver_BCn(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
    riemannSolver_BCn(vl_m, Tl_m, zs_m, r, vl_hat_m, Tl_hat_m);  
    
  }
  
        
  if (faceIndex == 4) {
    double r = 0.;
    
    riemannSolver_BC0(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
    riemannSolver_BC0(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
    riemannSolver_BC0(vl_p, Tl_p, zs_p, r, vl_hat_p, Tl_hat_p);  
    
    riemannSolver_BC0(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
    riemannSolver_BC0(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
    riemannSolver_BC0(vl_m, Tl_m, zs_m, r, vl_hat_m, Tl_hat_m);  
    
    
  }
  
  if (faceIndex == 5) {
    
    double r = 0.;
    
    riemannSolver_BCn(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
    riemannSolver_BCn(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
    riemannSolver_BCn(vl_p, Tl_p, zs_p, r, vl_hat_p, Tl_hat_p);  
    
    riemannSolver_BCn(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
    riemannSolver_BCn(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
    riemannSolver_BCn(vl_m, Tl_m, zs_m, r, vl_hat_m, Tl_hat_m);  
    
  }
  
      }
      else {
  
  riemannSolver_Nodal(vn_p,vn_m, Tn_p,Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
  riemannSolver_Nodal(vm_p,vm_m, Tm_p,Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
  riemannSolver_Nodal(vl_p,vl_m, Tl_p,Tl_m, zs_p , zs_m, vl_hat_p , vl_hat_m, Tl_hat_p, Tl_hat_m);
      }
      
      //std::cout<< n[0]- n_p[0] << "  " <<  n[1]- n_p[1] << std::endl;
      // std::cout<< n[0]- n_m[0] << "  " <<  n[1]- n_m[1] << std::endl;
      
      // std::cout << std::endl;
      
      // generate fluctuations in the local basis coordinates: n, m
      double FLn = 0.5*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
      double FLm = 0.5*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));
      double FLl = 0.5*(zs_m*(vl_m-vl_hat_m) + (Tl_m-Tl_hat_m));    
      
      double FRn = 0.5*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
      double FRm = 0.5*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));
      double FRl = 0.5*(zs_p*(vl_p-vl_hat_p) - (Tl_p-Tl_hat_p));    
      
      
      double FL_n = 0.5/zp_m*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
      double FL_m=0;
      double FL_l=0;
      
      if(zs_m > 0){
  FL_m = 0.5/zs_m*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));
  FL_l = 0.5/zs_m*(zs_m*(vl_m-vl_hat_m) + (Tl_m-Tl_hat_m));    
      }
      
      
      double FR_n = 0.5/zp_p*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
      double FR_m=0;
    double FR_l=0;
    
    if(zs_p > 0){    
      FR_m = 0.5/zs_p*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));
      FR_l = 0.5/zs_p*(zs_p*(vl_p-vl_hat_p) - (Tl_p-Tl_hat_p));
    }

    // rotate back to the physical coordinates x, y
    double FLx = n_m[0]*FLn + m_m[0]*FLm + l_m[0]*FLl;
    double FLy = n_m[1]*FLn + m_m[1]*FLm + l_m[1]*FLl;
    double FLz = n_m[2]*FLn + m_m[2]*FLm + l_m[2]*FLl;    

    double FRx = n_p[0]*FRn + m_p[0]*FRm + l_p[0]*FRl;
    double FRy = n_p[1]*FRn + m_p[1]*FRm + l_p[1]*FRl;
    double FRz = n_p[2]*FRn + m_p[2]*FRm + l_p[2]*FRl;

    double FL_x = n_m[0]*FL_n + m_m[0]*FL_m + l_m[0]*FL_l;
    double FL_y = n_m[1]*FL_n + m_m[1]*FL_m + l_m[1]*FL_l;
    double FL_z = n_m[2]*FL_n + m_m[2]*FL_m + l_m[2]*FL_l;

    double FR_x = n_p[0]*FR_n + m_p[0]*FR_m + l_p[0]*FR_l;
    double FR_y = n_p[1]*FR_n + m_p[1]*FR_m + l_p[1]*FR_l;
    double FR_z = n_p[2]*FR_n + m_p[2]*FR_m + l_p[2]*FR_l;
     
    // construct flux fluctuation vectors obeying the eigen structure of the PDE
    // and choose physically motivated penalties such that we can prove
    // numerical stability.

    
    FR[idx_FLR(i,j, 0)] = norm_p_qr/rho_p*FRx;
    FR[idx_FLR(i,j, 1)] = norm_p_qr/rho_p*FRy;
    FR[idx_FLR(i,j, 2)] = norm_p_qr/rho_p*FRz;
    
    FL[idx_FLR(i,j, 0)] = norm_m_qr/rho_m*FLx;    
    FL[idx_FLR(i,j, 1)] = norm_m_qr/rho_m*FLy;    
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
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_generic2optimised(FL,FL_);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_generic2optimised(FR,FR_);

  delete[] FL;  
  delete[] FR;  
  delete[] QL;  
  delete[] QR;
#endif  

  
}


void ElasticWaveEquation3D::ElasticWaveEquation::Gram_Schmidt(double* y, double* z){
  //Gram Schmidt orthonormalization
 
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
  }

  double norm_z = std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  
  for (int i = 0; i< 3; i++){
    z[i] =  z[i]/norm_z;
  }

}

void ElasticWaveEquation3D::ElasticWaveEquation::localBasis(double* n, double * m, double* l, int d){

  if (d == 2)
    {
      l[0] = 0.;
      l[1] = 0.;
      l[2] = 1.0;
      
      m[0] = n[1]*l[2]-n[2]*l[1];
      m[1] = -(n[0]*l[2]-n[2]*l[0]);
      m[2] = n[0]*l[1]-n[1]*l[0];
      
    }
  
  
  if (d == 3)
    {
      double tol, diff_norm1, diff_norm2;
      
      tol = 1e-12;
      m[0] = 0.;
      m[1] = 1.;
      m[2] = 0.;
      
      diff_norm1 =  std::sqrt(pow(n[0]-m[0],2) + pow(n[1]-m[1], 2) + pow(n[2]-m[2], 2));
      
      diff_norm2 =  std::sqrt(pow(n[0]+m[0],2) + pow(n[1]+m[1], 2) + pow(n[2]+m[2], 2));

      
      if (diff_norm1 >= tol && diff_norm2 >= tol){
        Gram_Schmidt(n, m);}  else
        {
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



void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){

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


void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_Nodal_rupture(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& v_hat ,double& sigma_hat_p, double& sigma_hat_m, double alpha){

  double p=0;
  double q=0;
  double phi=0;
  double sigma_hat=0;
  //  double v_hat=0;
  double eta=0;

   p=z_m*v_p + sigma_p;
   q=z_p*v_m - sigma_m;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= eta*(p/z_p - q/z_m);

   v_hat = 1/(eta+alpha) * phi;
   sigma_hat = alpha/(eta+alpha) * phi;
   
   sigma_hat_p= sigma_hat;
   sigma_hat_m= sigma_hat;

   v_hat_m=(p-sigma_hat_m)/z_p - v_hat;
   v_hat_p=(q-sigma_hat_p)/z_m + v_hat;

 }



void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   if(z > 0){
     v_hat = (1+r)/z*p;
     sigma_hat = (1-r)*p;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }

}


void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   if(z > 0){
     v_hat = (1+r)/z*q;
     sigma_hat = -(1-r)*q;
      }else{
     v_hat = v;
     sigma_hat = sigma;
   }


}



void ElasticWaveEquation3D::ElasticWaveEquation::get_normals(int normalNonZeroIndex,double& norm, double* n,const double* Q)

{

  double q_x=Q[13];
  double q_y=Q[14];
  double q_z=Q[15];
  
  double r_x=Q[16];
  double r_y=Q[17];
  double r_z=Q[18];
  
  double s_x=Q[19];
  double s_y=Q[20];
  double s_z=Q[21];

  if (normalNonZeroIndex == 0){
  
    norm = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z);
    n[0] = q_x/norm;
    n[1] = q_y/norm;
    n[2] = q_z/norm;  

  }
      
  if (normalNonZeroIndex == 1){

    norm = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
    n[0] = r_x/norm;
    n[1] = r_y/norm;
    n[2] = r_z/norm;  
  }

  if (normalNonZeroIndex == 2){

    norm = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z);
    n[0] = s_x/norm;
    n[1] = s_y/norm;
    n[2] = s_z/norm;  
  }

}


void ElasticWaveEquation3D::ElasticWaveEquation::extract_tractions_and_particle_velocity(double* n,const double* Q, double& Tx,double& Ty,double& Tz,double& vx,double& vy,double& vz ){

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

void ElasticWaveEquation3D::ElasticWaveEquation::rotate_into_orthogonal_basis(double* n,double* m,double* l, double Tx,double Ty,double Tz, double& Tn, double& Tm, double& Tl){
    Tn= Tx*n[0] + Ty*n[1] + Tz*n[2];
    Tm= Tx*m[0] + Ty*m[1] + Tz*m[2];
    Tl= Tx*l[0] + Ty*l[1] + Tz*l[2];
}

void ElasticWaveEquation3D::ElasticWaveEquation::rotate_into_physical_basis(double* n,double* m,double* l, double Fn,double Fm,double Fl, double& Fx, double& Fy, double& Fz){

  Fx = n[0]*Fn + m[0]*Fm + l[0]*Fl;
  Fy = n[1]*Fn + m[1]*Fm + l[1]*Fl;
  Fz = n[2]*Fn + m[2]*Fm + l[2]*Fl;
  
}

void ElasticWaveEquation3D::ElasticWaveEquation::generate_fluctuations_left(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) + (T-T_hat));
}

void ElasticWaveEquation3D::ElasticWaveEquation::generate_fluctuations_right(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) - (T-T_hat));
}



void ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_boundary(int faceIndex,double r, double vn , double vm , double vl, double Tn , double Tm ,double Tl , double zp, double zs , double& vn_hat , double& vm_hat ,double& vl_hat , double& Tn_hat , double& Tm_hat ,double& Tl_hat)
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


exahype::solvers::Solver::RefinementControl ElasticWaveEquation3D::ElasticWaveEquation::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}
