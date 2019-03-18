#include "MyLinearSolver.h"

#include "MyLinearSolver_Variables.h"

#include "../../ExaHyPE/kernels/KernelUtils.h"


tarch::logging::Log Linear::MyLinearSolver::_log( "Linear::MyLinearSolver" );


void Linear::MyLinearSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Linear::MyLinearSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void Linear::MyLinearSolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  double x_0[3]={2.0, 2.0, 2.0};
  double exponent=0.0;

  for(int i=0 ; i < 3 ; i++){
    exponent = exponent + (x[i]-x_0[i])*(x[i]-x_0[i]);
  }

  double x0 = x[0]-x_0[0];
  double x1 = x[1]-x_0[1];
  double x2 = x[2]-x_0[2];
  
  Q[ 0] = 0.0; //std::exp(-exponent);
  Q[ 1] = 0.0;
  Q[ 2] = 0.0;
  Q[ 3] = 0.0;
  Q[ 4] = 0.0;
  Q[ 5] = 0.0;   // Material parameters:
  Q[ 6] = 1.0;   // density (g/cm^3)
  Q[ 7] = 1.484;   // acoustic wave speed (km/s)
  Q[ 8] = 0.0;  //*(x[0]-x_0[0])*(x[0]-x_0[0]);
  Q[ 9] = 0.0;
  Q[10] = 0.0;

  double c = Q[7];

  double dpml = 0.55;
  int n = 1;
  double tol = 1e-3; 
  
  double d0 = 0.5*(n+1)*c/(2*dpml)* log(1/tol);

  double xa = 0.55;
  double xb = 4.45;

  // std::cout<<d0<<std::endl;

  // std::exit(-1);
  //double d0 = 0.0;
  
  if (x[0] < xa){
    Q[ 8] = d0*pow((xa-x[0])/dpml, n);
  }

  if (x[0] > xb){
    Q[8] = d0*pow((x[0]-xb)/dpml, n);
  }


  if (x[1] < xa){
    Q[ 9] = d0*pow((xa-x[1])/dpml, n);
  }

  if (x[1] > xb){
     Q[ 9]= d0*pow((x[1]-xb)/dpml, n);
  }

  if (x[2] < xa){
     Q[10] = d0*pow((xa-x[2])/dpml, n);
  }

  if (x[2] > xb){
    
     Q[10] = d0*pow((x[2]-xb)/dpml, n);
  }
  

  // std::cout<<x[0]<< "  " <<  Q[ 8]<<std::endl;
  // std::cout<<x[1]<< "  " <<  Q[ 9]<<std::endl;
  // std::cout<<x[2]<< "  " <<  Q[ 10]<<std::endl;

  // std::cout<<std::endl;
   
}

void Linear::MyLinearSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required
  double c   =  Q[7];
  
  lambda[ 0] =  c;
  lambda[ 1] = -c;
  lambda[ 2] = 0.0;
  lambda[ 3] = 0.0;
  lambda[ 4] = 0.0;
  lambda[ 5] = 0.0;
}


void Linear::MyLinearSolver::flux(const double* const Q,double** const F) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required
  double u = Q[1];
  double v = Q[2];
  double w = Q[3];
  
  F[0][ 0] = -u;
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;

  F[1][ 0] = -v;
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;
  F[1][ 4] = -v;
  F[1][ 5] = 0.0;

  F[2][ 0] = -w;
  F[2][ 1] = 0.0;
  F[2][ 2] = 0.0;
  F[2][ 3] = 0.0;
  F[2][ 4] = 0.0;
  F[2][ 5] = -w;
}


void Linear::MyLinearSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {

  double d_x = Q[8];
  double d_y = Q[9];
  double d_z = Q[10];
  
  double p_x = gradQ[0];  
  double u_x = gradQ[1];
  double v_x = gradQ[2];
  double w_x = gradQ[3];

  double p_y = gradQ[6];  
  double u_y = gradQ[7];
  double v_y = gradQ[8];
  double w_y = gradQ[9];

  double p_z = gradQ[12];  
  double u_z = gradQ[13];
  double v_z = gradQ[14];
  double w_z = gradQ[15];    
  
  //BgradQ[0]= -u_x;
  BgradQ[0]= 0;  
  BgradQ[1]= -p_x;
  BgradQ[2]= 0;
  BgradQ[3]= 0;
  BgradQ[4]= 0;
  BgradQ[5]= 0;

  //BgradQ[4]=-v_y;
  BgradQ[6]= 0;  
  BgradQ[7]= 0;
  BgradQ[8]= -p_y;
  BgradQ[9]= 0;
  BgradQ[10]= 0;
  BgradQ[11]= 0;

  //BgradQ[8]=-w_z;
  BgradQ[12]= 0;  
  BgradQ[13]= 0;
  BgradQ[14]= 0;
  BgradQ[15]= -p_z;
  BgradQ[16]= 0;
  BgradQ[17]= 0;

}


void Linear::MyLinearSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters

  // @todo Please implement/augment if required

  for(int i = 0 ; i< 11; i++){
    stateOut[ i] = stateIn[i];
  }


  fluxOut[ 0] = fluxIn[ 0];
  fluxOut[ 1] = fluxIn[ 1];
  fluxOut[ 2] = fluxIn[ 2];
  fluxOut[ 3] = fluxIn[ 3];
  fluxOut[ 4] = fluxIn[ 4];
  fluxOut[ 5] = fluxIn[ 5];
  
}


void Linear::MyLinearSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn){

  double d_x = Q[8];
  double d_y = Q[9];
  double d_z = Q[10];

  double n[3]={0.0,0.0,0.0};
  n[d]=1.0;

  double B1[6][6];
  double B2[6][6];
  double B3[6][6];

  kernels::idx2 idx_Bn(6,6);

  for (int i =0; i< 6 ; i++){
    for (int j =0; j< 6 ; j++){
      B1[i][j]=0;
      B2[i][j]=0;
      B3[i][j]=0;	
    }
  }

  B1[0][1] = -1.0;
  B1[1][0] = -1.0;

  B2[0][2] = -1.0;
  B2[2][0] = -1.0;
  B2[4][2] = -(d_x-d_y);

  B3[0][3] = -1.0;
  B3[3][0] = -1.0;
  B3[5][3] = -(d_x-d_z);

  for (int i =0; i< 6 ; i++){
    for (int j =0; j< 6 ; j++){
      Bn[idx_Bn(i,j)] = n[0]*B1[i][j]+n[1]*B2[i][j]+n[2]*B3[i][j];
    }
  }
}


void Linear::MyLinearSolver::algebraicSource(const double* const Q,double* const S) {

  double p = Q[0];
  double u = Q[1];
  double v = Q[2];
  double w = Q[3];

  double sigma = Q[4];
  double psi   = Q[5];


  double rho = Q[6];  
  double c   = Q[7];
  double mu  = rho*c*c;
  
  double d_x = Q[8];
  double d_y = Q[9];
  double d_z = Q[10];

  S[0] = d_x*p - mu*sigma - mu*psi;
  S[1] = d_x*u;
  S[2] = d_y*v;
  S[3] = d_z*w;
  S[4] = d_y*sigma;
  S[5] = d_z*psi;
  
}

void Linear::MyLinearSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0){

  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  double f = 0.0;
  double c = 1.484;
  double rho = 1.0;
  double M0 = 1.;

  
  
  
  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

  x0[0] = 1.5;
  x0[1] = 2.5;
  x0[2] = 2.5;
  
  forceVector[0] = 1.*f;
  forceVector[1] = 0.0;
  forceVector[2] = 0.0;
  forceVector[3] = 0.0;
  forceVector[4] = 0.0;
  forceVector[5] = 0.0;

}

void Linear::MyLinearSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs){

  double rho = Q[6];  
  double c   = Q[7];
  double mu  = rho*c*c;


  double d_x = Q[8];
  double d_y = Q[9];
  double d_z = Q[10];

  rhs[0]=   mu * rhs[0];
  rhs[1]=1/rho * rhs[1];
  rhs[2]=1/rho * rhs[2];
  rhs[3]=1/rho * rhs[3];

  rhs[4]= 0*rhs[4];
  rhs[5]= 0*rhs[5];

  rhs[6]=  mu * rhs[6];
  rhs[7]=1/rho * rhs[7];
  rhs[8]=1/rho * rhs[8];
  rhs[9]=1/rho * rhs[9];

  rhs[10]= (d_x-d_y)*rhs[10];
  rhs[11]= 0*rhs[11];

  rhs[12]=   mu * rhs[12];
  rhs[13]=1/rho * rhs[13];
  rhs[14]=1/rho * rhs[14];
  rhs[15]=1/rho * rhs[15];

  rhs[16]= 0*rhs[16];
  rhs[17]= (d_x-d_z)*rhs[17];
  
}

void Linear::MyLinearSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,double* const tempFaceUnknownsArray,double** const tempStateSizedVectors,double** const tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex,bool isBoundaryFace, int faceIndex){

  constexpr int numberOfVariables  = MyLinearSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyLinearSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyLinearSolver::Order+1;
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

  double norm_p_qr=1.0;
  double norm_m_qr=1.0;

  //std::cout<<isBoundaryFace<<std::endl;
  

  //std::cout<<isBoundaryFace<<std::endl;

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {

      double c_p=QR[idx_QLR(i,j,7)];
      double rho_p=QR[idx_QLR(i,j,6)];
      double mu_p=c_p*c_p*rho_p;

      double c_m=QL[idx_QLR(i,j,7)];
      double rho_m=QL[idx_QLR(i,j,6)];
      double mu_m=c_m*c_m*rho_m;

      

      for(int k = 0 ; k< 3 ;k++){
	n_m[k] = n[k];
	n_p[k] = n[k];
      }

      localBasis(n_p, m_p, l_p, 3);
      localBasis(n_m, m_m, l_m, 3);     
       
      
      
      double v_m=QL[idx_QLR(i,j,1)]*n_m[0]+QL[idx_QLR(i,j,2)]*n_m[1]+QL[idx_QLR(i,j,3)]*n_m[2];
      double v_p=QR[idx_QLR(i,j,1)]*n_p[0]+QR[idx_QLR(i,j,2)]*n_p[1]+QR[idx_QLR(i,j,3)]*n_p[2];

      double sigma_m = QL[idx_QLR(i,j,0)];
      double sigma_p = QR[idx_QLR(i,j,0)];

      double dm_x = QL[idx_QLR(i,j,8)];
      double dm_y = QL[idx_QLR(i,j,9)];
      double dm_z = QL[idx_QLR(i,j,10)];

      double dp_x = QR[idx_QLR(i,j,8)];
      double dp_y = QR[idx_QLR(i,j,9)];
      double dp_z = QR[idx_QLR(i,j,10)];

      double z_p=rho_p*c_p;
      double z_m=rho_m*c_m;

      double v_hat_p=0;
      double v_hat_m=0;
      double sigma_hat_p=0;
      double sigma_hat_m=0;

      if (isBoundaryFace) {
	if (faceIndex == 0) {
	
	  double r = 0.;
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	}
      
      
	if (faceIndex == 1) {
	  double r = 0.;
	
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	}
      
      
	if (faceIndex == 2) {
	
	  double r = 0.;
	
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	
	}
      
	if (faceIndex == 3) {
	
	  double r = 0.;
	
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
        
	}

	if (faceIndex == 4) {
	
	  double r = 0.;
	
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	
	}
      
	if (faceIndex == 5) {
	
	  double r = 0.;
	
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
        
	}


	
      }
      else {
	riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);      
      }
      




      FR[idx_FLR(i,j,0)] = -norm_p_qr*0.5*mu_p*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p);
      FL[idx_FLR(i,j,0)] =  norm_m_qr*0.5*mu_m*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m);

      FR[idx_FLR(i,j,1)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[0];
      FL[idx_FLR(i,j,1)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[0];

      FR[idx_FLR(i,j,2)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[1];
      FL[idx_FLR(i,j,2)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[1];

      FR[idx_FLR(i,j,3)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[2];
      FL[idx_FLR(i,j,3)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[2];

      FR[idx_FLR(i,j,4)] = -(dp_x-dp_y)*norm_p_qr*0.5*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p)*n_p[1];
      FL[idx_FLR(i,j,4)] =  (dm_x-dm_y)*norm_m_qr*0.5*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m)*n_m[1];

      FR[idx_FLR(i,j,5)] = -(dp_x-dp_z)*norm_p_qr*0.5*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p)*n_p[2];
      FL[idx_FLR(i,j,5)] =  (dm_x-dm_z)*norm_m_qr*0.5*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m)*n_m[2];

      
      

    }
  }
}



void Linear::MyLinearSolver::Gram_Schmidt(double* const y, double* const z){
  //Gram Schmidt orthonormalization
 
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
    z[i] = 1.0/std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2])*z[i];
  }

}

void Linear::MyLinearSolver::localBasis(double* const n, double* const m, double* const l, int d){

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
	Gram_Schmidt(n, m);}	else
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



void Linear::MyLinearSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
   double p=0;
   double q=0;
   double phi=0;
   double v_hat=0;
   double eta=0;

   p=z_m*v_p + sigma_p;
   q=z_p*v_m - sigma_m;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= eta*(p/z_p - q/z_m);

   sigma_hat_p=phi;
   sigma_hat_m=phi;

   v_hat_m=(p-phi)/z_p;
   v_hat_p=(q+phi)/z_m;

 }


void Linear::MyLinearSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

}


void Linear::MyLinearSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

}




exahype::solvers::Solver::RefinementControl Linear::MyLinearSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}
