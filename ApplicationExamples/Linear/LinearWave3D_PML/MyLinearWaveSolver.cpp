#include "MyLinearWaveSolver.h"

#include "MyLinearWaveSolver_Variables.h"

#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"



tarch::logging::Log Linear::MyLinearWaveSolver::_log( "Linear::MyLinearWaveSolver" );


void Linear::MyLinearWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Linear::MyLinearWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 5
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    //initialise luh
    constexpr int basisSize = MyLinearWaveSolver::Order+1;
    int numberOfData=MyLinearWaveSolver::NumberOfParameters+MyLinearWaveSolver::NumberOfVariables;
    
    kernels::idx4 id_4(basisSize,basisSize,basisSize,numberOfData);
    
    double offset_x=center[0]-0.5*dx[0];
    double offset_y=center[1]-0.5*dx[1];
    double offset_z=center[2]-0.5*dx[2];
    
    
    double width_x=dx[0];
    double width_y=dx[1];
    double width_z=dx[2];

    double dpml =  1.5;
    int n = 2;
    double tol = 1e-3; 
   
    
    double xa = 1.5;
    double xb = 8.5;
    
    double ya = 1.5;
    double yb = 8.5;
    
    double za = 1.5;
    double zb = 8.5;

    double d_x = 0.0;
    double d_y = 0.0;
    double d_z = 0.0;
      
    
     for (int k=0; k< basisSize; k++){
      for (int j=0; j< basisSize; j++){
	for (int i=0; i< basisSize; i++){
     
	  double x  =  (offset_x+width_x*kernels::gaussLegendreNodes[basisSize-1][i]);
	  double y  =  (offset_y+width_y*kernels::gaussLegendreNodes[basisSize-1][j]);
	  double z  =  (offset_z+width_z*kernels::gaussLegendreNodes[basisSize-1][k]);

	   
	  //pressure
	  luh[id_4(k,j,i,0)]  = 0*std::exp(-100*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));
	   
	  // particle velocities
	  luh[id_4(k,j,i,1)]  = 0;
	  luh[id_4(k,j,i,2)]  = 0;
	  luh[id_4(k,j,i,3)]  = 0;

	  // auxiliary variables
	  luh[id_4(k,j,i,4)]  = 0;
	  luh[id_4(k,j,i,5)]  = 0;
	   
	  // material parameters
	  luh[id_4(k,j,i,6)]  = 1.0;  // density [g/cm^3]
	  luh[id_4(k,j,i,7)]  = 1.484;   // wavespeed [km/s]

	  double c =  luh[id_4(k,j,i,7)];
	  double d0 = 0.5*(n+1)*c/(2*dpml)* log(1/tol);
	  // PML damping profiles
	  if (x < xa){
	    d_x = d0*pow((xa-x)/dpml, n);
	  }

	  if (x > xb){
	    d_x = d0*pow((x-xb)/dpml, n);
	  }

	  // if (y < ya){
	  //   d_y = d0*pow((ya-y)/dpml, n);
	  // }

	  if (y > yb){
	    d_y = d0*pow((y-yb)/dpml, n);
	  }

	  if (z < za){
	    d_z = d0*pow((za-z)/dpml, n);
	  }

	  if (z > zb){
	    d_z = d0*pow((z-zb)/dpml, n);
	  }

 	   // PML damping functions
	  luh[id_4(k,j,i,8)]   = d_x;  // density [g/cm^3]
	  luh[id_4(k,j,i,9)]   = d_y;   // wavespeed [km/s]
	  luh[id_4(k,j,i,10)]  = d_z;   // wavespeed [km/s]
	   

	   
	   
	}
      }
    }

  }
}

void Linear::MyLinearWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 5

  // @todo Please implement/augment if required
  constexpr int numberOfVariables  = MyLinearWaveSolver::NumberOfVariables;
  constexpr int numberOfParameters = MyLinearWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  

  for (int i = 0; i<numberOfData; i++){
    stateOut[i] = stateIn[i] ;
  }
 
  for (int i = 0; i< numberOfVariables; i++){
  fluxOut[i] =  fluxIn[i];
 }
}

exahype::solvers::Solver::RefinementControl Linear::MyLinearWaveSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Linear::MyLinearWaveSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 6 + 5
  
  // @todo Please implement/augment if required
  double c   =  Q[7];
  
  lambda[ 0] =  c;
  lambda[ 1] = -c;
  lambda[ 2] = 0.0;
  lambda[ 3] = 0.0;
  lambda[ 4] = 0.0;
  lambda[ 5] = 0.0;
}


void Linear::MyLinearWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 5
  
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



void  Linear::MyLinearWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  
  double p_x = gradQ[0];  
  double p_y = gradQ[6];  
  double p_z = gradQ[12];  
   
  
  BgradQ[0]= 0;  
  BgradQ[1]= -p_x;
  BgradQ[2]= 0;
  BgradQ[3]= 0;
  BgradQ[4]= 0;
  BgradQ[5]= 0;

  BgradQ[6]= 0;  
  BgradQ[7]= 0;
  BgradQ[8]= -p_y;
  BgradQ[9]= 0;
  BgradQ[10]= 0;
  BgradQ[11]= 0;

  BgradQ[12]= 0;  
  BgradQ[13]= 0;
  BgradQ[14]= 0;
  BgradQ[15]= -p_z;
  BgradQ[16]= 0;
  BgradQ[17]= 0;
}

void  Linear::MyLinearWaveSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0, int n) {
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  
  //return;
  
  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

   forceVector[4] = 0.0;
   forceVector[5] = 0.0;

  if(n == 0){
  
    x0[0] = 5.0;
    x0[1] = 5.;
    x0[2] = 5.;
    
    forceVector[0] = f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
    
  }else if(n == 1){
    
    x0[0] = 7.5;
    x0[1] = 5.;
    x0[2] = 5.;
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;

  }else if(n == 2){
  
  x0[0] = 5.;
  x0[1] = 2.5;
  x0[2] = 5.;
    
  forceVector[0] = 0.0;
  forceVector[1] = 0.0;
  forceVector[2] = 0.0;
  forceVector[3] = 0.0;

 }else if(n == 3){
  
   x0[0] = 5.;
   x0[1] = 7.5;
   x0[2] = 5.;
    
   forceVector[0] = 0.0;
   forceVector[1] = 0.0;
   forceVector[2] = 0.0;
   forceVector[3] = 0.0;

 }    
}


void Linear::MyLinearWaveSolver::algebraicSource(const double* const Q,double* const S) {

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


    /**
     * @TODO LR : document
     */
void Linear::MyLinearWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  // @todo Please implement/augment if required

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

void Linear::MyLinearWaveSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

}


void Linear::MyLinearWaveSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

}


void Linear::MyLinearWaveSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
   double p=0;
   double q=0;
   double phi=0;
   double v_hat=0;
   double eta=0;

   p = z_m*v_p + sigma_p;
   q = z_p*v_m - sigma_m;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= eta*(p/z_p - q/z_m);

   sigma_hat_p=phi;
   sigma_hat_m=phi;

   v_hat_m=(p-phi)/z_p;
   v_hat_p=(q+phi)/z_m;

 }


void Linear::MyLinearWaveSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex){

  constexpr int numberOfVariables  = MyLinearWaveSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyLinearWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyLinearWaveSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx3 idx_QLR(basisSize, basisSize, numberOfData);

  kernels::idx3 idx_FLR(basisSize, basisSize, numberOfVariables);

  double n[3]={0,0,0};
  n[normalNonZeroIndex]=1;

  double cp_m;
  double rho_m;
  double lam_m;

  double cp_p;
  double rho_p;
  double lam_p;

  double n_p[3]={0,0,0};
  double n_m[3]={0,0,0};

  double m_p[3]={0,0,0};
  double m_m[3]={0,0,0};

  double l_p[3]={0,0,0};
  double l_m[3]={0,0,0};

  double norm_p_qr;
  double norm_m_qr;

  //std::cout<<" riemann solver called"<<std::endl;
  //std::exit(-1);


  for (int k = 0; k < 3; k++){

    n_m[k] = n[k];
    n_p[k] = n[k];
    
    
    
  }

  norm_m_qr = 1.0;
  norm_p_qr = 1.0;   

  
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {

      rho_m = QL[idx_QLR(i,j,6)];
      cp_m  = QL[idx_QLR(i,j,7)];
      
      lam_m = cp_m*cp_m*rho_m;

      rho_p = QR[idx_QLR(i,j,6)];
      cp_p  = QR[idx_QLR(i,j,7)];
      
      lam_p = cp_p*cp_p*rho_p;


      double dm_x = QL[idx_QLR(i,j,8)];
      double dm_y = QL[idx_QLR(i,j,9)];
      double dm_z = QL[idx_QLR(i,j,10)];

      double dp_x = QR[idx_QLR(i,j,8)];
      double dp_y = QR[idx_QLR(i,j,9)];
      double dp_z = QR[idx_QLR(i,j,10)];

      
      
      
      double v_m=QL[idx_QLR(i,j,1)]*n_m[0]+QL[idx_QLR(i,j,2)]*n_m[1]+QL[idx_QLR(i,j,3)]*n_m[2];
      double v_p=QR[idx_QLR(i,j,1)]*n_p[0]+QR[idx_QLR(i,j,2)]*n_p[1]+QR[idx_QLR(i,j,3)]*n_p[2];;
      
      double sigma_m = QL[idx_QLR(i,j,0)];
      double sigma_p = QR[idx_QLR(i,j,0)];
      
      
      double z_p=rho_p*cp_p;
      double z_m=rho_m*cp_m;
      
      double v_hat_p=0;
      double v_hat_m=0;
      double sigma_hat_p=0;
      double sigma_hat_m=0;

      double r = 0.;

      if (faceIndex  == 2) { r = 1.0;} 
      
      
      if (isBoundaryFace) { // external boundaries
	
	 if (faceIndex % 2  == 0) {
	  
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	}
	
	
	if  (faceIndex % 2  == 1) {
	 
	  
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	  
	}
	
      }
      else {// interelment boundaries
	riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);
      }
      
      
      
      FR[idx_FLR(i, j, 0)] = -norm_p_qr*0.5*lam_p*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p);
      FL[idx_FLR(i, j, 0)] =  norm_m_qr*0.5*lam_m*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m);
      

      FR[idx_FLR(i, j, 1)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[0];
      FL[idx_FLR(i, j, 1)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[0];
      
      FR[idx_FLR(i, j, 2)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[1];
      FL[idx_FLR(i, j, 2)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[1];
      
      FR[idx_FLR(i, j, 3)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[2];
      FL[idx_FLR(i, j, 3)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[2];

      FR[idx_FLR(i, j, 4)] = -(dp_x-dp_y)*norm_p_qr*0.5*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p)*n_p[1];
      FL[idx_FLR(i, j, 4)] =  (dm_x-dm_y)*norm_m_qr*0.5*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m)*n_m[1];
      
      FR[idx_FLR(i, j, 5)] = -(dp_x-dp_z)*norm_p_qr*0.5*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p)*n_p[2];
      FL[idx_FLR(i, j, 5)] =  (dm_x-dm_z)*norm_m_qr*0.5*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m)*n_m[2];

      
    }
    
  }
}

