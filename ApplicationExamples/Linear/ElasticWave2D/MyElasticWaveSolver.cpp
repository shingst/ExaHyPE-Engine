#include "MyElasticWaveSolver.h"

#include "MyElasticWaveSolver_Variables.h"

#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"


tarch::logging::Log Elastic::MyElasticWaveSolver::_log( "Elastic::MyElasticWaveSolver" );


void Elastic::MyElasticWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Elastic::MyElasticWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 3
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    //initialise luh

    constexpr int basisSize = MyElasticWaveSolver::Order+1;
    int numberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

     kernels::idx3 id_xyf(basisSize,basisSize,numberOfData);

     double offset_x=center[0]-0.5*dx[0];
     double offset_y=center[1]-0.5*dx[1];
    

     double width_x=dx[0];
     double width_y=dx[1];
   

    for (int j=0; j< basisSize; j++){
      for (int i=0; i< basisSize; i++){
     
	double x  =  (offset_x+width_x*kernels::legendre::nodes[basisSize-1][i]);
	double y  =  (offset_y+width_y*kernels::legendre::nodes[basisSize-1][j]);


	//particle velocities
	luh[id_xyf(j,i,0)]  = std::exp(-100*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
	luh[id_xyf(j,i,1)]  = std::exp(-100*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));

	// stress vector
	luh[id_xyf(j,i,2)]  = 0;
	luh[id_xyf(j,i,3)]  = 0;
	luh[id_xyf(j,i,4)]  = 0;

	// parameters
	luh[id_xyf(j,i,5)]  = 2.7;     // density [g/cm^3]
	luh[id_xyf(j,i,6)]  = 6.0;     // pwavespeed [km/s]
	luh[id_xyf(j,i,7)]  = 3.343;   // swavespeed [km/s]


      }
    }
  }
}

void Elastic::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 3

  // @todo Please implement/augment if required
  constexpr int numberOfVariables  = MyElasticWaveSolver::NumberOfVariables;
  constexpr int numberOfParameters = MyElasticWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
 

  for (int i = 0; i<numberOfData; i++){
    stateOut[i] = stateIn[i];
  }
 
  for (int i = 0; i< numberOfVariables; i++){
  fluxOut[i] =  fluxIn[i];
 }
  
}

exahype::solvers::Solver::RefinementControl Elastic::MyElasticWaveSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Elastic::MyElasticWaveSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 3
  
  // @todo Please implement/augment if required
  double cp = Q[6];
  double cs = Q[7];
   
  lambda[0] = cp;
  lambda[1] = cs;
  lambda[2] = -cp;
  lambda[3] = -cs;
  lambda[4] = 0.0;
}


void Elastic::MyElasticWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 3
  
  // @todo Please implement/augment if required

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



void  Elastic::MyElasticWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required

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

// void Elastic::MyElasticWaveSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn){
//   static tarch::logging::Log _log("MyElasticWaveSolver::coefficientMatrix");

//   const double rho  = Q[5];   // km/s
//   const double cp   = Q[6];   // km/s
//   const double cs   = Q[7];   // km/s
    
//   double mu = rho*cs*cs;
//   double lam = rho*cp*cp - 2*mu;
  
//   double nv[2] = {0.0, 0.0};

//   nv[d] = 1.0;
  
//   double B1[5][5];
//   double B2[5][5];
   
//   for(int i=0; i<5; i++) {
//     for(int j=0; j<5; j++) {
//       B1[i][j] = 0.0;
//       B2[i][j] = 0.0;
//     }
//   }
  
//   B1[0][2] = -1/rho; 
//   B1[1][4] = -1/rho;
//   B1[2][0] = -(2*mu + lam); 
//   B1[4][1] = -mu;

//   B2[0][4] = -1/rho; 
//   B2[1][3] = -1/rho;
//   B2[3][1] = -(2*mu + lam); 
//   B2[4][0] = -mu;
  
//   for(int i=0; i<5; i++) {
//     for(int j=0; j<5; j++) {
      
//       Bn[i*5 + j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
//     }
//   }

// }

void  Elastic::MyElasticWaveSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0, int n) {
  // @todo Please implement/augment if required

  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

  if(n == 0){
  
    x0[0] = 5.0;
    x0[1] = 5.0;
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = f;
    forceVector[3] = f;
    forceVector[4] = 0.0;
    
  }else if(n == 1){
  
    x0[0] = 5.0;
    x0[1] = 2.5;
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.*f;;
    forceVector[3] = 0.*f;;
    forceVector[4] = 0.0;

  }    
}

    /**
     * @TODO LR : document
     */
void Elastic::MyElasticWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  // @todo Please implement/augment if required

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





void Elastic::MyElasticWaveSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex) {
  
  constexpr int numberOfVariables  = MyElasticWaveSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyElasticWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyElasticWaveSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx2 idx_QLR(basisSize, numberOfData);
  
  kernels::idx2 idx_FLR(basisSize, numberOfVariables);
  
  double n[2]={0,0};
  n[normalNonZeroIndex]=1;
  
  double rho_p;
  double cp_p;
  double cs_p;
  
  double mu_p;
  double lam_p;

  double rho_m;
  double cp_m;
  double cs_m;
  
  double mu_m;
  double lam_m;
  
  double n_p[2]={0,0};
  double n_m[2]={0,0};
  
  double m_p[2]={0,0};
  double m_m[2]={0,0};

  double l_p[2]={0,0};
  double l_m[2]={0,0};
  
  double norm_p_qr;
  double norm_m_qr;
  
  for (int k = 0; k < 2; k++){
    
    n_m[k] = n[k];
    n_p[k] = n[k];
    
  }
  
  norm_m_qr = 1.0;
  norm_p_qr = 1.0;   
  
  
  for (int i = 0; i < basisSize; i++) {

    // extract parameters
    rho_m = QL[idx_QLR(i,5)];
    cp_m  = QL[idx_QLR(i,6)];
    cs_m  = QL[idx_QLR(i,7)];
    
    mu_m = cs_m*cs_m*rho_m;
    lam_m = cp_m*cp_m*rho_m-2*mu_m;

    // extract parameters
    rho_p = QR[idx_QLR(i,5)];
    cp_p  = QR[idx_QLR(i,6)];
    cs_p  = QR[idx_QLR(i,7)];
    
    mu_p = cs_p*cs_p*rho_p;
    lam_p = cp_p*cp_p*rho_p-2*mu_p;
    

    //get_normals(normalNonZeroIndex, norm_p_qr, n_p, QR + idx_QLR(i,j,0));
    //get_normals(normalNonZeroIndex, norm_m_qr, n_m, QL + idx_QLR(i,j,0));    
    
    double Tx_m,Ty_m;
    double Tx_p,Ty_p;
    
    double vx_m,vy_m;
    double vx_p,vy_p;

    extract_tractions_and_particle_velocity(n_p,QR+idx_QLR(i,0),Tx_p,Ty_p,vx_p,vy_p);
    extract_tractions_and_particle_velocity(n_m,QL+idx_QLR(i,0),Tx_m,Ty_m,vx_m,vy_m); 
    
    localBasis(n_p, m_p);
    localBasis(n_m, m_m);
    
    double Tn_m,Tm_m;
    double vn_m,vm_m;
    double Tn_p,Tm_p;
    double vn_p,vm_p;
    
    // rotate fields into l, m, n basis
    rotate_into_orthogonal_basis(n_m,m_m,Tx_m,Ty_m,Tn_m,Tm_m);
    rotate_into_orthogonal_basis(n_m,m_m,vx_m,vy_m,vn_m,vm_m);
    rotate_into_orthogonal_basis(n_p,m_p,Tx_p,Ty_p,Tn_p,Tm_p);
    rotate_into_orthogonal_basis(n_p,m_p,vx_p,vy_p,vn_p,vm_p);      
    
    
    // extract local s-wave and p-wave impedances
    double zs_p=rho_p*cs_p;
    double zs_m=rho_m*cs_m;
    
    double zp_p=rho_p*cp_p;
    double zp_m=rho_m*cp_m;
    
    

    // generate interface data preserving the amplitude of the outgoing charactertritics
    // and satisfying interface conditions exactly.
    
    double vn_hat_p,vm_hat_p;    
    double vn_hat_m,vm_hat_m;    
    double Tn_hat_p,Tm_hat_p;    
    double Tn_hat_m,Tm_hat_m;

    double FLn, FLm; 
    double FRn,FRm;
    double FL_n,FL_m;
    double FR_n,FR_m;
    double FLx,FLy;
    double FRx,FRy;
    double FL_x,FL_y;
    double FR_x,FR_y;
  
    
    
    if (isBoundaryFace) {
      
      double r= faceIndex==2 ? 1 : 0;
      riemannSolver_boundary(faceIndex,r,vn_m,vm_m,Tn_m,Tm_m,zp_m,zs_m,vn_hat_m,vm_hat_m,Tn_hat_m,Tm_hat_m);
      riemannSolver_boundary(faceIndex,r,vn_p,vm_p,Tn_p,Tm_p,zp_p,zs_p,vn_hat_p,vm_hat_p,Tn_hat_p,Tm_hat_p);
      
    }else{
      riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
      riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
    }
    
    
    // generate fluctuations in the local basis coordinates: n, m
    generate_fluctuations_left(zp_m,Tn_m,Tn_hat_m,vn_m,vn_hat_m,FLn);
    generate_fluctuations_left(zs_m,Tm_m,Tm_hat_m,vm_m,vm_hat_m,FLm);
    
    generate_fluctuations_right(zp_p,Tn_p,Tn_hat_p,vn_p,vn_hat_p,FRn);
    generate_fluctuations_right(zs_p,Tm_p,Tm_hat_p,vm_p,vm_hat_p,FRm);
    
    
    FL_n = FLn/zp_m;
    FL_m=0;
    if(zs_m > 0){
      FL_m = FLm/zs_m;
    }
    
    FR_n = FRn/zp_p;
    FR_m=0;
    if(zs_p > 0){    
      FR_m = FRm/zs_p;
    }
    
    // rotate back to the physical coordinates x, y, z
    rotate_into_physical_basis(n_m,m_m,FLn,FLm,FLx,FLy);
    rotate_into_physical_basis(n_p,m_p,FRn,FRm,FRx,FRy);
    rotate_into_physical_basis(n_m,m_m,FL_n,FL_m,FL_x,FL_y);
    rotate_into_physical_basis(n_p,m_p,FR_n,FR_m,FR_x,FR_y);
     
    // construct flux fluctuation vectors obeying the eigen structure of the PDE
    // and choose physically motivated penalties such that we can prove
    // numerical stability.
    
    FR[idx_FLR(i, 0)] = norm_p_qr/rho_p*FRx;
    FL[idx_FLR(i, 0)] = norm_m_qr/rho_m*FLx;
    
    FR[idx_FLR(i, 1)] = norm_p_qr/rho_p*FRy;
    FL[idx_FLR(i, 1)] = norm_m_qr/rho_m*FLy;
    
    
    FL[idx_FLR(i, 2)] = norm_m_qr*((2*mu_m+lam_m)*n_m[0]*FL_x+lam_m*n_m[1]*FL_y);
    FL[idx_FLR(i, 3)] = norm_m_qr*((2*mu_m+lam_m)*n_m[1]*FL_y+lam_m*n_m[0]*FL_x);
    
    FR[idx_FLR(i, 2)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[0]*FR_x+lam_p*n_p[1]*FR_y);
    FR[idx_FLR(i, 3)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[1]*FR_y+lam_p*n_p[0]*FR_x);
    
    
    FL[idx_FLR(i, 4)] =  norm_m_qr*mu_m*(n_m[1]*FL_x + n_m[0]*FL_y);
    FR[idx_FLR(i, 4)] = -norm_p_qr*mu_p*(n_p[1]*FR_x + n_p[0]*FR_y);
  }
}


void Elastic::MyElasticWaveSolver::extract_tractions_and_particle_velocity(double* const n,const double* const Q, double& Tx,double& Ty,double& vx,double& vy){

  double sigma_xx = Q[2];
  double sigma_yy = Q[3];
  double sigma_xy = Q[4];
  
  
  Tx = n[0]*sigma_xx + n[1]*sigma_xy;
  Ty = n[0]*sigma_xy + n[1]*sigma_yy;   
  
  vx = Q[0];
  vy = Q[1];   
}

void Elastic::MyElasticWaveSolver::rotate_into_orthogonal_basis(double* const n,double* const m, double Tx,double Ty, double& Tn, double& Tm){
    Tn= Tx*n[0] + Ty*n[1];
    Tm= Tx*m[0] + Ty*m[1];
}

void Elastic::MyElasticWaveSolver::rotate_into_physical_basis(double* const n,double* const m, double Fn,double Fm, double& Fx, double& Fy){

  Fx = n[0]*Fn + m[0]*Fm;
  Fy = n[1]*Fn + m[1]*Fm;
  
}

void Elastic::MyElasticWaveSolver::generate_fluctuations_left(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) + (T-T_hat));
}

void Elastic::MyElasticWaveSolver::generate_fluctuations_right(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) - (T-T_hat));
}

void Elastic::MyElasticWaveSolver::riemannSolver_boundary(int faceIndex,double r, double vn , double vm, double Tn , double Tm, double zp, double zs , double& vn_hat , double& vm_hat, double& Tn_hat , double& Tm_hat)
{

  if (faceIndex % 2  == 0) {

    riemannSolver_BC0(vn, Tn, zp, r, vn_hat, Tn_hat);
    riemannSolver_BC0(vm, Tm, zs, r, vm_hat, Tm_hat);
  }
      
      
  if (faceIndex % 2 == 1) {

    riemannSolver_BCn(vn, Tn, zp, r, vn_hat, Tn_hat);
    riemannSolver_BCn(vm, Tm, zs, r, vm_hat, Tm_hat);	
  }

}

void Elastic::MyElasticWaveSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

}


void Elastic::MyElasticWaveSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

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


void Elastic::MyElasticWaveSolver::localBasis(double* const n, double* const m){
      m[0] = n[1];
      m[1] = -n[0];
}
