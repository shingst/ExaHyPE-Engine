#include "MyElastodynamicsSolver.h"

#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"

#include "MyElastodynamicsSolver_Variables.h"

#include <algorithm>


tarch::logging::Log Elastodynamics::MyElastodynamicsSolver::_log( "Elastodynamics::MyElastodynamicsSolver" );


void Elastodynamics::MyElastodynamicsSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Elastodynamics::MyElastodynamicsSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions             = 2
  // Number of variables    = 8 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  static tarch::logging::Log _log("MyElastodynamicsSolver::adjustPointSolution");
  
  if tarch::la::equals(t,0.0) {
    Variables vars(Q);

    vars.v(0,0);
    vars.sigma(0,0,0);

    vars.rho() = 2.6;     // gm/cm^3
    vars.cs() = 2.0;      // km/s
    vars.cp() = 4.0;      // km/s

    if (x[0] > 1.0) {
    vars.rho() = 2.7;     // gm/cm^3
    vars.cs() = 3.464;    // km/s
    vars.cp() = 6.0;      // km/s
    }
    double x0 = 2.0;
    double y0 = 15.0;

    Q[0] = std::exp(-((x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0))/1.0);
    Q[1] = std::exp(-((x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0))/1.0);
    Q[2] = 0;
    Q[3] = 0;
    Q[4] = 0;
  }
}

void Elastodynamics::MyElastodynamicsSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 8 + #parameters
  
  // @todo Please implement/augment if required

  double rho  = Q[5];   // km/s
  double cs   = Q[6];   // km/s
  double cp   = Q[7];   // km/s

  //std::cout << rho << cs << cp <<'\n';

  //std::exit(-1);

  //
  lambda[0] = cp;
  lambda[1] = cs;
  lambda[2] = -cp;
  lambda[3] = -cs;
  lambda[4] = 0.0;
}


void Elastodynamics::MyElastodynamicsSolver::flux(const double* const Q,double** const F) {
  // Dimensions             = 2
  // Number of variables    = 8 + #parameters
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;

  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
}


void Elastodynamics::MyElastodynamicsSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 2
  // Number of variables    = 8 + #parameters

  // @todo Please implement/augment if required

  double rho  = stateIn[5];   // km/s
  double cs   = stateIn[6];   // km/s
  double cp   = stateIn[7];   // km/s

  


  double n[2] = {0,0};
  n[normalNonZero] = 1.;

  // extract local s-wave and p-wave impedances
  double zp=rho*cp;
  double zs=rho*cs;
  
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];

  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];

  stateOut[5] = stateIn[5];
  stateOut[6] = stateIn[6];
  stateOut[7] = stateIn[7];

  if (faceIndex == 0) {
    
    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = stateIn[4];

    stateOut[5] = stateIn[5];
    stateOut[6] = stateIn[6];
    stateOut[7] = stateIn[7];
    
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];
   
    
    double v_x =  stateIn[0];
    double v_y =  stateIn[1];
    
    double sigma_xx =  stateIn[2];
    double sigma_yy =  stateIn[3];
    
    double sigma_xy =  stateIn[4];
    
    // extract traction: Tx, Ty, Tz
    // double Tx =  0.;
    // double Ty =  0.;
    // double Tz =  0.;
    
    double vx_hat =  0.;
    double vy_hat =  0.;
    
    double sigma_xx_hat =  0.;
    double sigma_xy_hat =  0.;
    
    double r = 1.;
    
    riemannSolver_BC0(v_x, sigma_xx, zp, r, vx_hat, sigma_xx_hat);
    riemannSolver_BC0(v_y, sigma_xy, zs, r, vy_hat, sigma_xy_hat);
    
    stateOut[0] = vx_hat;
    stateOut[1] = vy_hat;
    
    stateOut[2] = sigma_xx_hat;
    stateOut[4] = sigma_xy_hat;
    
  }
  
  
  if (faceIndex == 1) {

    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = stateIn[4];
    
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];
    
    
    double v_x =  stateIn[0];
    double v_y =  stateIn[1];
    
    double sigma_xx =  stateIn[2];
    double sigma_yy =  stateIn[3];
    
    double sigma_xy =  stateIn[4];
    
    // extract traction: Tx, Ty, Tz
    // double Tx =  0.;
    // double Ty =  0.;
    // double Tz =  0.;
    
    double vx_hat =  0.;
    double vy_hat =  0.;
    
    double sigma_xx_hat =  0.;
    double sigma_xy_hat =  0.;
    
    double r = 0.;
    
    riemannSolver_BCn(v_x, sigma_xx, zp, r, vx_hat, sigma_xx_hat);
    riemannSolver_BCn(v_y, sigma_xy, zs, r, vy_hat, sigma_xy_hat);
    
    stateOut[0] = vx_hat;
    stateOut[1] = vy_hat;
    
    stateOut[2] = sigma_xx_hat;
    stateOut[4] = sigma_xy_hat;
    
      
  }


  if (faceIndex == 2) {

    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = stateIn[4];
    
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];
    
    
    double v_x =  stateIn[0];
    double v_y =  stateIn[1];
    
    double sigma_xx =  stateIn[2];
    double sigma_yy =  stateIn[3];
    
    double sigma_xy =  stateIn[4];
    
    // extract traction: Tx, Ty, Tz
    // double Tx =  0.;
    // double Ty =  0.;
    // double Tz =  0.;
    
    double vx_hat =  0.;
    double vy_hat =  0.;
    
    double sigma_yy_hat =  0.;
    double sigma_xy_hat =  0.;
    
    double r = 0.;
    
    riemannSolver_BC0(v_x, sigma_xy, zs, r, vx_hat, sigma_xy_hat);
    riemannSolver_BC0(v_y, sigma_yy, zp, r, vy_hat, sigma_yy_hat);
    
    stateOut[0] = vx_hat;
    stateOut[1] = vy_hat;
    
    stateOut[3] = sigma_yy_hat;
    stateOut[4] = sigma_xy_hat;
    
  }

  if (faceIndex == 3) {
    
    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = stateIn[4];
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];
    
    double v_x =  stateIn[0];
    double v_y =  stateIn[1];
    
    double sigma_xx =  stateIn[2];
    double sigma_yy =  stateIn[3];
    
    double sigma_xy =  stateIn[4];
    
    // extract traction: Tx, Ty, Tz
    // double Tx =  0.;
    // double Ty =  0.;
    // double Tz =  0.;
    
    double vx_hat =  0.;
    double vy_hat =  0.;
    
    double sigma_yy_hat =  0.;
    double sigma_xy_hat =  0.;
    
    double r = 0.;
    
    riemannSolver_BCn(v_x, sigma_xy, zs, r, vx_hat, sigma_xy_hat);
    riemannSolver_BCn(v_y, sigma_yy, zp, r, vy_hat, sigma_yy_hat);
    
    stateOut[0] = vx_hat;
    stateOut[1] = vy_hat;
    
    stateOut[3] = sigma_yy_hat;
    stateOut[4] = sigma_xy_hat;
    
  }
}


exahype::solvers::Solver::RefinementControl Elastodynamics::MyElastodynamicsSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Elastodynamics::MyElastodynamicsSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ){
  kernels::idx2 idx(DIMENSIONS, NumberOfVariables);

  static tarch::logging::Log _log("Elastodynamics::MyElastodynamicsSolver::nonConservativeProduct");
  ReadOnlyVariables vars(Q);

  const double* Qx = &gradQ[0]; // (:,1)
  const double* Qy = &gradQ[5]; // (:,2)

  const double rho  = Q[5];   // km/s
  const double cs   = Q[6];   // km/s
  const double cp   = Q[7];   // km/s

  
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp-2.0*mu;
  
  BgradQ[0] = -1/rho*Qx[2];
  BgradQ[1] = -1/rho*Qx[4];
  BgradQ[2] = -(lam + 2.0*mu)*Qx[0];
  BgradQ[3] = -0.0;
  BgradQ[4] = -mu*Qx[1];

  BgradQ[5] = -1/rho*Qy[4];
  BgradQ[6] = -1/rho*Qy[3];
  BgradQ[7] = -0.0;
  BgradQ[8] = -(lam + 2*mu)*Qy[1];
  BgradQ[9] = -mu*Qy[0];

}



void Elastodynamics::MyElastodynamicsSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn){
  static tarch::logging::Log _log("MyElastodynamicsSolver::coefficientMatrix");

  const double rho  = Q[5];   // km/s
  const double cs   = Q[6];   // km/s
  const double cp   = Q[7];   // km/s
    
  double mu = rho*cs*cs;
  double lam = rho*cp*cp - 2*mu;
  
  double nv[2] = {0.0, 0.0};

  nv[d] = 1.0;
  
  double B1[5][5];
  double B2[5][5];
   
  for(int i=0; i<5; i++) {
    for(int j=0; j<5; j++) {
      B1[i][j] = 0.0;
      B2[i][j] = 0.0;
    }
  }
  
  B1[0][2] = -1/rho; 
  B1[1][4] = -1/rho;
  B1[2][0] = -(2*mu + lam); 
  B1[4][1] = -mu;

  B2[0][4] = -1/rho; 
  B2[1][3] = -1/rho;
  B2[3][1] = -(2*mu + lam); 
  B2[4][0] = -mu;
  
  for(int i=0; i<5; i++) {
    for(int j=0; j<5; j++) {
      
      Bn[i*5 + j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
    }
  }

}


void Elastodynamics::MyElastodynamicsSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double t,const double dt, const tarch::la::Vector<DIMENSIONS, double>& lengthScale,const int direction, bool isBoundaryFace, int faceIndex){

  constexpr int numberOfVariables  = MyElastodynamicsSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyElastodynamicsSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyElastodynamicsSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx2 idx_QLR(basisSize, numberOfData);

  kernels::idx2 idx_FLR(basisSize, NumberOfVariables);

  double n[2]={0,0};
  n[normalNonZeroIndex]=1;

  
  for (int i = 0; i < basisSize; i++) {


    double rho_m  = QL[idx_QLR(i,5)];   // km/s
    double cs_m   = QL[idx_QLR(i,6)];   // km/s
    double cp_m   = QL[idx_QLR(i,7)];   // km/s

    double rho_p  = QR[idx_QLR(i,5)];   // km/s
    double cs_p   = QR[idx_QLR(i,6)];   // km/s
    double cp_p   = QR[idx_QLR(i,7)];   // km/s
    
    //std::cout << rho_p << cs_p << cp_p <<numberOfData<<'\n';

    //std::exit(-1);
    
    double mu_m = rho_m*cs_m*cs_m;
    double lam_m = rho_m*cp_m*cp_m - 2*mu_m;

    double mu_p = rho_p*cs_p*cs_p;
    double lam_p = rho_p*cp_p*cp_p - 2*mu_p;
  


    // extract tractions and particle velocities
    double sigma_m_xx =  QL[idx_QLR(i,2)];
    double sigma_m_yy =  QL[idx_QLR(i,3)];
    double sigma_m_xy =  QL[idx_QLR(i,4)];

    double sigma_p_xx =  QR[idx_QLR(i,2)];
    double sigma_p_yy =  QR[idx_QLR(i,3)];
    double sigma_p_xy =  QR[idx_QLR(i,4)];

    double Tx_m = n[0]*sigma_m_xx + n[1]*sigma_m_xy;
    double Ty_m = n[0]*sigma_m_xy + n[1]*sigma_m_yy;

    double Tx_p = n[0]*sigma_p_xx + n[1]*sigma_p_xy;
    double Ty_p = n[0]*sigma_p_xy + n[1]*sigma_p_yy;

    double vx_m = QL[idx_QLR(i,0)];
    double vy_m = QL[idx_QLR(i,1)];

    double vx_p = QR[idx_QLR(i,0)];
    double vy_p = QR[idx_QLR(i,1)];

    
    // generate local orthonormal basis: n, m
    double m[3] = {0,0,0};
    double l[3] = {0,0,0};
    
    localBasis(n, m, l, 2);
      //m[0] = n[1];
      //m[1] = -n[0];

    // rotate tractions and particle velocities into orthogonal coordinates: n, m
    double Tn_m= Tx_m*n[0] + Ty_m*n[1];
    double Tm_m= Tx_m*m[0] + Ty_m*m[1];
    
    double Tn_p= Tx_p*n[0] + Ty_p*n[1];
    double Tm_p= Tx_p*m[0] + Ty_p*m[1];
    
    double vn_m= vx_m*n[0] + vy_m*n[1];
    double vm_m= vx_m*m[0] + vy_m*m[1];
    
    double vn_p= vx_p*n[0] + vy_p*n[1];
    double vm_p= vx_p*m[0] + vy_p*m[1];

   
    // extract local s-wave and p-wave impedances
    double zs_p=rho_p*cs_p;
    double zs_m=rho_m*cs_m;

    double zp_p=rho_p*cp_p;
    double zp_m=rho_m*cp_m;

    // impedance must be greater than zero !
    if (zs_p <= 0.0 || zs_m <= 0.0 || zp_p <= 0.0 || zp_m <= 0.0){
      std::cout<<zs_p<<' '<<zs_m<<' '<<zp_p<<' '<<zp_m<<'\n';
      std::cout<<' Impedance must be greater than zero ! '<<'\n';
      std::exit(-1);
    }

    
    // generate interface data preserving the amplitude of the outgoing charactertritics
    // and satisfying interface conditions exactly.
    double vn_hat_p=0;
    double vm_hat_p=0;

    double vn_hat_m=0;
    double vm_hat_m=0;

    double Tn_hat_p=0;
    double Tm_hat_p=0;

    double Tn_hat_m=0;
    double Tm_hat_m=0;

    // data is generated by solving a local Riemann problem and contraining the solutions against
    // physical interface conditions
    riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
    riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);

    // generate fluctuations in the local basis coordinates: n, m
    double FLn = 0.5*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
    double FLm = 0.5*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));

    double FRn = 0.5*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
    double FRm = 0.5*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));

    double FL_n = 0.5/zp_m*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
    double FL_m = 0.5/zs_m*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));

    double FR_n = 0.5/zp_p*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
    double FR_m = 0.5/zs_p*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));

    // rotate back to the physical coordinates x, y
    double FLx = n[0]*FLn + m[0]*FLm;
    double FLy = n[1]*FLn + m[1]*FLm;

    double FRx = n[0]*FRn + m[0]*FRm;
    double FRy = n[1]*FRn + m[1]*FRm;

    double FL_x = n[0]*FL_n + m[0]*FL_m;
    double FL_y = n[1]*FL_n + m[1]*FL_m;

    double FR_x = n[0]*FR_n + m[0]*FR_m;
    double FR_y = n[1]*FR_n + m[1]*FR_m;
     
    // construct flux fluctuation vectors obeying the eigen structure of the PDE
    // and choose physically motivated penalties such that we can prove
    // numerical stability.

    FR[idx_FLR(i, 0)] = 1/rho_p*FRx;
    FL[idx_FLR(i, 0)] = 1/rho_m*FLx;
    
    FR[idx_FLR(i, 1)] = 1/rho_p*FRy;
    FL[idx_FLR(i, 1)] = 1/rho_m*FLy;

    FL[idx_FLR(i, 2)] = (2*mu_m+lam_m)*n[0]*FL_x + lam_m*n[1]*FL_y;
    FL[idx_FLR(i, 3)] = (2*mu_m+lam_m)*n[1]*FL_y + lam_m*n[0]*FL_x;

    FR[idx_FLR(i, 2)] = -((2*mu_p+lam_p)*n[0]*FR_x + lam_p*n[1]*FR_y);
    FR[idx_FLR(i, 3)] = -((2*mu_p+lam_p)*n[1]*FR_y + lam_p*n[0]*FR_x);


    FL[idx_FLR(i, 4)] =  mu_m*(n[1]*FL_x + n[0]*FL_y);
    FR[idx_FLR(i, 4)] = -mu_p*(n[1]*FR_x + n[0]*FR_y);
    
  }
  
}


// void Elastodynamics::MyElastodynamicsSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0){
//      double pi = 3.14159265359;
//   double sigma = 0.1149;
//   double t0 = 0.7;
//   //double t0 = 0.1;
//   double f = 0.0;
//   double M0 = 1000.0;
  

//   //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

//   f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


//   //f = M0*t/(t0*t0)*std::exp(-t/t0);

//   x0[0] = 0.5;
//   x0[1] = 0.5;

//   forceVector[0] = 1.*f;
//   forceVector[1] = 0.0;
//   forceVector[2] = 0.0;

// }

void Elastodynamics::MyElastodynamicsSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0){
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));
  
  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
  
  f = M0*t/(t0*t0)*std::exp(-t/t0);
  
  x0[0] = 2.0;
  x0[1] = 15.0;
  
  forceVector[0] = 0.0;
  forceVector[1] = 0.0;
  forceVector[2] = 0;
  forceVector[3] = 0;
  forceVector[4] = f;

}


void Elastodynamics::MyElastodynamicsSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
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


void Elastodynamics::MyElastodynamicsSolver::Gram_Schmidt(double* const y, double* const z){
  //Gram Schmidt orthonormalization
 
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
    z[i] = 1.0/std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2])*z[i];
  }

}

void Elastodynamics::MyElastodynamicsSolver::localBasis(double* const n, double* const m, double* const l, int d){

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


void Elastodynamics::MyElastodynamicsSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

 }


void Elastodynamics::MyElastodynamicsSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

 }
