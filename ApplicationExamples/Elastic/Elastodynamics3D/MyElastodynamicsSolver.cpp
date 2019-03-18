#include "MyElastodynamicsSolver.h"

#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"


#include "MyElastodynamicsSolver_Variables.h"

#include <algorithm>


tarch::logging::Log Elastodynamics::MyElastodynamicsSolver::_log( "Elastodynamics::MyElastodynamicsSolver" );


void Elastodynamics::MyElastodynamicsSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Elastodynamics::MyElastodynamicsSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void Elastodynamics::MyElastodynamicsSolver::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
  // Dimensions             = 3
  // Number of variables    = 9 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  static tarch::logging::Log _log("MyElastodynamicsSolver::adjustPointSolution");

  if (tarch::la::equals(t,0.0)) {
    Variables vars(Q);

    Q[0] = 0*std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/0.01);
    Q[1] = 0*std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/0.01);
    Q[2] = 0*std::exp(-((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/0.01);

    Q[3] = 0.0;
    Q[4] = 0.0;
    Q[5] = 0.0;
    Q[6] = 0.0;
    Q[7] = 0.0;
    Q[8] = 0.0;
  }
}

void Elastodynamics::MyElastodynamicsSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions             = 3
  // Number of variables    = 9 + #parameters
  
  // @todo Please implement/augment if required
  double cp = 6.0;
  double cs = 3.343;
  
  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = -cs;
  lambda[3] = 0.0;
  lambda[4] = 0.0;
  lambda[5] = 0.0;
  lambda[6] = +cs;
  lambda[7] = +cs;
  lambda[8] = +cp;
  
}


void Elastodynamics::MyElastodynamicsSolver::flux(const double* const Q,double** const F) {
  // Dimensions             = 3
  // Number of variables    = 9 + #parameters
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  F[0][5] = 0.0;
  F[0][6] = 0.0;
  F[0][7] = 0.0;
  F[0][8] = 0.0;

  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  F[1][5] = 0.0;
  F[1][6] = 0.0;
  F[1][7] = 0.0;
  F[1][8] = 0.0;

  F[2][0] = 0.0;
  F[2][1] = 0.0;
  F[2][2] = 0.0;
  F[2][3] = 0.0;
  F[2][4] = 0.0;
  F[2][5] = 0.0;
  F[2][6] = 0.0;
  F[2][7] = 0.0;
  F[2][8] = 0.0;
}

void Elastodynamics::MyElastodynamicsSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions             = 3
  // Number of variables    = 9 + #parameters

  // @todo Please implement/augment if required

  double cp = 6.0;
  double cs = 3.343;
  
  double rho = 2.7;

  double n[3] = {0,0,0};
  n[normalNonZero] = 1.;

  // extract local s-wave and p-wave impedances
  double zp=rho*cp;
  double zs=rho*cs;
  
  stateOut[0] = 0.0;
  stateOut[1] = 0.0;
  stateOut[2] = 0.0;
  stateOut[3] = 0.0;
  stateOut[4] = 0.0;
  stateOut[5] = 0.0;
  stateOut[6] = 0.0;
  stateOut[7] = 0.0;
  stateOut[8] = 0.0;

  fluxOut[0] = 0.0;
  fluxOut[1] = 0.0;
  fluxOut[2] = 0.0;
  fluxOut[3] = 0.0;
  fluxOut[4] = 0.0;
  fluxOut[5] = 0.0;
  fluxOut[6] = 0.0;
  fluxOut[7] = 0.0;
  fluxOut[8] = 0.0;

  if (faceIndex == 0) {

    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    stateOut[3] = stateIn[3];
    stateOut[4] = stateIn[4];
    stateOut[5] = stateIn[5];
    stateOut[6] = stateIn[6];
    stateOut[7] = stateIn[7];
    stateOut[8] = stateIn[8];
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
    fluxOut[3] = fluxIn[3];
    fluxOut[4] = fluxIn[4];
    fluxOut[5] = fluxIn[5];
    fluxOut[6] = fluxIn[6];
    fluxOut[7] = fluxIn[7];
    fluxOut[8] = fluxIn[8];
    
    double v_x =  stateIn[0];
    double v_y =  stateIn[1];
    double v_z =  stateIn[2];
    
    double sigma_xx =  stateIn[3];
    double sigma_yy =  stateIn[4];
    double sigma_zz =  stateIn[5];

    double sigma_xy =  stateIn[6];
    double sigma_xz =  stateIn[7];
    double sigma_yz =  stateIn[8];

    // extract traction: Tx, Ty, Tz
    // double Tx =  0.;
    // double Ty =  0.;
    // double Tz =  0.;

    double vx_hat =  0.;
    double vy_hat =  0.;
    double vz_hat =  0.;

    double sigma_xx_hat =  0.;
    double sigma_xy_hat =  0.;
    double sigma_xz_hat =  0.;

    double r = 1.;
    
    riemannSolver_BC0(v_x, sigma_xx, zp, r, vx_hat, sigma_xx_hat);
    riemannSolver_BC0(v_y, sigma_xy, zs, r, vy_hat, sigma_xy_hat);
    riemannSolver_BC0(v_z, sigma_xz, zs, r, vz_hat, sigma_xz_hat);

    stateOut[0] = vx_hat;
    stateOut[1] = vy_hat;
    stateOut[2] = vz_hat;

    stateOut[3] = sigma_xx_hat;
    stateOut[6] = sigma_xy_hat;
    stateOut[7] = sigma_xz_hat;

    
    
      
  }
}

exahype::solvers::Solver::RefinementControl Elastodynamics::MyElastodynamicsSolver::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Elastodynamics::MyElastodynamicsSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Dimensions             = 3
  // Number of variables    = 9 (#unknowns + #parameters)
  ReadOnlyVariables vars(Q);

  static tarch::logging::Log _log("Elastodynamics::MyElastodynamicsSolver::nonConservativeProduct");

  const double* Qx = &gradQ[0]; // (:,1)
  const double* Qy = &gradQ[9]; // (:,2)
  const double* Qz = &gradQ[18]; // (:,3)

  //  std::cout << "ncp " << vars.rho() << "," << vars.cs() << "," << vars.cp() << std::endl;
  //  double cs = 3.464;  // km/s
  //  double cp = 6.0;    // km/s
  //  double rho = 2.67;  // gm/cm^3

  double cp = 6.0;
  double cs = 3.343;
  
  double rho = 2.7;
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp-2.0*mu;
  
  BgradQ[0] = -1.0/rho*Qx[3];
  BgradQ[1] = -1.0/rho*Qx[6];
  BgradQ[2] = -1.0/rho*Qx[7];
  BgradQ[3] = -(2.0*mu + lam)*Qx[0];
  BgradQ[4] = -lam*Qx[0];
  BgradQ[5] = -lam*Qx[0];
  BgradQ[6] = -mu*Qx[1];
  BgradQ[7] = -mu*Qx[2];
  BgradQ[8] =  0.0;

  BgradQ[9]  = -1.0/rho*Qy[6];
  BgradQ[10] = -1.0/rho*Qy[4];
  BgradQ[11] = -1.0/rho*Qy[8];
  BgradQ[12] = -lam*Qy[1];
  BgradQ[13] = -(lam+2.0*mu)*Qy[1];
  BgradQ[14] = -lam*Qy[1];
  BgradQ[15] = -mu*Qy[0];
  BgradQ[16] = -0.0;
  BgradQ[17] = -mu*Qy[2];

  BgradQ[18] = -1.0/rho*Qz[7];
  BgradQ[19] = -1.0/rho*Qz[8];
  BgradQ[20] = -1.0/rho*Qz[5];
  BgradQ[21] = -lam*Qz[2];
  BgradQ[22] = -lam*Qz[2];
  BgradQ[23] = -(lam+2.0*mu)*Qz[2];
  BgradQ[24] = - 0.0;
  BgradQ[25] = -mu*Qz[0];
  BgradQ[26] = -mu*Qz[1];
}


void Elastodynamics::MyElastodynamicsSolver::coefficientMatrix(const double* const Q,const int normalNonZero,double* const Bn) {
  ReadOnlyVariables vars(Q);
  
  //  std::cout << "matrixb " << vars.rho() << "," << vars.cs() << "," << vars.cp() << std::endl;
  double cs = 3.343;  // km/s
  double cp = 6.0;    // km/s
  double rho = 2.7;   // gm/cm^3
  
  double lam = rho * (cp * cp - 2 * cs * cs);  // lambda (MPa)
  double mu= rho * cs * cs;

  double nv[3] = {0.0};
  nv[normalNonZero] = 1.0;
  
  double B1[9][9];
  double B2[9][9];
  double B3[9][9];

  for(int i=0; i<9; i++) {
   for(int j=0; j<9; j++) {
     B1[i][j] = 0.0;
     B2[i][j] = 0.0;
     B3[i][j] = 0.0;
   }
 }

  B1[3][0] = -1.0/rho; 
  B1[6][1] = -1.0/rho;
  B1[7][2] = -1.0/rho;

  B1[0][3] = -(2.0*mu + lam); 
  B1[0][4] = -lam;
  B1[0][5] = -lam;

  B1[1][6] = -mu; 
  B1[2][7] = -mu;

  B2[6][0] = -1.0/rho; 
  B2[4][1] = -1.0/rho;
  B2[8][2] = -1.0/rho;

  B2[1][3] = -lam;
  B2[1][4] = -(2.0*mu + lam);
  B2[1][5] = -lam;

  B2[0][6] = -mu; 
  B2[2][8] = -mu;

  B3[7][0] = -1.0/rho; 
  B3[8][1] = -1.0/rho;
  B3[5][2] = -1.0/rho;

  B3[2][3] = -lam;
  B3[2][4] = -lam;
  B3[2][5] = -(2.0*mu + lam);

  B3[0][7] = -mu; 
  B3[1][8] = -mu;


  for(int i=0; i<9; i++) {
   for(int j=0; j<9; j++) {
      Bn[i*9 + j] = nv[0]*B1[i][j] + nv[1]*B2[i][j] + nv[2]*B3[i][j];
   }
 }
}






void Elastodynamics::MyElastodynamicsSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,double* const tempFaceUnknownsArray,double** const tempStateSizedVectors,double** const tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex){

  constexpr int numberOfVariables  = MyElastodynamicsSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyElastodynamicsSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyElastodynamicsSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx3 idx_QLR(basisSize, basisSize, numberOfVariables);

  kernels::idx3 idx_FLR(basisSize, basisSize, numberOfVariables);

  double n[3]={0,0,0};
  n[normalNonZeroIndex]=1.;

  double cp = 6.0;
  double cs = 3.343;
  
  double rho = 2.7;
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp - 2*mu;
  
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {

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


      double Tx_m = n[0]*sigma_m_xx + n[1]*sigma_m_xy + n[2]*sigma_m_xz;
      double Ty_m = n[0]*sigma_m_xy + n[1]*sigma_m_yy + n[2]*sigma_m_yz;
      double Tz_m = n[0]*sigma_m_xz + n[1]*sigma_m_yz + n[2]*sigma_m_zz;

      double Tx_p = n[0]*sigma_p_xx + n[1]*sigma_p_xy + n[2]*sigma_p_xz;
      double Ty_p = n[0]*sigma_p_xy + n[1]*sigma_p_yy + n[2]*sigma_p_yz;
      double Tz_p = n[0]*sigma_p_xz + n[1]*sigma_p_yz + n[2]*sigma_p_zz;


      double vx_m = QL[idx_QLR(i,j,0)];
      double vy_m = QL[idx_QLR(i,j,1)];
      double vz_m = QL[idx_QLR(i,j,2)];

      double vx_p = QR[idx_QLR(i,j,0)];
      double vy_p = QR[idx_QLR(i,j,1)];
      double vz_p = QR[idx_QLR(i,j,2)];

    
    // generate local orthonormal basis: n, m, l
    double m[3] = {0,0,0};
    double l[3] = {0,0,0};
    
    localBasis(n, m, l, 3);
      //m[0] = n[1];
      //m[1] = -n[0];

    // rotate tractions and particle velocities into orthogonal coordinates: n, m
    double Tn_m= Tx_m*n[0] + Ty_m*n[1] + Tz_m*n[2];
    double Tm_m= Tx_m*m[0] + Ty_m*m[1] + Tz_m*m[2];
    double Tl_m= Tx_m*l[0] + Ty_m*l[1] + Tz_m*l[2];

    double Tn_p= Tx_p*n[0] + Ty_p*n[1] + Tz_p*n[2];
    double Tm_p= Tx_p*m[0] + Ty_p*m[1] + Tz_p*m[2];
    double Tl_p= Tx_p*l[0] + Ty_p*l[1] + Tz_p*l[2];
    
    
    double vn_m= vx_m*n[0] + vy_m*n[1] + vz_m*n[2];
    double vm_m= vx_m*m[0] + vy_m*m[1] + vz_m*m[2];
    double vl_m= vx_m*l[0] + vy_m*l[1] + vz_m*l[2];


    double vn_p= vx_p*n[0] + vy_p*n[1] + vz_p*n[2];
    double vm_p= vx_p*m[0] + vy_p*m[1] + vz_p*m[2];
    double vl_p= vx_p*l[0] + vy_p*l[1] + vz_p*l[2];
  
   
    // extract local s-wave and p-wave impedances
    double zs_p=rho*cs;
    double zs_m=rho*cs;

    double zp_p=rho*cp;
    double zp_m=rho*cp;

    
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
    // physical interface conditions. Data must preserve the amplitude of the outgoing characteris
    // and satisfy the physical interface conditions exactly
    riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
    riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
    riemannSolver_Nodal(vl_p,vl_m, Tl_p, Tl_m, zs_p , zs_m, vl_hat_p , vl_hat_m, Tl_hat_p, Tl_hat_m);

    // generate fluctuations in the local basis coordinates: n, m, l.
    double FLn = 0.5*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
    double FLm = 0.5*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));
    double FLl = 0.5*(zs_m*(vl_m-vl_hat_m) + (Tl_m-Tl_hat_m));

    double FRn = 0.5*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
    double FRm = 0.5*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));
    double FRl = 0.5*(zs_p*(vl_p-vl_hat_p) - (Tl_p-Tl_hat_p));

    double FL_n = 0.5/zp_m*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
    double FL_m = 0.5/zs_m*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));
    double FL_l = 0.5/zs_m*(zs_m*(vl_m-vl_hat_m) + (Tl_m-Tl_hat_m));

    double FR_n = 0.5/zp_p*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
    double FR_m = 0.5/zs_p*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));
    double FR_l = 0.5/zs_p*(zs_p*(vl_p-vl_hat_p) - (Tl_p-Tl_hat_p));

    
    // rotate back to the physical coordinates x, y, z
    double FLx = n[0]*FLn + m[0]*FLm + l[0]*FLl;
    double FLy = n[1]*FLn + m[1]*FLm + l[1]*FLl;
    double FLz = n[2]*FLn + m[2]*FLm + l[2]*FLl;

    double FRx = n[0]*FRn + m[0]*FRm + l[0]*FRl;
    double FRy = n[1]*FRn + m[1]*FRm + l[1]*FRl;
    double FRz = n[2]*FRn + m[2]*FRm + l[2]*FRl;


    double FL_x = n[0]*FL_n + m[0]*FL_m + l[0]*FL_l;
    double FL_y = n[1]*FL_n + m[1]*FL_m + l[1]*FL_l;
    double FL_z = n[2]*FL_n + m[2]*FL_m + l[2]*FL_l;

    double FR_x = n[0]*FR_n + m[0]*FR_m + l[0]*FR_l;
    double FR_y = n[1]*FR_n + m[1]*FR_m + l[1]*FR_l;
    double FR_z = n[2]*FR_n + m[2]*FR_m + l[2]*FR_l;
     
    // construct flux fluctuation vectors obeying the eigen structure of the PDE
    // and choose physically motivated penalties such that we can prove
    // numerical stability.

    FR[idx_FLR(i, j, 0)] = 1/rho*FRx;
    FL[idx_FLR(i, j, 0)] = 1/rho*FLx;
    
    FR[idx_FLR(i, j, 1)] = 1/rho*FRy;
    FL[idx_FLR(i, j, 1)] = 1/rho*FLy;

    FR[idx_FLR(i, j, 2)] = 1/rho*FRz;
    FL[idx_FLR(i, j, 2)] = 1/rho*FLz;

    FL[idx_FLR(i, j, 3)] = (2*mu+lam)*n[0]*FL_x + lam*n[1]*FL_y + lam*n[2]*FL_z;
    FL[idx_FLR(i, j, 4)] = (2*mu+lam)*n[1]*FL_y + lam*n[0]*FL_x + lam*n[2]*FL_z;
    FL[idx_FLR(i, j, 5)] = (2*mu+lam)*n[2]*FL_z + lam*n[0]*FL_x + lam*n[1]*FL_y;

    FR[idx_FLR(i, j, 3)] = -((2*mu+lam)*n[0]*FR_x + lam*n[1]*FR_y + lam*n[2]*FR_z);
    FR[idx_FLR(i, j, 4)] = -((2*mu+lam)*n[1]*FR_y + lam*n[0]*FR_x + lam*n[2]*FR_z);
    FR[idx_FLR(i, j, 5)] = -((2*mu+lam)*n[2]*FR_z + lam*n[0]*FR_x + lam*n[1]*FR_y);


    FL[idx_FLR(i, j, 6)] =  mu*(n[1]*FL_x + n[0]*FL_y);
    FL[idx_FLR(i, j, 7)] =  mu*(n[2]*FL_x + n[0]*FL_z);
    FL[idx_FLR(i, j, 8)] =  mu*(n[1]*FL_z + n[2]*FL_y);

    FR[idx_FLR(i, j, 6)] =  -mu*(n[1]*FR_x + n[0]*FR_y);
    FR[idx_FLR(i, j, 7)] =  -mu*(n[2]*FR_x + n[0]*FR_z);
    FR[idx_FLR(i, j, 8)] =  -mu*(n[1]*FR_z + n[2]*FR_y);
    
    }
  }
}



void Elastodynamics::MyElastodynamicsSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0) {
  //TODO KD // @todo Please implement/augment if required and set bool function
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //  double t0 = 0.1;
  double f;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

  x0[0] = 0.5;
  x0[1] = 0.5;
  x0[2] = 0.5;

  //double exp0 = std::exp(-sqrt((x[0] - x0[0]) * (x[0] - x0[0]) + (x[1] - x0[1]) * (x[1] - x0[1]) +
  //                 (x[2] - x0[2]) * (x[2] - x0[2])) /((0.25)*(0.25)));

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
