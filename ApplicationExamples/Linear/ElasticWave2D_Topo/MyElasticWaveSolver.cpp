#include "MyElasticWaveSolver.h"
#include "MyElasticWaveSolver_Variables.h"

#include "CurvilinearTransformation.h"



tarch::logging::Log Elastic::MyElasticWaveSolver::_log( "Elastic::MyElasticWaveSolver" );


void Elastic::MyElasticWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Elastic::MyElasticWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 10
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    constexpr int basisSize = MyElasticWaveSolver::Order+1;
    int numberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

    kernels::idx3 id_xyf(basisSize,basisSize,numberOfData);
    kernels::idx2 id_xy(basisSize,basisSize);

    double offset_x=center[0]-0.5*dx[0];
    double offset_y=center[1]-0.5*dx[1];
    double width_x=dx[0];
    double width_y=dx[1];
    
    int num_nodes = basisSize;

    //Number of elements on the unit square
    int ne_x = std::round(1/dx[0]);
    int ne_y = std::round(1/dx[1]);

    //Number of nodes on the unit square
    int nx = ne_x *(num_nodes-1) + 1;
    int ny = ne_y *(num_nodes-1) + 1;

    //first local node indeces in global domain
    int i_m;
    int j_m;

    // Coordiantes of the undeformed domain (rectangular)
    double a_x = 0;
    double a_y = 0;
    double b_x = 10.0;
    double b_y = 5.0;

    // Position of interface in x
    double fault_position = 5.0;

    // Width of the current block
    double blockWidth_x;
    double blockWidth_y=(b_y-a_y); //equal to width of the rectangle
    
    double recwidth_x  =(b_x-a_x);  //width in x of the rectangle

    //relative position of the fault on the unit square
    double fault_ref = (fault_position-a_x)/(recwidth_x);

    // number of the current block
    int n =  offset_x >= fault_ref ? 1 : 0 ;
     
    if(n == 0){
      blockWidth_x=(fault_position-a_x);
      ne_x = std::round((ne_x+1)*fault_ref);
      nx   = ne_x*(num_nodes-1)+1;
      
      i_m =  std::round((offset_x)/width_x) *(num_nodes-1);
      j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
     }else{
      double ne_x0 = std::round((ne_x+1)*fault_ref);
      
      ne_x = std::round(1/dx[0])-ne_x0;
      nx =  ne_x *(num_nodes-1)+1;
      blockWidth_x= (b_x-fault_position);
       
      i_m =  std::floor((offset_x-fault_ref)/width_x) *(num_nodes-1);
      j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    }

     double* left_bnd_x = new double[ny];
     double* left_bnd_y = new double[ny];
     double* right_bnd_x = new double[ny];
     double* right_bnd_y = new double[ny];
     double* bottom_bnd_x = new double[nx];
     double* bottom_bnd_y = new double[nx];
     double* top_bnd_x = new double[nx];
     double* top_bnd_y = new double[nx];

     getBoundaryCurvesForBlock(num_nodes,nx, ny, n, fault_position, a_x, a_y, b_x, b_y,blockWidth_x,blockWidth_y ,left_bnd_x,left_bnd_y,right_bnd_x,right_bnd_y,bottom_bnd_x,bottom_bnd_y,top_bnd_x,top_bnd_y);

     double* curvilinear_x = new double[num_nodes*num_nodes];
     double* curvilinear_y = new double[num_nodes*num_nodes];
     
     int i_p = i_m + num_nodes;
     int j_p = j_m + num_nodes; 
     
     transFiniteInterpolation(nx, ny, j_m, j_p, i_m, i_p, num_nodes,  left_bnd_x,  left_bnd_y,  right_bnd_x,  right_bnd_y,  bottom_bnd_x,  bottom_bnd_y,  top_bnd_x,  top_bnd_y,  curvilinear_x ,  curvilinear_y );
     
     delete left_bnd_x;
     delete left_bnd_y;     
     delete right_bnd_x;
     delete right_bnd_y;
     delete bottom_bnd_x;
     delete bottom_bnd_y;
     delete top_bnd_x;
     delete top_bnd_y; 
     
     double* gl_vals_x = new double[num_nodes*num_nodes];
     double* gl_vals_y = new double[num_nodes*num_nodes];
     double* jacobian = new double[num_nodes*num_nodes];
     double* q_x = new double[num_nodes*num_nodes];
     double* q_y = new double[num_nodes*num_nodes];
     double* r_x = new double[num_nodes*num_nodes];
     double* r_y = new double[num_nodes*num_nodes];

     metricDerivativesAndJacobian(num_nodes,curvilinear_x,curvilinear_y,gl_vals_x,gl_vals_y,q_x,q_y,r_x,r_y,jacobian, width_x, width_y);
     
     for (int j=0; j< basisSize; j++){
       for (int i=0; i< basisSize; i++){
	 double x  =  gl_vals_x[id_xy(i,j)];
	 double y  =  gl_vals_y[id_xy(i,j)];
	 
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

	// Jacobian
	luh[id_xyf(i,j,8)] = jacobian[id_xy(i,j)];
	
	// metric derivatives
	luh[id_xyf(i,j,9)]  = q_x[id_xy(i,j)];
	luh[id_xyf(i,j,10)] = q_y[id_xy(i,j)];
	luh[id_xyf(i,j,11)] = r_x[id_xy(i,j)];
	luh[id_xyf(i,j,12)] = r_y[id_xy(i,j)];
	
	// curvilinear mesh
	luh[id_xyf(i,j,13)] = gl_vals_x[id_xy(i,j)];
	luh[id_xyf(i,j,14)] = gl_vals_y[id_xy(i,j)];
      }
    }
     
     delete gl_vals_x;
     delete gl_vals_y;
     delete jacobian;
     delete q_x; 
     delete q_y;
     delete r_x; 
     delete r_y; 
  }
}

void Elastic::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 10

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
  // Number of variables + parameters  = 5 + 10
  
  double cp = Q[6];
  double cs = Q[7];

  double q_x;
  double q_y;
  double r_x;
  double r_y;
  
  extractTransformation(Q,q_x,q_y,r_x,r_y);
   
  lambda[0] = std::sqrt(q_x*q_x + q_y*q_y)*cp;
  lambda[1] = std::sqrt(q_x*q_x + q_y*q_y)*cs;
  lambda[2] = std::sqrt(r_x*r_x + r_y*r_y)*cp;
  lambda[3] = std::sqrt(r_x*r_x + r_y*r_y)*cs;
  lambda[4] = 0.0;
}


void Elastic::MyElasticWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 10
  
  double sigma_xx=Q[2];
  double sigma_yy=Q[3];
  double sigma_xy=Q[4];

  double jacobian=Q[8];
  double q_x;
  double q_y;
  double r_x;
  double r_y;
  
  extractTransformation(Q,q_x,q_y,r_x,r_y);
  
  F[0][0] = -jacobian*(q_x*sigma_xx+q_y*sigma_xy);
  F[0][1] = -jacobian*(q_x*sigma_xy+q_y*sigma_yy);
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;

  F[1][0] = -jacobian*(r_x*sigma_xx+r_y*sigma_xy);
  F[1][1] = -jacobian*(r_x*sigma_xy+r_y*sigma_yy);
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
}



void  Elastic::MyElasticWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  
  double q_x;
  double q_y;
  double r_x;
  double r_y;

  double u_q=gradQ[0];    
  double v_q=gradQ[1];
  double u_r=gradQ[5];
  double v_r=gradQ[6];

  extractTransformation(Q,q_x,q_y,r_x,r_y);
  
  BgradQ[0] = 0;
  BgradQ[1] = 0;
  BgradQ[2] = -q_x*u_q;
  BgradQ[3] = -q_y*v_q;
  BgradQ[4] = -(q_y*u_q+q_x*v_q);
  
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = -r_x*u_r;
  BgradQ[8] = -r_y*v_r;
  BgradQ[9] = -(r_y*u_r+r_x*v_r);
}

void  Elastic::MyElasticWaveSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0, int n) {
  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;

  double a_x = 0.0;
  double a_y = 0.0;
  double b_x = 10.0;
  double b_y = 5.0;
 
  double blockWidth_y = (b_y-a_y);
  double blockWidth_x = (b_x-a_x);

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

  if(n == 0){
    double x1 = 5.0;
    double y1 = 2.5;
    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);

    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;

    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = f;
    forceVector[3] = f;
    forceVector[4] = 0.0;
    
  }else if(n == 1){
    double x1 = 5.0;
    double y1 = 2.5;
    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);

    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.*f;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.0;
  }    
}


    /**
     * @TODO LR : document
     */
void Elastic::MyElasticWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  const double rho  = Q[5];   // km/s
  const double cp   = Q[6];   // km/s
  const double cs   = Q[7];   // km/s
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp-2.0*mu;
  double rho_jacobian = Q[8]*Q[5];

  rhs[0] = 1.0/rho_jacobian*rhs[0];
  rhs[1] = 1.0/rho_jacobian*rhs[1];
  
  double rhs_2= (2*mu+lam)*rhs[2]+lam*rhs[3];
  double rhs_3= (2*mu+lam)*rhs[3]+lam*rhs[2];
  
  rhs[2]=rhs_2;
  rhs[3]=rhs_3;  
  rhs[4]=mu*rhs[4];


  rhs[5] = 1.0/rho_jacobian*rhs[5];
  rhs[6] = 1.0/rho_jacobian*rhs[6];
  
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
  
  double rho_p,cp_p,cs_p, mu_p,lam_p;
  double rho_m,cp_m,cs_m, mu_m,lam_m;  

  double norm_p_qr, norm_m_qr;
  
  double n_p[2]={0,0};
  double n_m[2]={0,0};
  double m_p[2]={0,0};
  double m_m[2]={0,0};
  
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
    
    get_normals(normalNonZeroIndex, norm_p_qr, n_p, QR + idx_QLR(i,0));
    get_normals(normalNonZeroIndex, norm_m_qr, n_m, QL + idx_QLR(i,0));    
    
    double Tx_m,Ty_m,vx_m,vy_m;
    double Tx_p,Ty_p,vx_p,vy_p;

    extract_tractions_and_particle_velocity(n_p,QR+idx_QLR(i,0),Tx_p,Ty_p,vx_p,vy_p);
    extract_tractions_and_particle_velocity(n_m,QL+idx_QLR(i,0),Tx_m,Ty_m,vx_m,vy_m); 
    
    localBasis(n_p, m_p);
    localBasis(n_m, m_m);
    
    double Tn_m,Tm_m,vn_m,vm_m;
    double Tn_p,Tm_p,vn_p,vm_p;
    
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
    
    double vn_hat_p,vm_hat_p,vn_hat_m,vm_hat_m;        
    double Tn_hat_p,Tm_hat_p,Tn_hat_m,Tm_hat_m;

    double FLn,FLm,FRn,FRm; 
    double FL_n,FL_m,FR_n,FR_m;
    double FLx,FLy,FRx,FRy;
    double FL_x,FL_y,FR_x,FR_y;
    
    if (isBoundaryFace) {
      double r= faceIndex==3 ? 1 : 0;
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
    if(zs_m > 0){
      FL_m = FLm/zs_m;
    }else{
      FL_m=0;
    }
    
    FR_n = FRn/zp_p;
    if(zs_p > 0){    
      FR_m = FRm/zs_p;
    }else{
      FR_m=0;
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

void Elastic::MyElasticWaveSolver::extractTransformation(const double* const Q, double& q_x,double& q_y,double& r_x,double& r_y) {
  q_x     =Q[9];
  q_y     =Q[10];
  r_x     =Q[11];
  r_y     =Q[12];
}

void Elastic::MyElasticWaveSolver::get_normals(int normalNonZeroIndex,double& norm, double* const n,const double* const Q){
  double q_x;
  double q_y;
  double r_x;
  double r_y;
  
  extractTransformation(Q,q_x,q_y,r_x,r_y);
  if (normalNonZeroIndex == 0){
    norm = std::sqrt(q_x*q_x + q_y*q_y);
    n[0] = q_x/norm;
    n[1] = q_y/norm;
  }
  if (normalNonZeroIndex == 1){
    norm = std::sqrt(r_x*r_x + r_y*r_y);
    n[0] = r_x/norm;
    n[1] = r_y/norm;
  }
}

