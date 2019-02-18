#include "MyLinearWaveSolver.h"

#include "MyLinearWaveSolver_Variables.h"

#include <algorithm>

#include "CurvilinearTransformation.h"


tarch::logging::Log Linear::MyLinearWaveSolver::_log( "Linear::MyLinearWaveSolver" );


void Linear::MyLinearWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Linear::MyLinearWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 9
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    //initialise luh

      constexpr int basisSize = MyLinearWaveSolver::Order+1;
    int numberOfData=MyLinearWaveSolver::NumberOfParameters+MyLinearWaveSolver::NumberOfVariables;


     kernels::idx3 id_xyf(basisSize,basisSize,numberOfData);
     kernels::idx2 id_xy(basisSize,basisSize);

     int num_nodes = basisSize;
     
    
     
     int ne_x = std::round(1/dx[0]);
     int ne_y = std::round(1/dx[1]);


     
     
     int nx = ne_x *(num_nodes-1) + 1;
     int ny = ne_y *(num_nodes-1) + 1;

     double offset_x=center[0]-0.5*dx[0];
     double offset_y=center[1]-0.5*dx[1];
    

     double width_x=dx[0];
     double width_y=dx[1];


     int i_m;
     int j_m;


     double a_x = 0;
     double a_y = 0;

     double b_x = 10.0;
     double b_y = 5.0;
     double fault_position = 5.0;


     double blockWidth_x;
     double blockWidth_y=(b_y-a_y);

     double block_width0_x = (b_x-a_x);

     
     
     double fault_ref = (fault_position-a_x)/(block_width0_x);

     //int n = center[0] >= fault_ref ? 1 : 0 ;

      int n =  offset_x >= fault_ref ? 1 : 0 ;

     
     if(n == 0){
       
       blockWidth_x=(fault_position-a_x);
       ne_x = std::round((ne_x+1)*fault_ref);
       nx =  ne_x *(num_nodes-1)+1;
       
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
     
     double* gl_vals_x = new double[num_nodes*num_nodes];
     double* gl_vals_y = new double[num_nodes*num_nodes];
     double* jacobian = new double[num_nodes*num_nodes];
     
     
     double* q_x = new double[num_nodes*num_nodes];
     double* q_y = new double[num_nodes*num_nodes];
     
     double* r_x = new double[num_nodes*num_nodes];
     double* r_y = new double[num_nodes*num_nodes];

     metricDerivativesAndJacobian(num_nodes,curvilinear_x,curvilinear_y,gl_vals_x,gl_vals_y,q_x,q_y,r_x,r_y,jacobian, width_x, width_y);

    for (int i=0; i< num_nodes; i++){
      for (int j=0; j< num_nodes; j++){
	

	double x= gl_vals_x[id_xy(i,j)];
	double y= gl_vals_y[id_xy(i,j)];
     

	//Pressure
	luh[id_xyf(i,j,0)]  = 0*std::exp(-100*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));

	// //Velocity
	luh[id_xyf(i,j,1)]  = 0;
	luh[id_xyf(i,j,2)]  = 0;

	// Parameters
	luh[id_xyf(i,j,3)]  = 2.7;     // density [g/cm^3]
	luh[id_xyf(i,j,4)]  = 6.0;     // wavespeed [km/s]

	// Jacobian
	luh[id_xyf(i,j,5)] = jacobian[id_xy(i,j)];
	
	// metric derivatives
	luh[id_xyf(i,j,6)] = q_x[id_xy(i,j)];
	luh[id_xyf(i,j,7)] = q_y[id_xy(i,j)];
	luh[id_xyf(i,j,8)] = r_x[id_xy(i,j)];
	luh[id_xyf(i,j,9)] = r_y[id_xy(i,j)];
	
	// curvilinear mesh
	luh[id_xyf(i,j,10)] = gl_vals_x[id_xy(i,j)];
	luh[id_xyf(i,j,11)] = gl_vals_y[id_xy(i,j)];

	//std::cout<<gl_vals_x[id_xy(i,j)]<< " " << gl_vals_y[id_xy(i,j)]<<std::endl;


      }
    }

  }
}

void Linear::MyLinearWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 9

  // @todo Please implement/augment if required
  constexpr int numberOfVariables  = MyLinearWaveSolver::NumberOfVariables;
  constexpr int numberOfParameters = MyLinearWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
 

  for (int i = 0; i<numberOfData; i++){
    stateOut[i] = stateIn[i];
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
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 9
  
  // @todo Please implement/augment if required

  double cp = Q[4];

  double q_x=Q[6];
  double q_y=Q[7];
  double r_x=Q[8];
  double r_y=Q[9];
  
  lambda[0] = std::sqrt(q_x*q_x + q_y*q_y)*cp;
  lambda[1] = std::sqrt(r_x*r_x + r_y*r_y)*cp;
  lambda[2] =  0;

  //std::cout<< lambda[0]<< " " <<  lambda[1]<< " " <<  lambda[2]<<std::endl;
  
}


void Linear::MyLinearWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 3 + 9
  
  // @todo Please implement/augment if required

  double jacobian=Q[5];
  
  double q_x=Q[6];
  double q_y=Q[7];
  double r_x=Q[8];
  double r_y=Q[9];

  double u=Q[1];
  double v=Q[2];
  
  F[0][0] = -jacobian*(q_x*u+q_y*v);
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  

  F[1][0] = -jacobian*(r_x*u+r_y*v);  
  F[1][1] = 0.0;
  F[1][2] = 0.0;
   
}



void  Linear::MyLinearWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required

  double jacobian=Q[5];
  
  double q_x=Q[6];
  double q_y=Q[7];
  double r_x=Q[8];
  double r_y=Q[9];

  double p_q = gradQ[0];
  double u_q = gradQ[1];
  double v_q = gradQ[2];

  double p_r = gradQ[3];
  double u_r = gradQ[4];
  double v_r = gradQ[5];

  //  BgradQ[0] = -lam*(u_q*q_x+v_q*q_y);
  BgradQ[0] = 0;
  BgradQ[1] = -q_x*p_q;
  BgradQ[2] = -q_y*p_q ;

  //  BgradQ[3] = -lam*(u_r*r_x+v_r*r_y);
  BgradQ[3] = 0;
  BgradQ[4] = -r_x*p_r;
  BgradQ[5] = -r_y*p_r;
  
  
}

void Linear::MyLinearWaveSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn){
  static tarch::logging::Log _log("MyLinearWaveSolver::coefficientMatrix");

  double cp = Q[4];
  double rho = Q[3];

  //double cp = 1.0;
  //double rho = 1.0;
  double lam = rho*cp*cp;

  double nv[2] = {0.0};

  nv[d] = 1.0;
  
  double B1[3][3];
  double B2[3][3];
   
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      B1[i][j] = 0.0;
      B2[i][j] = 0.0;
    }
  }
  
  B1[0][1] = -lam; 
  B1[1][0] = -1/rho;
  
  B2[0][2] = -lam; 
  B2[2][0] = -1/rho;
  
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      Bn[i*3+ j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
    }
  }

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
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);
  double a_x = 0.0;
  double a_y = 0.0;
  double b_x = 10.0;
  double b_y = 5.0;
 
  double blockWidth_y = (b_y-a_y);
  double blockWidth_x = (b_x-a_x);

  if(n == 0){
    
    double x1 = 5.0;
    double y1 = 2.5;
     
     
    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);

    //std::cout << fault_ref_x << " " << fault_ref_y << std::endl;

    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    
    forceVector[0] = f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    
  }else if(n == 1){


    double x1 = 5.0;
    double y1 = 2.5;
     
     
    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);

    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    
    forceVector[0] = 0.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;

  }    
}

    /**
     * @TODO LR : document
     */
void Linear::MyLinearWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  // @todo Please implement/augment if required

  double rho = Q[3];
  double cp = Q[4];
  double jacobian= Q[5];

  double lam=cp*cp*rho;

  

  rhs[0] = lam/jacobian *rhs[0];
  rhs[1] = 1   /rho*rhs[1];
  rhs[2] = 1   /rho*rhs[2];

  rhs[3] = lam /jacobian *rhs[3];
  rhs[4] = 1   /rho*rhs[4];
  rhs[5] = 1   /rho*rhs[5];
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
  
  kernels::idx2 idx_QLR(basisSize, numberOfData);
  
  kernels::idx2 idx_FLR(basisSize, numberOfVariables);
  
  double n[2]={0,0};
  n[normalNonZeroIndex]=1;
  
  double cp_p;
  double rho_p;
  double lam_p;

  double cp_m;
  double rho_m;
  double lam_m;
  
  double n_p[2]={0,0};
  double n_m[2]={0,0};

  double m_p[2]={0,0};
  double m_m[2]={0,0};
  
  double norm_p_qr;
  double norm_m_qr;
  
  //std::cout<<" riemann solver called"<<std::endl;
  //std::exit(-1);
  

  // for (int k = 0; k < 2; k++){
    
  //   n_m[k] = n[k];
  //   n_p[k] = n[k];
    
  // }
  
  // norm_m_qr = 1.0;
  // norm_p_qr = 1.0;   
  
  
  for (int i = 0; i < basisSize; i++) {
    
    rho_m = QL[idx_QLR(i,3)];
    cp_m  = QL[idx_QLR(i,4)];
    lam_m = cp_m*cp_m*rho_m;

    rho_p = QR[idx_QLR(i,3)];
    cp_p  = QR[idx_QLR(i,4)];
    lam_p = cp_p*cp_p*rho_p;


    double qm_x=QL[idx_QLR(i,6)];
    double qm_y=QL[idx_QLR(i,7)];
    double rm_x=QL[idx_QLR(i,8)];
    double rm_y=QL[idx_QLR(i,9)];
    
    double qp_x=QR[idx_QLR(i,6)];
    double qp_y=QR[idx_QLR(i,7)];
    double rp_x=QR[idx_QLR(i,8)];
    double rp_y=QR[idx_QLR(i,9)];

    if (normalNonZeroIndex == 0){
	
      norm_m_qr = std::sqrt(qm_x*qm_x + qm_y*qm_y);
      n_m[0] = qm_x/norm_m_qr;
      n_m[1] = qm_y/norm_m_qr;
      
      norm_p_qr = std::sqrt(qp_x*qp_x + qp_y*qp_y);
      n_p[0] = qp_x/norm_p_qr;
      n_p[1] = qp_y/norm_p_qr;
      
      m_m[0] = n_m[1];
      m_m[1] =-n_m[0];
      
      m_p[0] = n_p[1];
      m_p[1] =-n_p[0];
      
    }

    if (normalNonZeroIndex == 1){
      
      norm_m_qr = std::sqrt(rm_x*rm_x + rm_y*rm_y);
      n_m[0] = rm_x/norm_m_qr;
      n_m[1] = rm_y/norm_m_qr;
      
      norm_p_qr = std::sqrt(rp_x*rp_x + rp_y*rp_y);
      n_p[0] = rp_x/norm_p_qr;
      n_p[1] = rp_y/norm_p_qr;
      
      m_m[0] = n_m[1];
      m_m[1] =-n_m[0];
      
      m_p[0] = n_p[1];
      m_p[1] =-n_p[0];
      
    }
    
    
    
    double v_m=QL[idx_QLR(i,1)]*n_m[0]+QL[idx_QLR(i,2)]*n_m[1];
    double v_p=QR[idx_QLR(i,1)]*n_p[0]+QR[idx_QLR(i,2)]*n_p[1];
    
    double sigma_m = QL[idx_QLR(i,0)];
    double sigma_p = QR[idx_QLR(i,0)];
    
    
    double z_p=rho_p*cp_p;
    double z_m=rho_m*cp_m;
    
    double v_hat_p=0;
    double v_hat_m=0;
    double sigma_hat_p=0;
    double sigma_hat_m=0;

    double r = 0.;

    if (faceIndex  == 3) { r = 1.0;} 
    
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
    else {// inter-element boundaries
      riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);
    }
        
    FR[idx_FLR(i,0)] = -norm_p_qr*0.5*lam_p*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p);
    FL[idx_FLR(i,0)] =  norm_m_qr*0.5*lam_m*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m);
    
    FR[idx_FLR(i,1)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[0];
    FL[idx_FLR(i,1)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[0];
    
    FR[idx_FLR(i,2)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[1];
    FL[idx_FLR(i,2)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[1];
        
  }
  
}
