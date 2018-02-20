#include "MyLinearSolver.h"
#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"

#include "MyLinearSolver_Variables.h"

#include <algorithm>

#include "CurvilinearTransformation.h"


tarch::logging::Log Linear::MyLinearSolver::_log( "Linear::MyLinearSolver" );


void Linear::MyLinearSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  static tarch::logging::Log _log("MyLinearSolver::init");
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Linear::MyLinearSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::AdjustSolutionValue");
  //  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PatchWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;  
}

void Linear::MyLinearSolver::adjustPatchSolution(
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt,
      double* luh) {

  constexpr int basisSize = MyLinearSolver::Order+1;
  int num_nodes = basisSize;
  int numberOfData=MyLinearSolver::NumberOfParameters+MyLinearSolver::NumberOfVariables;

  kernels::idx3 id_xyz(basisSize,basisSize,numberOfData);
  kernels::idx2 id_xy(basisSize,basisSize);  

  double* left_bnd_x = new double[num_nodes];

  double* left_bnd_y = new double[num_nodes];

  double* right_bnd_x = new double[num_nodes];
  double* right_bnd_y = new double[num_nodes];

  double* bottom_bnd_x = new double[num_nodes];
  double* bottom_bnd_y = new double[num_nodes];

  double* top_bnd_x = new double[num_nodes];
  double* top_bnd_y = new double[num_nodes];

  getBoundaryCurves(num_nodes,left_bnd_x,left_bnd_y,right_bnd_x,right_bnd_y,bottom_bnd_x,bottom_bnd_y,top_bnd_x,top_bnd_y);


  double* curvilinear_x = new double[num_nodes*num_nodes];
  double* curvilinear_y = new double[num_nodes*num_nodes];


  transFiniteInterpolation( num_nodes,  left_bnd_x,  left_bnd_y,  right_bnd_x,  right_bnd_y,  bottom_bnd_x,  bottom_bnd_y,  top_bnd_x,  top_bnd_y,  curvilinear_x ,  curvilinear_y );

  
  double* gl_vals_x = new double[num_nodes*num_nodes];
  double* gl_vals_y = new double[num_nodes*num_nodes];
  double* jacobian = new double[num_nodes*num_nodes];

  double* q_x = new double[num_nodes*num_nodes];
  double* q_y = new double[num_nodes*num_nodes];
  
  double* r_x = new double[num_nodes*num_nodes];
  double* r_y = new double[num_nodes*num_nodes];  
  
  metricDerivativesAndJacobian(num_nodes,curvilinear_x,curvilinear_y,gl_vals_x,gl_vals_y,q_x,q_y,r_x,r_y,jacobian);
  

  for (int i=0; i< num_nodes; i++){
      for (int j=0; j< num_nodes; j++){
	double x= gl_vals_x[id_xy(j,i)];
	double y= gl_vals_y[id_xy(j,i)];

	
	//Pressure
	luh[id_xyz(i,j,0)] = 10*std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01);
	//luh[id_xyz(j,i,0)] = x*0.1;
	
	//Velocity
	luh[id_xyz(j,i,1)] = 0;
	luh[id_xyz(j,i,2)] = 0;
	
	  
	luh[id_xyz(j,i,3)] = jacobian[id_xy(j,i)];

	luh[id_xyz(j,i,4)] = q_x[id_xy(j,i)];
	luh[id_xyz(j,i,5)] = q_y[id_xy(j,i)];
	luh[id_xyz(j,i,6)] = r_x[id_xy(j,i)];
	luh[id_xyz(j,i,7)] = r_y[id_xy(j,i)];	
	
	luh[id_xyz(j,i,8)] = gl_vals_x[id_xy(j,i)];
	luh[id_xyz(j,i,9)] = gl_vals_y[id_xy(j,i)];	

      }
  }


  
}

void Linear::MyLinearSolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  static tarch::logging::Log _log("MyLinearSolver::adjustPointSolution");

  Variables vars(Q);

  //  vars.p() = std::exp(-((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5))/0.01);
  vars.p() = 1;
  //  vars.p()=0;
  //  vars.v(0,0);
  vars.v(1,1);  


 //  kernels::idx2 id_xy(num_nodes,num_nodes);

 //  for(int j = 0 ; j < num_nodes ; j ++){
 //    //    printf("%f \n",top_bnd_y[j]);
 //    for(int i = 0 ; i< num_nodes ; i ++){

 //      // printf("%d \n",i);
 //      // printf("%d \n",j);

 //      // printf("%f\n",curvilinear_x[id_xy(j,i)]);
 //      //       printf("%f\n",curvilinear_y[id_xy(i,j)]);
 //    } 
 //  } 

 //  for(int j = 0 ; j < num_nodes ; j ++){
 //    for(int i = 0 ; i< num_nodes ; i ++){

 //      // printf("%d \n",i);
 //      // printf("%d \n",j);

 //      // printf("%f\n",curvilinear_x[id_xy(j,i)]);
 //    } 
 //  } 

  
 // double* gl_vals_x = new double[num_nodes*num_nodes];
 // double* gl_vals_y = new double[num_nodes*num_nodes];
 // double* jacobian = new double[num_nodes*num_nodes];
 // double* q_x = new double[num_nodes*num_nodes];
 // double* q_y = new double[num_nodes*num_nodes];

 // double* r_x = new double[num_nodes*num_nodes];
 // double* r_y = new double[num_nodes*num_nodes];  
 
 // metricDerivativesAndJacobian(num_nodes,curvilinear_x,curvilinear_y,gl_vals_x,gl_vals_y,q_x,q_y,r_x,r_y,jacobian);

 
 // for(int j = 0 ; j < num_nodes ; j ++){
 //    for(int i = 0 ; i< num_nodes ; i ++){
 //      printf("\n");
 //      // printf("%f\n",jacobian[id_xy(i,j)]);
 //      printf("%d\n",(0<jacobian[id_xy(i,j)]));      
 //      //printf("%f\n",q_x[id_xy(i,j)]);
 //      //printf("%f\n",q_y[id_xy(i,j)]);
 //      //printf("%f\n",r_x[id_xy(i,j)]);
 //      //printf("%f\n",r_y[id_xy(i,j)]);      
 //    }

 // } 

 

 // std::exit(-1);
}

void Linear::MyLinearSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters
  
  // @todo Please implement/augment if required
  
  double cp = 6.0;

  static tarch::logging::Log _log("MyLinearSolver::eigenvalues");
  lambda[0] =  cp;
  lambda[1] = -cp;
  lambda[2] =  0;

}


void Linear::MyLinearSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters
  
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::flux");

  //  std::cout << "called Flux" << std::endl;
  double q_x,q_y,r_x,r_y,jacobian,p_q,p_r;

  //  double cp=6.0;
  //  double rho=2.7;
  //  double lam=cp*cp*rho;      

  jacobian=Q[3];
  q_x=Q[4];
  q_y=Q[5];
  r_x=Q[6];
  r_y=Q[7];
  
  double u=Q[1];
  double v=Q[2];
  
  F[0][0] = -jacobian*(q_x*u+q_y*v);
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  

  F[1][0] = -jacobian*(r_x*u+r_y*v);  
  F[1][1] = 0.0;
  F[1][2] = 0.0;
 
  

}

void Linear::MyLinearSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
					      const double * const fluxIn,const double* const stateIn,
					      double *fluxOut,double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 3 + #parameters

  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::boundaryValues");

  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  fluxOut[0] =  fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];

  
  double cp=6.0;
  double rho=2.7;
  double lam=cp*cp*rho;      

  double z = rho*cp;

  double q_x=stateIn[4];
  double q_y=stateIn[5];
  double r_x=stateIn[6];
  double r_y=stateIn[7];

  double v_x =  stateIn[1];
  double v_y =  stateIn[2];

  double norm_q=sqrt(q_x*q_x+q_y*q_y);
  double norm_r=sqrt(r_x*r_x+r_y*r_y);

  double n_x;
  double n_y;  

  if (faceIndex == 0) {
      
    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
   
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
 
    n_x=q_x/norm_q;
    n_y=q_y/norm_q;

    double v_n = v_x*n_x + v_y *n_y;    

    double p  =  stateIn[0];


    
    double vx_hat =  0.;
    double p_hat =  0.;
   
    
    double r = 0.;
    
    riemannSolver_BC0(v_n, p, z, r, vx_hat, p_hat);
   
    
    stateOut[0] = p_hat;
    stateOut[1] = vx_hat;
    
  }
  
  
  if (faceIndex == 1) {

    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    
    
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
 
    n_x=q_x/norm_q;
    n_y=q_y/norm_q;
    double v_n = v_x*n_x + v_y *n_y;    
    

    double p  =  stateIn[0];
    
    
    double vx_hat =  0.;
    double p_hat =  0.;
   
    
    double r = 0.;
    
    riemannSolver_BCn(v_n, p, z, r, vx_hat, p_hat);
   
    stateOut[0] = p_hat;
    stateOut[1] = vx_hat;

      
  }


  if (faceIndex == 2) {
    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    
    
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
 
    n_x=r_x/norm_r;
    n_y=r_y/norm_r;
    double v_n = v_x*n_x + v_y *n_y;    


    double p  =  stateIn[0];
    
    
    
    double vy_hat =  0.;
    
    double p_hat =  0.;
   
    
    double r = 0.;
    
    riemannSolver_BC0(v_n, p, z, r, vy_hat, p_hat);
   
    
    stateOut[0] = p_hat;
    stateOut[2] = vy_hat;

   
  }

  if (faceIndex == 3) {
    
    stateOut[0] = stateIn[0];
    stateOut[1] = stateIn[1];
    stateOut[2] = stateIn[2];
    
    
    
    fluxOut[0] = fluxIn[0];
    fluxOut[1] = fluxIn[1];
    fluxOut[2] = fluxIn[2];
 
    n_x=r_x/norm_r;
    n_y=r_y/norm_r;
    double v_n = v_x*n_x + v_y *n_y;    


    double p  =  stateIn[0];
    
    
    double vy_hat =  0.;
    double p_hat =  0.;
   
    
    double r = 1.;
    
    riemannSolver_BCn(v_n, p, z, r, vy_hat, p_hat);
   
    
    stateOut[0] = p_hat;
    stateOut[2] = vy_hat;
  }

  //  std::cout << faceIndex << std::endl;
  //  std::cout << n_x << std::endl;
  //  std::cout << n_y << std::endl;
  //  std::cout << std::endl;    
}


exahype::solvers::Solver::RefinementControl Linear::MyLinearSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  static tarch::logging::Log _log("MyLinearSolver::refinementCriterion");
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Linear::MyLinearSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ){
  kernels::idx2 idx(DIMENSIONS, NumberOfVariables);

  static tarch::logging::Log _log("MyLinearSolver::nonConservativeProduct");

  double q_x,q_y,r_x,r_y,jacobian;

  //double cp=6.0;
  //double rho=2.7;
  //double lam=cp*cp*rho;      
  
  jacobian=Q[3];
  q_x=Q[4];
  q_y=Q[5];
  r_x=Q[6];
  r_y=Q[7];
  
  double p_q=gradQ[0];
  double p_r=gradQ[3];

  double u_q=gradQ[1];
  double u_r=gradQ[4];
  
  double v_q=gradQ[2];  
  double v_r=gradQ[5];


  // std::cout << gradQ[0] << std::endl;
  // std::cout << gradQ[3] << std::endl;
  // std::cout << gradQ[1] << std::endl;
  // std::cout << gradQ[5] << std::endl;

  //BgradQ[0] = 0;
  
  //  std::cout << "metic coefficients" << std::endl;
  //std::cout << q_x << std::endl;
  // std::cout << q_y << std::endl;
  // std::cout << r_x << std::endl;
  // std::cout << r_y << std::endl;
  // std::cout << jacobian << std::endl;    


  // BgradQ[0] = -lam*gradQ[1];
  // BgradQ[1] = -1/rho*gradQ[0];
  // BgradQ[2] = 0;
  
  // BgradQ[3]= -lam*gradQ[5];
  // BgradQ[4]= 0;
  // BgradQ[5]= -1/rho*gradQ[3];
  
  
  //  BgradQ[0] = -lam*(u_q*q_x+v_q*q_y);
  BgradQ[0] = 0;
  BgradQ[1] = -p_q*q_x;
  BgradQ[2] = -q_y*p_q ;

  //  BgradQ[3] = -lam*(u_r*r_x+v_r*r_y);
  BgradQ[3] = 0;
  BgradQ[4] = -p_r*r_x;
  BgradQ[5] = -p_r*r_y;
  
  //  BgradQ[0] = -gradQ[1];
  // BgradQ[0] = 0;
  // BgradQ[1] = q_x*p_q;
  // BgradQ[2] = q_y*p_q;
  // BgradQ[3] = 0;
  // BgradQ[4] = 0;
  // BgradQ[5] = 0;
  // BgradQ[6] = 0;
  // BgradQ[7] = 0;
  // BgradQ[8] = 0;
  // BgradQ[9] = 0;

 // //BgradQ[10]= -gradQ[12];
  // BgradQ[10]= 0;
  // BgradQ[11]= r_x*p_r;
  // BgradQ[12]= r_y*p_r;
  // BgradQ[13] = 0;
  // BgradQ[14] = 0;
  // BgradQ[15] = 0;
  // BgradQ[16] = 0;
  // BgradQ[17] = 0;
  // BgradQ[18] = 0;
  // BgradQ[19] = 0;

}


void Linear::MyLinearSolver::coefficientMatrix(const double* const Q,const int d,double* Bn){
  static tarch::logging::Log _log("MyLinearSolver::coefficientMatrix");

  double cp = 6.0;
  double rho = 2.7;
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


void Linear::MyLinearSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex,bool isBoundaryFace, int faceIndex){

  constexpr int numberOfVariables  = MyLinearSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyLinearSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyLinearSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx2 idx_QLR(basisSize, numberOfData);

  kernels::idx2 idx_FLR(basisSize, numberOfVariables);

  double n[2]={0,0};
  n[normalNonZeroIndex]=1;

  double cp=6.0;
  double rho=2.7;
  double lam=cp*cp*rho;      


    for (int i = 0; i < basisSize; i++) {
      double v_m=QL[idx_QLR(i,1)]*n[0]+QL[idx_QLR(i,2)]*n[1];
      double v_p=QR[idx_QLR(i,1)]*n[0]+QR[idx_QLR(i,2)]*n[1];

      double sigma_m = QL[idx_QLR(i,0)];
      double sigma_p = QR[idx_QLR(i,0)];


      double z_p=rho*cp;
      double z_m=rho*cp;

      double v_hat_p=0;
      double v_hat_m=0;
      double sigma_hat_p=0;
      double sigma_hat_m=0;

      riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);


      FR[idx_FLR(i, 0)] = -0.5*lam*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p);
      FL[idx_FLR(i, 0)] =  0.5*lam*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m);


      FR[idx_FLR(i, 1)] = 0.5/rho*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n[0];
      FL[idx_FLR(i, 1)] = 0.5/rho*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n[0];

      FR[idx_FLR(i, 2)] = 0.5/rho*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n[1];
      FL[idx_FLR(i, 2)] = 0.5/rho*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n[1];

    }

}



void Linear::MyLinearSolver::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0){
     double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

  x0[0] = 0.5;
  x0[1] = 0.5;

  forceVector[0] = 1.*f;
  forceVector[1] = 0.0;
  forceVector[2] = 0.0;
}



 void Linear::MyLinearSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double v_hat_p , double v_hat_m, double sigma_hat_p, double sigma_hat_m){
   double p=0;
   double q=0;
   double phi=0;
   double v_hat=0;
   double eta=0;

   p=z_p*v_m + sigma_m;
   q=z_m*v_p - sigma_p;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= p/z_p - q/z_m;

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


void Linear::MyLinearSolver::multiplyMaterialParameterMatrix(const double* const Q, double* rhs){

  double jacobian= Q[4];
  double cp=6.0;
  double rho=2.7;
  double lam=cp*cp*rho;      

  rhs[0] = lam/jacobian *rhs[0];
  rhs[1] = 1   /rho*rhs[1];
  rhs[2] = 1   /rho*rhs[2];

  rhs[3] = lam /jacobian *rhs[3];
  rhs[4] = 1   /rho*rhs[4];
  rhs[5] = 1   /rho*rhs[5];

}
