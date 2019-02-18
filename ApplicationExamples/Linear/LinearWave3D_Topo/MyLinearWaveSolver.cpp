#include "MyLinearWaveSolver.h"

#include "MyLinearWaveSolver_Variables.h"

#include <algorithm>

#include "CurvilinearTransformation.h"


tarch::logging::Log Linear::MyLinearWaveSolver::_log( "Linear::MyLinearWaveSolver" );


void Linear::MyLinearWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Linear::MyLinearWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 15
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    //initialise luh

    constexpr int basisSize = MyLinearWaveSolver::Order+1;
  int numberOfData=MyLinearWaveSolver::NumberOfParameters+MyLinearWaveSolver::NumberOfVariables;

  kernels::idx4 id_xyzf(basisSize,basisSize,basisSize,numberOfData);
  kernels::idx3 id_xyz(basisSize,basisSize,basisSize);

  int num_nodes = basisSize;


 
  
  int ne_x = std::round(1/dx[0]);
  int ne_y = std::round(1/dx[1]);
  int ne_z = std::round(1/dx[2]);			    

  int nx = ne_x *(num_nodes-1) + 1;
  int ny = ne_y *(num_nodes-1) + 1;
  int nz = ne_z *(num_nodes-1) + 1;


  double offset_x=center[0]-0.5*dx[0];
  double offset_y=center[1]-0.5*dx[1];
  double offset_z=center[2]-0.5*dx[2];  

  double width_x=dx[0];
  double width_y=dx[1];
  double width_z=dx[2];


  int i_m;
  int j_m;
  int k_m;

  double a_x = 0;
  double a_y = 0;
  double a_z = 0;
  
  double b_x = 10.0;
  double b_y = 10.0;
  double b_z = 10.0;
  
  double fault_position = 5.0;
  
  double blockWidth_x;
  double blockWidth_y = b_y-a_y;
  double blockWidth_z = b_z-a_z;

  double block_width0_x = b_x-a_x;

  double fault_ref = (fault_position-a_x)/(block_width0_x);

  int n = center[0] > fault_ref ? 1 : 0 ;
  
  if(n == 0){

    blockWidth_x=(fault_position-a_x);
    ne_x = std::round((ne_x+1)*fault_ref);
    nx =  ne_x *(num_nodes-1)+1;
    

    i_m =  std::round((offset_x)/width_x) *(num_nodes-1);
    j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    k_m =  std::round((offset_z)/width_z) *(num_nodes-1);
    
  }else{
    
    double ne_x0 = std::round((ne_x+1)*fault_ref);
    
    ne_x = std::round(1/dx[0])-ne_x0;
    nx =  ne_x *(num_nodes-1)+1;
    blockWidth_x= (b_x-fault_position);


    i_m =  std::floor((offset_x-fault_ref)/width_x) *(num_nodes-1);
    j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    k_m =  std::round((offset_z)/width_z) *(num_nodes-1);    
  }

  double* left_bnd_x = new double[ny*nz];
  double* left_bnd_y = new double[ny*nz];
  double* left_bnd_z = new double[ny*nz];  

  double* right_bnd_x = new double[ny*nz];
  double* right_bnd_y = new double[ny*nz];
  double* right_bnd_z = new double[ny*nz];  

  double* bottom_bnd_x = new double[nx*nz];
  double* bottom_bnd_y = new double[nx*nz];
  double* bottom_bnd_z = new double[nx*nz];  

  double* top_bnd_x = new double[nx*nz];
  double* top_bnd_y = new double[nx*nz];
  double* top_bnd_z = new double[nx*nz];  

  double* front_bnd_x = new double[nx*ny];
  double* front_bnd_y = new double[nx*ny];
  double* front_bnd_z = new double[nx*ny];  

  double* back_bnd_x = new double[nx*ny];
  double* back_bnd_y = new double[nx*ny];
  double* back_bnd_z = new double[nx*ny];  
  
  
  // getBoundaryCurves3D( num_nodes,
  // 		       offset_x,  offset_y,  offset_z,
  // 		       width_x,  width_y ,  width_z ,
  // 		       left_bnd_x,  left_bnd_y,  left_bnd_z,
  // 		       right_bnd_x,  right_bnd_y,  right_bnd_z,
  // 		       bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
  // 		       top_bnd_x,  top_bnd_y,  top_bnd_z,
  // 		       front_bnd_x,  front_bnd_y,  front_bnd_z,
  // 		       back_bnd_x,  back_bnd_y,  back_bnd_z);

  

  
  // getBoundaryCurves3D_fixedTopFace_forBlock( num_nodes,
  // 					     nx,ny,nz,n	,	     
  // 					     block_width_x,  block_width_y ,  block_width_z ,
  // 					     left_bnd_x,  left_bnd_y,  left_bnd_z,
  // 					     right_bnd_x,  right_bnd_y,  right_bnd_z,
  // 					     bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
  // 					     top_bnd_x,  top_bnd_y,  top_bnd_z,
  // 					     front_bnd_x,  front_bnd_y,  front_bnd_z,
  // 					     back_bnd_x,  back_bnd_y,  back_bnd_z);

  getBoundaryCurves3D_cutOffTopography_withFault( num_nodes,
						  nx,ny,nz,n,fault_position,
						  a_x, a_y, a_z,
						  b_x, b_y, b_z,
						  blockWidth_x,  blockWidth_y,  blockWidth_z,
						  left_bnd_x,  left_bnd_y,  left_bnd_z,
						  right_bnd_x,  right_bnd_y,  right_bnd_z,
						  bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
						  top_bnd_x,  top_bnd_y,  top_bnd_z,
						  front_bnd_x,  front_bnd_y,  front_bnd_z,
						  back_bnd_x,  back_bnd_y,  back_bnd_z);


   
  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); //left right


  double* curvilinear_x = new double[num_nodes*num_nodes*num_nodes];
  double* curvilinear_y = new double[num_nodes*num_nodes*num_nodes];
  double* curvilinear_z = new double[num_nodes*num_nodes*num_nodes];  

 	
  int i_p = i_m + num_nodes;
  int j_p = j_m + num_nodes;
  int k_p = k_m + num_nodes;   

  
  transFiniteInterpolation3D( nx,  ny,  nz,
  			      k_m,  k_p ,
  			      j_m,  j_p ,
  			      i_m,  i_p ,
  			      num_nodes,
			      width_x,width_y,width_z,
  			      left_bnd_x,
  			      right_bnd_x,
  			      bottom_bnd_x,
  			      top_bnd_x,
  			      front_bnd_x,
  			      back_bnd_x,
  			      curvilinear_x
  			      );

  transFiniteInterpolation3D( nx,  ny,  nz,
  			      k_m,  k_p ,
  			      j_m,  j_p ,
  			      i_m,  i_p ,
  			      num_nodes,
			      width_x,width_y,width_z,
  			      left_bnd_y,
  			      right_bnd_y,
  			      bottom_bnd_y,
  			      top_bnd_y,
  			      front_bnd_y,
  			      back_bnd_y,
  			      curvilinear_y
  			      );

  
  transFiniteInterpolation3D( nx,  ny,  nz,
  			      k_m,  k_p ,
  			      j_m,  j_p ,
  			      i_m,  i_p ,
  			      num_nodes,
			      width_x,width_y,width_z,			      
  			      left_bnd_z,
  			      right_bnd_z,
  			      bottom_bnd_z,
  			      top_bnd_z,
  			      front_bnd_z,
  			      back_bnd_z,
  			      curvilinear_z
  			      );

  
  double* gl_vals_x = new double[num_nodes*num_nodes*num_nodes];
  double* gl_vals_y = new double[num_nodes*num_nodes*num_nodes];
  double* gl_vals_z = new double[num_nodes*num_nodes*num_nodes];

  double* jacobian = new double[num_nodes*num_nodes*num_nodes];

  double* q_x = new double[num_nodes*num_nodes*num_nodes];
  double* q_y = new double[num_nodes*num_nodes*num_nodes];
  double* q_z = new double[num_nodes*num_nodes*num_nodes];
  
  double* r_x = new double[num_nodes*num_nodes*num_nodes];
  double* r_y = new double[num_nodes*num_nodes*num_nodes];
  double* r_z = new double[num_nodes*num_nodes*num_nodes];

  double* s_x = new double[num_nodes*num_nodes*num_nodes];
  double* s_y = new double[num_nodes*num_nodes*num_nodes];
  double* s_z = new double[num_nodes*num_nodes*num_nodes];  

  
  metricDerivativesAndJacobian3D(num_nodes,
  				 curvilinear_x,  curvilinear_y,  curvilinear_z,
  				 gl_vals_x,  gl_vals_y,  gl_vals_z,
  				 q_x,  q_y,  q_z,
  				 r_x,  r_y,  r_z,
  				 s_x,  s_y,  s_z,				  
  				 jacobian,
  				 width_x,  width_y,  width_z
  				 );


  for (int k=0; k< num_nodes; k++){
    for (int j=0; j< num_nodes; j++){
      for (int i=0; i< num_nodes; i++){


	
	double x= gl_vals_x[id_xyz(k,j,i)];
	double y= gl_vals_y[id_xyz(k,j,i)];
	double z= gl_vals_z[id_xyz(k,j,i)];

	// pressure
	luh[id_xyzf(k,j,i,0)]  = std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	
	// velocity
	luh[id_xyzf(k,j,i,1)]  = 0;
	luh[id_xyzf(k,j,i,2)]  = 0;
	luh[id_xyzf(k,j,i,3)]  = 0;

	// material parameters
	luh[id_xyzf(k,j,i,4)]  = 2.7;    //rho [g/cm^3]
	luh[id_xyzf(k,j,i,5)]  = 6.0;  //c   [km/s]	

	// Jacobian
	luh[id_xyzf(k,j,i,6)]  = jacobian[id_xyz(k,j,i)];

	// metric derivatives:
	// q_x, q_y, q_z
	luh[id_xyzf(k,j,i,7)]  = q_x[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,8)]  = q_y[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,9)]  = q_z[id_xyz(k,j,i)];

	// r_x, r_y, r_z
	luh[id_xyzf(k,j,i,10)] = r_x[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,11)] = r_y[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,12)] = r_z[id_xyz(k,j,i)];

	// s_x, s_y, s_z
	luh[id_xyzf(k,j,i,13)] = s_x[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,14)] = s_y[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,15)] = s_z[id_xyz(k,j,i)];

	//
	luh[id_xyzf(k,j,i,16)] = gl_vals_x[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,17)] = gl_vals_y[id_xyz(k,j,i)];
	luh[id_xyzf(k,j,i,18)] = gl_vals_z[id_xyz(k,j,i)];

      }
    }
  }
  }
}

void Linear::MyLinearWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 15

  // @todo Please implement/augment if required
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
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 15
  
  // @todo Please implement/augment if required
  double c = Q[5];

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];
  
  lambda[ 0] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*c;
  lambda[ 1] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*c;
  lambda[ 2] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*c;
  lambda[ 3] = 0;
}


void Linear::MyLinearWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 4 + 15
  
  // @todo Please implement/augment if required
  double u = Q[1];
  double v = Q[2];
  double w = Q[3];


  double jacobian = Q[6];

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];

  
  
  F[0][ 0] = -jacobian*(q_x*u+q_y*v+q_z*w);
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;

  F[1][ 0] = -jacobian*(r_x*u+r_y*v+r_z*w);
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;

  F[2][ 0] = -jacobian*(s_x*u+s_y*v+s_z*w);
  F[2][ 1] = 0.0;
  F[2][ 2] = 0.0;
  F[2][ 3] = 0.0;
 
}



void  Linear::MyLinearWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  double p_q = gradQ[0];  
  double u_q = gradQ[1];
  double v_q = gradQ[2];
  double w_q = gradQ[3];

  double p_r = gradQ[4];  
  double u_r = gradQ[5];
  double v_r = gradQ[6];
  double w_r = gradQ[7];

  double p_s = gradQ[8];  
  double u_s = gradQ[9];
  double v_s = gradQ[10];
  double w_s = gradQ[11];


  double jacobian = Q[6];

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];
  
  BgradQ[0]= 0;  
  BgradQ[1]= -q_x*p_q;
  BgradQ[2]= -q_y*p_q;
  BgradQ[3]= -q_z*p_q;

  BgradQ[4]=0;  
  BgradQ[5]=-r_x*p_r;
  BgradQ[6]=-r_y*p_r;
  BgradQ[7]=-r_z*p_r;

  BgradQ[8]=0;  
  BgradQ[9] =-s_x*p_s;
  BgradQ[10]=-s_y*p_s;
  BgradQ[11]=-s_z*p_s;
  
}

void Linear::MyLinearWaveSolver::coefficientMatrix(const double* const Q,const int d,double* const Bn){
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

void  Linear::MyLinearWaveSolver::pointSource(const double* const x,const double t,const double dt, double* const forceVector, double* const x0, int n) {
  // @todo Please implement/augment if required

  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;

  double a_x = 0.0;
  double a_y = 0.0;
  double a_z = 0.0;
  
  double b_x = 10.0;
  double b_y = 10.0;
  double b_z = 10.0;
 
  double blockWidth_y = (b_y-a_y);
  double blockWidth_x = (b_x-a_x);
  double blockWidth_z = (b_z-a_z);

  if(n == 0){
    
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

    double x1 = 5.0;
    double y1 = 5.0;
    double z1 = 5.0;

    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);
    double fault_ref_z = (z1-a_z)/(blockWidth_z);
    
    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    x0[2] = fault_ref_z;
    
    forceVector[0] = 1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
    
  }else if(n == 1){
    
    f = 0*M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

    double x1 = 7.5;
    double y1 = 5.0;
    double z1 = 5.0;

    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);
    double fault_ref_z = (z1-a_z)/(blockWidth_z);
    
    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    x0[2] = fault_ref_z;
    
    forceVector[0] = 1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
    
  }else if(n == 2){
    
    f = 0*M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

    double x1 = 5.0;
    double y1 = 2.5;
    double z1 = 5.0;

    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);
    double fault_ref_z = (z1-a_z)/(blockWidth_z);
    
    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    x0[2] = fault_ref_z;
    
    
    forceVector[0] = 1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
  }else if(n == 3){
    
    f = 0*M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

    double x1 = 5.0;
    double y1 = 7.5;
    double z1 = 5.0;

    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);
    double fault_ref_z = (z1-a_z)/(blockWidth_z);
    
    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    x0[2] = fault_ref_z;
    
    forceVector[0] = 1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
  }
}

    /**
     * @TODO LR : document
     */
void Linear::MyLinearWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  // @todo Please implement/augment if required
  double rho = Q[4];  
  double c   = Q[5];
  double jacobian = Q[6];  
  double mu  = rho*c*c;

  rhs[0]=mu/jacobian * rhs[0];
  rhs[1]=1/rho * rhs[1];
  rhs[2]=1/rho * rhs[2];
  rhs[3]=1/rho * rhs[3];

  rhs[4]=mu/jacobian * rhs[4];
  rhs[5]=1/rho * rhs[5];
  rhs[6]=1/rho * rhs[6];
  rhs[7]=1/rho * rhs[7];

  rhs[8] =mu/jacobian * rhs[8];
  rhs[9] =1/rho * rhs[9];
  rhs[10]=1/rho * rhs[10];
  rhs[11]=1/rho * rhs[11];
}

void Linear::MyLinearWaveSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex){
  
  constexpr int numberOfVariables  = MyLinearWaveSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyLinearWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyLinearWaveSolver::Order+1;
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
  
  
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      
      double c_p=QR[idx_QLR(i,j,5)];
      double rho_p=QR[idx_QLR(i,j,4)];
      double mu_p=c_p*c_p*rho_p;
      
      double c_m=QL[idx_QLR(i,j,5)];
      double rho_m=QL[idx_QLR(i,j,4)];
      double mu_m=c_m*c_m*rho_m;
      
      get_normals(normalNonZeroIndex, norm_p_qr, n_p, QR + idx_QLR(i,j,0));
      get_normals(normalNonZeroIndex, norm_m_qr, n_m, QL + idx_QLR(i,j,0));  
      
      localBasis(n_p, m_p, l_p, 3);
      localBasis(n_m, m_m, l_m, 3);     
      
      double v_m=QL[idx_QLR(i,j,1)]*n_m[0]+QL[idx_QLR(i,j,2)]*n_m[1]+QL[idx_QLR(i,j,3)]*n_m[2];
      double v_p=QR[idx_QLR(i,j,1)]*n_p[0]+QR[idx_QLR(i,j,2)]*n_p[1]+QR[idx_QLR(i,j,3)]*n_p[2];
      
      double sigma_m = QL[idx_QLR(i,j,0)];
      double sigma_p = QR[idx_QLR(i,j,0)];
      
      double z_p=rho_p*c_p;
      double z_m=rho_m*c_m;
      
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
      else { // inter-element boundaries
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
      
    }
  }
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

 void Linear::MyLinearWaveSolver::Gram_Schmidt(double* const y, double* const z){
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

void Linear::MyLinearWaveSolver::localBasis(double* const n, double* const m, double* const l, int d){

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

void Linear::MyLinearWaveSolver::get_normals(int normalNonZeroIndex,double& norm, double* const n,const double* const Q)

{

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];

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
