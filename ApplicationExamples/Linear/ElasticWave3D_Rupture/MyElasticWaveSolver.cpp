#include "MyElasticWaveSolver.h"

#include "MyElasticWaveSolver_Variables.h"

#include "CurvilinearTransformation.h"

tarch::logging::Log Elastic::MyElasticWaveSolver::_log( "Elastic::MyElasticWaveSolver" );

void Elastic::MyElasticWaveSolver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Elastic::MyElasticWaveSolver::adjustSolution(double* const luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 16
  // @todo Please implement/augment if required
  int level=std::round(log(_domainSize[0]/dx[0])/log(3.)) + 1;
  if (tarch::la::equals(t,0.0)) {
    
    constexpr int basisSize = MyElasticWaveSolver::Order+1;
    int numberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

    kernels::idx4 id_xyzf(basisSize,basisSize,basisSize,numberOfData);
    kernels::idx3 id_xyz(basisSize,basisSize,basisSize);

    int num_nodes = basisSize;

    double offset_x=center[0]-0.5*dx[0];
    double offset_y=center[1]-0.5*dx[1];
    double offset_z=center[2]-0.5*dx[2];  

    double width_x=dx[0];
    double width_y=dx[1];
    double width_z=dx[2];


    // //Number of elements on the unit square
    // int ne_x = std::round(1/dx[0]);
    // int ne_y = std::round(1/dx[1]);
    // int ne_z = std::round(1/dx[2]);			    

    // //Number of nodes on the unit square
    // int nx = ne_x *(num_nodes-1) + 1;
    // int ny = ne_y *(num_nodes-1) + 1;
    // int nz = ne_z *(num_nodes-1) + 1;

    //first local node indeces in global domain
    int i_m;
    int j_m;
    int k_m;

    // Coordiantes of the undeformed domain (cuboidal)
    double a_x = 0;
    double a_y = 0;
    double a_z = 0;
  
    double b_x = 10.0;
    double b_y = 10.0;
    double b_z = 10.0;

    int ne_x = std::round((b_x-a_x)/dx[0]); //number of elements in x direction
    int ne_y = std::round((b_y-a_y)/dx[1]); //number of elements in y direction
    int ne_z = std::round((b_z-a_z)/dx[2]); //number of elements in z direction

    //Number of nodes on the unit square
    int nx = ne_x *(num_nodes-1) + 1;
    int ny = ne_y *(num_nodes-1) + 1;
    int nz = ne_z *(num_nodes-1) + 1;

    // nx = ne_x *(num_nodes-1) + 1; //global number of nodes in x direction considering collocation
    // ny = ne_y *(num_nodes-1) + 1; //global number of nodes in y direction
    // nz = ne_z *(num_nodes-1) + 1; //global number of nodes in z direction
    
    // Position of interface in x  
    double fault_position = 5.0;

    // Width of the current block    
    double blockWidth_x;
    double blockWidth_y = b_y-a_y; //equal to width of the rectangle
    double blockWidth_z = b_z-a_z; //equal to width of the rectangle


    double recWidth_x = b_x-a_x; //width in x of the rectangle

    //relative position of the fault on the unit square
    double fault_ref = (fault_position-a_x)/(recWidth_x);
    
    // number of the current block
    //int n = center[0] > fault_ref ? 1 : 0 ;
    int n = offset_x >  fault_position ? 1 : 0 ; //1: Water Column 0: Solid

    if(n == 0){
    blockWidth_x=fault_position;
    ne_x = std::round((ne_x+1)*fault_position/recWidth_x);
    nx =  ne_x *(num_nodes-1)+1;


    i_m =  std::round((offset_x)/width_x) *(num_nodes-1);
    j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    k_m =  std::round((offset_z)/width_z) *(num_nodes-1);    
  }else{
    blockWidth_x=b_x-fault_position;
    ne_x = ne_x-std::round((ne_x+1)*fault_position/recWidth_x);
    nx =  ne_x *(num_nodes-1)+1;


    i_m =  std::floor((offset_x-fault_position)/width_x) *(num_nodes-1);
    j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    k_m =  std::round((offset_z)/width_z) *(num_nodes-1);    
  }
  
  
    // if(n == 0){
    //   blockWidth_x=(fault_position-a_x);
    //   ne_x = std::round((ne_x+1)*fault_ref);
    //   nx =  ne_x *(num_nodes-1)+1;
      
    //   i_m =  std::round((offset_x)/width_x) *(num_nodes-1);
    //   j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    //   k_m =  std::round((offset_z)/width_z) *(num_nodes-1);
    // }else{
    //   double ne_x0 = std::round((ne_x+1)*fault_ref);
    //   ne_x = std::round(1/dx[0])-ne_x0;
    //   nx =  ne_x *(num_nodes-1)+1;
    //   blockWidth_x= (b_x-fault_position);

    //   i_m =  std::floor((offset_x-fault_ref)/width_x) *(num_nodes-1);
    //   j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    //   k_m =  std::round((offset_z)/width_z) *(num_nodes-1);    
    // }

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
    if( level <= getCoarsestMeshLevel()){

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

    }
    double* curvilinear_x = new double[num_nodes*num_nodes*num_nodes];
    double* curvilinear_y = new double[num_nodes*num_nodes*num_nodes];
    double* curvilinear_z = new double[num_nodes*num_nodes*num_nodes];  

 	
    int i_p = i_m + num_nodes;
    int j_p = j_m + num_nodes;
    int k_p = k_m + num_nodes;   
    if( level <= getCoarsestMeshLevel()){

  
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
    delete left_bnd_x;
    delete left_bnd_y;
    delete left_bnd_z;
    delete right_bnd_x;
    delete right_bnd_y;
    delete right_bnd_z;
    delete bottom_bnd_x;
    delete bottom_bnd_y;
    delete bottom_bnd_z;
    delete top_bnd_x;
    delete top_bnd_y;
    delete top_bnd_z;
    delete front_bnd_x;
    delete front_bnd_y;
    delete front_bnd_z;
    delete back_bnd_x;
    delete back_bnd_y;
    delete back_bnd_z;
    }

  
  
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

    if( level <= getCoarsestMeshLevel()){
  
    metricDerivativesAndJacobian3D(num_nodes,
				   curvilinear_x,  curvilinear_y,  curvilinear_z,
				   gl_vals_x,  gl_vals_y,  gl_vals_z,
				   q_x,  q_y,  q_z,
				   r_x,  r_y,  r_z,
				   s_x,  s_y,  s_z,				  
				   jacobian,
				   width_x,  width_y,  width_z
				   );
    }

    for (int k=0; k< num_nodes; k++){
      for (int j=0; j< num_nodes; j++){
	for (int i=0; i< num_nodes; i++){
	
	  double x= gl_vals_x[id_xyz(k,j,i)];
	  double y= gl_vals_y[id_xyz(k,j,i)];
	  double z= gl_vals_z[id_xyz(k,j,i)];


	  // velocity //Gaussian
	  luh[id_xyzf(k,j,i,0)]  = 0*std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	  luh[id_xyzf(k,j,i,1)]  = 0;
	  luh[id_xyzf(k,j,i,2)]  = 0;
	
	  // stress field
	  luh[id_xyzf(k,j,i,3)]  = 0;
	  luh[id_xyzf(k,j,i,4)]  = 0;
	  luh[id_xyzf(k,j,i,5)]  = 0;
	  luh[id_xyzf(k,j,i,6)]  = 0;
	  luh[id_xyzf(k,j,i,7)]  = 0;
	  luh[id_xyzf(k,j,i,8)]  = 0;

	  //Slip
	  luh[id_xyzf(k,j,i,9)]  = 0;

	  // material parameters
	  luh[id_xyzf(k,j,i,10)] = 2.7; //rho
	  luh[id_xyzf(k,j,i,11)] = 6.0; //cp
	  luh[id_xyzf(k,j,i,12)] = 3.343; //cs

	  if( level <= getCoarsestMeshLevel()){

	  // Jacobian
	  luh[id_xyzf(k,j,i,13)]  = jacobian[id_xyz(k,j,i)];	  

	  // q_x, q_y, q_z
	  luh[id_xyzf(k,j,i,14)]  = q_x[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,15)]  = q_y[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,16)] =  q_z[id_xyz(k,j,i)];
	  
	  // r_x, r_y, r_z
	  luh[id_xyzf(k,j,i,17)] = r_x[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,18)] = r_y[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,19)] = r_z[id_xyz(k,j,i)];
	  
	  // s_x, s_y, s_z
	  luh[id_xyzf(k,j,i,20)] = s_x[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,21)] = s_y[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,22)] = s_z[id_xyz(k,j,i)];
	  
	  // x,y,z
	  luh[id_xyzf(k,j,i,23)] = gl_vals_x[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,24)] = gl_vals_y[id_xyz(k,j,i)];
	  luh[id_xyzf(k,j,i,25)] = gl_vals_z[id_xyz(k,j,i)];	  
	  }
	}
      }
    }
    
    delete gl_vals_x;
    delete gl_vals_y;
    delete gl_vals_z;
    delete jacobian;
    delete q_x;
    delete q_y;
    delete q_z;
    delete r_x;
    delete r_y;
    delete r_z;
    delete s_x;
    delete s_y;
    delete s_z;
    
  }
}

void Elastic::MyElasticWaveSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 16
  
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
  constexpr int basisSize=Order+1;
  int numberOfData=MyElasticWaveSolver::NumberOfParameters+MyElasticWaveSolver::NumberOfVariables;

  if (!tarch::la::equals(t,0.0)) {
    return exahype::solvers::Solver::RefinementControl::Keep;
  }

  kernels::idx4 idx_Q(basisSize,basisSize,basisSize,numberOfData);



  //check if element is near to fault

  bool is_fault=false;

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
      // TODO: Change position to static parameter
	double jacobian=luh[idx_Q(i,j,k,23)];
	is_fault = is_fault || std::abs(luh[idx_Q(i,j,k,23)] - 5) < dx[0]*1.0/(Order) + 1.0e-5;
      }
    }
  }

  if(is_fault){
    return exahype::solvers::Solver::RefinementControl::Refine;
  }
  
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void Elastic::MyElasticWaveSolver::eigenvalues(const double* const Q,const int d,double* const lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 16
  double cp = Q[11];
  double cs = Q[12];

  double q_x;
  double q_y;
  double q_z;  
  double r_x;
  double r_y;
  double r_z;
  double s_x;
  double s_y;
  double s_z;  
  
  extractTransformation(Q,q_x,q_y,q_z,r_x,r_y,r_z,s_x,s_y,s_z);
   
  lambda[0] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*cp;
  lambda[1] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*cs;
  lambda[2] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*cs;
  
  lambda[3] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*cp;
  lambda[4] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*cs;
  lambda[5] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*cs;
  
  lambda[6] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*cp;
  lambda[7] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*cs;
  lambda[8] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*cs;
}


void Elastic::MyElasticWaveSolver::flux(const double* const Q,double** const F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 9 + 16
  
  double sigma_xx=Q[3];
  double sigma_yy=Q[4];
  double sigma_zz=Q[5];  
  double sigma_xy=Q[6];
  double sigma_xz=Q[7];
  double sigma_yz=Q[8];

  double jacobian=Q[13];

  double q_x;
  double q_y;
  double q_z;  
  double r_x;
  double r_y;
  double r_z;
  double s_x;
  double s_y;
  double s_z;  
  
  extractTransformation(Q,q_x,q_y,q_z,r_x,r_y,r_z,s_x,s_y,s_z);

  
  F[0][ 0] = -jacobian*(q_x*sigma_xx+q_y*sigma_xy+q_z*sigma_xz);
  F[0][ 1] = -jacobian*(q_x*sigma_xy+q_y*sigma_yy+q_z*sigma_yz);
  F[0][ 2] = -jacobian*(q_x*sigma_xz+q_y*sigma_yz+q_z*sigma_zz);
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;
  F[0][ 6] = 0.0;
  F[0][ 7] = 0.0;
  F[0][ 8] = 0.0;
  F[0][ 9] = 0.0;  
  
  F[1][ 0] = -jacobian*(r_x*sigma_xx+r_y*sigma_xy+r_z*sigma_xz);
  F[1][ 1] = -jacobian*(r_x*sigma_xy+r_y*sigma_yy+r_z*sigma_yz);
  F[1][ 2] = -jacobian*(r_x*sigma_xz+r_y*sigma_yz+r_z*sigma_zz);
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
  F[1][ 5] = 0.0;
  F[1][ 6] = 0.0;
  F[1][ 7] = 0.0;
  F[1][ 8] = 0.0;
  F[1][ 9] = 0.0;  

  F[2][ 0] = -jacobian*(s_x*sigma_xx+s_y*sigma_xy+s_z*sigma_xz);
  F[2][ 1] = -jacobian*(s_x*sigma_xy+s_y*sigma_yy+s_z*sigma_yz);
  F[2][ 2] = -jacobian*(s_x*sigma_xz+s_y*sigma_yz+s_z*sigma_zz);
  F[2][ 3] = 0.0;
  F[2][ 4] = 0.0;
  F[2][ 5] = 0.0;
  F[2][ 6] = 0.0;
  F[2][ 7] = 0.0;
  F[2][ 8] = 0.0;
  F[2][ 9] = 0.0;  
}



void  Elastic::MyElasticWaveSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // @todo Please implement/augment if required
  
  constexpr int numberOfVariables  = MyElasticWaveSolver::NumberOfVariables;
  kernels::idx2 idx_BgradQ(3,numberOfVariables);
  double u_q = gradQ[idx_BgradQ(0,0)];
  double v_q = gradQ[idx_BgradQ(0,1)];
  double w_q = gradQ[idx_BgradQ(0,2)];

  double u_r = gradQ[idx_BgradQ(1,0)];
  double v_r = gradQ[idx_BgradQ(1,1)];
  double w_r = gradQ[idx_BgradQ(1,2)];

  double u_s = gradQ[idx_BgradQ(2,0)];
  double v_s = gradQ[idx_BgradQ(2,1)];
  double w_s = gradQ[idx_BgradQ(2,2)];

  double q_x;
  double q_y;
  double q_z;  
  double r_x;
  double r_y;
  double r_z;
  double s_x;
  double s_y;
  double s_z;  
  
  extractTransformation(Q,q_x,q_y,q_z,r_x,r_y,r_z,s_x,s_y,s_z);


  BgradQ[idx_BgradQ(0,0)] = 0;
  BgradQ[idx_BgradQ(0,1)] = 0;
  BgradQ[idx_BgradQ(0,2)] = 0;  
  BgradQ[idx_BgradQ(0,3)] = -q_x*u_q;
  BgradQ[idx_BgradQ(0,4)] = -q_y*v_q;
  BgradQ[idx_BgradQ(0,5)] = -q_z*w_q;
  BgradQ[idx_BgradQ(0,6)] = -(q_y*u_q+q_x*v_q); //sigma_xy
  BgradQ[idx_BgradQ(0,7)] = -(q_z*u_q+q_x*w_q); //sigma_xz
  BgradQ[idx_BgradQ(0,8)] = -(q_z*v_q+q_y*w_q); //sigma_yz
  BgradQ[idx_BgradQ(0,9)] = 0;

  BgradQ[idx_BgradQ(1,0)] = 0;
  BgradQ[idx_BgradQ(1,1)] = 0;
  BgradQ[idx_BgradQ(1,2)] = 0;  
  BgradQ[idx_BgradQ(1,3)] = -r_x*u_r;
  BgradQ[idx_BgradQ(1,4)] = -r_y*v_r;
  BgradQ[idx_BgradQ(1,5)] = -r_z*w_r;
  BgradQ[idx_BgradQ(1,6)] = -(r_y*u_r+r_x*v_r); //sigma_xy
  BgradQ[idx_BgradQ(1,7)] = -(r_z*u_r+r_x*w_r); //sigma_xz
  BgradQ[idx_BgradQ(1,8)] = -(r_z*v_r+r_y*w_r); //sigma_yz
  BgradQ[idx_BgradQ(1,9)] = 0;
      
  BgradQ[idx_BgradQ(2,0)] = 0;
  BgradQ[idx_BgradQ(2,1)] = 0;
  BgradQ[idx_BgradQ(2,2)] = 0;  
  BgradQ[idx_BgradQ(2,3)] = -s_x*u_s;
  BgradQ[idx_BgradQ(2,4)] = -s_y*v_s;
  BgradQ[idx_BgradQ(2,5)] = -s_z*w_s;
  BgradQ[idx_BgradQ(2,6)] = -(s_y*u_s+s_x*v_s); //sigma_xy
  BgradQ[idx_BgradQ(2,7)] = -(s_z*u_s+s_x*w_s); //sigma_xz
  BgradQ[idx_BgradQ(2,8)] = -(s_z*v_s+s_y*w_s); //sigma_yz    
  BgradQ[idx_BgradQ(2,9)] =0;
}

void Elastic::MyElasticWaveSolver::initPointSources() {

  double a_x = 0.0;
  double a_y = 0.0;
  double a_z = 0.0;
  
  double b_x = 10.0;
  double b_y = 10.0;
  double b_z = 10.0;
 
  double blockWidth_y = (b_y-a_y);
  double blockWidth_x = (b_x-a_x);
  double blockWidth_z = (b_z-a_z);

  double x1 = 5.0;
  double y1 = 5.0;
  double z1 = 5.0;
  
  double fault_ref_x = (x1-a_x)/(blockWidth_x);
  double fault_ref_y = (y1-a_y)/(blockWidth_y);
  double fault_ref_z = (z1-a_z)/(blockWidth_z);
  
  pointSourceLocation[0][0] = fault_ref_x;
  pointSourceLocation[0][1] = fault_ref_y;
  pointSourceLocation[0][2] =  fault_ref_z;

  x1 = 7.5;
  y1 = 5.0;
  z1 = 5.0;
  
  fault_ref_x = (x1-a_x)/(blockWidth_x);
  fault_ref_y = (y1-a_y)/(blockWidth_y);
  fault_ref_z = (z1-a_z)/(blockWidth_z);
  
  pointSourceLocation[1][0] = fault_ref_x;
  pointSourceLocation[1][1] = fault_ref_y;
  pointSourceLocation[1][2] = fault_ref_z;

  x1 = 5.0;
  y1 = 2.5;
  z1 = 5.0;
  
  fault_ref_x = (x1-a_x)/(blockWidth_x);
  fault_ref_y = (y1-a_y)/(blockWidth_y);
  fault_ref_z = (z1-a_z)/(blockWidth_z);
  
  pointSourceLocation[2][0] = fault_ref_x;
  pointSourceLocation[2][1] = fault_ref_y;
  pointSourceLocation[2][2] = fault_ref_z;
  
  x1 = 5.0;
  y1 = 7.5;
  z1 = 5.0;
  
  fault_ref_x = (x1-a_x)/(blockWidth_x);
  fault_ref_y = (y1-a_y)/(blockWidth_y);
  fault_ref_z = (z1-a_z)/(blockWidth_z);
  
  pointSourceLocation[3][0] = fault_ref_x;
  pointSourceLocation[3][1] = fault_ref_y;
  pointSourceLocation[3][2] = fault_ref_z;
}


void Elastic::MyElasticWaveSolver::pointSource(const double* const Q,const double* const x,const double t,const double dt, double* const forceVector,int n) {

  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 0.0; // Zeroed out !

 
  if(n == 0){
    
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
   
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 1.*f;
    forceVector[4] = 1.*f;
    forceVector[5] = 1.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
    forceVector[9] = 0.0;
    
  }else if(n == 1){
    
    f = 0*M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.*f;
    forceVector[5] = 0.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
    forceVector[9] = 0.0;
    
  }else if(n == 2){
    
    f = 0*M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

   
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.*f;
    forceVector[5] = 0.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
    forceVector[9] = 0.0;
    
  }else if(n == 3){
    
    f = 0*M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));

    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.*f;
    forceVector[5] = 0.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
    forceVector[9] = 0.0;
  }
}

    /**
     * @TODO LR : document
     */
void Elastic::MyElasticWaveSolver::multiplyMaterialParameterMatrix(const double* const Q, double* const rhs) {
  // double rho = Q[10];  
  // double c_p = Q[11];
  // double c_s = Q[12];  
  // double mu     = rho*c_s*c_s;
  // double lambda = rho*c_p*c_p-2*mu;
  // double jacobian = Q[13];  
  // double rho_jacobian_inv=1.0/(rho*jacobian);


  constexpr int numberOfVariables  = MyElasticWaveSolver::NumberOfVariables;  
  kernels::idx2 idx_rhs(3,numberOfVariables);
  
  double rho   = Q[10];  
  double c_p   = Q[11];
  double c_s   = Q[12];

   double mu     = rho*c_s*c_s;
   double lambda = rho*c_p*c_p-2*mu;

  double jacobian = Q[13];  
 

  double rho_jacobian_inv=1.0/(rho*jacobian);

  rhs[idx_rhs(0,0)]=rho_jacobian_inv * rhs[idx_rhs(0,0)];
  rhs[idx_rhs(0,1)]=rho_jacobian_inv * rhs[idx_rhs(0,1)];
  rhs[idx_rhs(0,2)]=rho_jacobian_inv * rhs[idx_rhs(0,2)];
  
  double lam_temp = lambda * (rhs[idx_rhs(0,3)] + rhs[idx_rhs(0,4)] + rhs[idx_rhs(0,5)]);

  rhs[idx_rhs(0,3)]=(2*mu) * rhs[idx_rhs(0,3)] +lam_temp;
  rhs[idx_rhs(0,4)]=(2*mu) * rhs[idx_rhs(0,4)] +lam_temp;
  rhs[idx_rhs(0,5)]=(2*mu) * rhs[idx_rhs(0,5)] +lam_temp;

  rhs[idx_rhs(0,6)]= mu*rhs[idx_rhs(0,6)];
  rhs[idx_rhs(0,7)]= mu*rhs[idx_rhs(0,7)];
  rhs[idx_rhs(0,8)]= mu*rhs[idx_rhs(0,8)];

  rhs[idx_rhs(1,0)] =rho_jacobian_inv * rhs[idx_rhs(1,0)];
  rhs[idx_rhs(1,1)]=rho_jacobian_inv * rhs[idx_rhs(1,1)];
  rhs[idx_rhs(1,2)]=rho_jacobian_inv * rhs[idx_rhs(1,2)];

  
  lam_temp = lambda * (rhs[idx_rhs(1,3)] + rhs[idx_rhs(1,4)] + rhs[idx_rhs(1,5)]);

  rhs[idx_rhs(1,3)]=(2*mu) * rhs[idx_rhs(1,3)] +lam_temp;
  rhs[idx_rhs(1,4)]=(2*mu) * rhs[idx_rhs(1,4)] +lam_temp;
  rhs[idx_rhs(1,5)]=(2*mu) * rhs[idx_rhs(1,5)] +lam_temp;

  rhs[idx_rhs(1,6)]= mu*rhs[idx_rhs(1,6)];
  rhs[idx_rhs(1,7)]= mu*rhs[idx_rhs(1,7)];
  rhs[idx_rhs(1,8)]= mu*rhs[idx_rhs(1,8)];

  rhs[idx_rhs(2,0)] =rho_jacobian_inv * rhs[idx_rhs(2,0)];
  rhs[idx_rhs(2,1)]=rho_jacobian_inv * rhs[idx_rhs(2,1)];
  rhs[idx_rhs(2,2)]=rho_jacobian_inv * rhs[idx_rhs(2,2)];

  
  lam_temp = lambda * (rhs[idx_rhs(2,3)] + rhs[idx_rhs(2,4)] + rhs[idx_rhs(2,5)]);

  rhs[idx_rhs(2,3)]=(2*mu) * rhs[idx_rhs(2,3)] +lam_temp;
  rhs[idx_rhs(2,4)]=(2*mu) * rhs[idx_rhs(2,4)] +lam_temp;
  rhs[idx_rhs(2,5)]=(2*mu) * rhs[idx_rhs(2,5)] +lam_temp;

  rhs[idx_rhs(2,6)]= mu*rhs[idx_rhs(2,6)];
  rhs[idx_rhs(2,7)]= mu*rhs[idx_rhs(2,7)];
  rhs[idx_rhs(2,8)]= mu*rhs[idx_rhs(2,8)];


  // rhs[0]=rho_jacobian_inv * rhs[0];
  // rhs[1]=rho_jacobian_inv * rhs[1];
  // rhs[2]=rho_jacobian_inv * rhs[2];
  
  // double lam_temp = lambda * (rhs[3] + rhs[4] + rhs[5]);
  // rhs[3]=(2*mu) * rhs[3] +lam_temp;
  // rhs[4]=(2*mu) * rhs[4] +lam_temp;
  // rhs[5]=(2*mu) * rhs[5] +lam_temp;

  // rhs[6]= mu*rhs[6];
  // rhs[7]= mu*rhs[7];
  // rhs[8]= mu*rhs[8];


  // rhs[9] =rho_jacobian_inv * rhs[9];
  // rhs[10]=rho_jacobian_inv * rhs[10];
  // rhs[11]=rho_jacobian_inv * rhs[11];

  
  // lam_temp = lambda * (rhs[12] + rhs[13] + rhs[14]);

  // rhs[12]=(2*mu) * rhs[12] +lam_temp;
  // rhs[13]=(2*mu) * rhs[13] +lam_temp;
  // rhs[14]=(2*mu) * rhs[14] +lam_temp;

  // rhs[15]= mu*rhs[15];
  // rhs[16]= mu*rhs[16];
  // rhs[17]= mu*rhs[17];

  // rhs[18] = rho_jacobian_inv * rhs[18];
  // rhs[19] = rho_jacobian_inv * rhs[19];
  // rhs[20] = rho_jacobian_inv * rhs[20];
  
  // lam_temp = lambda * (rhs[21] + rhs[22] + rhs[23]);

  // rhs[21]=(2*mu) * rhs[21] +lam_temp;
  // rhs[22]=(2*mu) * rhs[22] +lam_temp;
  // rhs[23]=(2*mu) * rhs[23] +lam_temp;

  // rhs[24]= mu*rhs[24];
  // rhs[25]= mu*rhs[25];
  // rhs[26]= mu*rhs[26];  
}

void Elastic::MyElasticWaveSolver::extractTransformation(const double* const Q,
							 double& q_x,double& q_y,double& q_z,
							 double& r_x,double& r_y,double& r_z,
							 double& s_x,double& s_y,double& s_z) {
  q_x=Q[14];
  q_y=Q[15];
  q_z=Q[16];
  r_x=Q[17];
  r_y=Q[18];
  r_z=Q[19];
  s_x=Q[20];
  s_y=Q[21];
  s_z=Q[22];
}

void Elastic::MyElasticWaveSolver::riemannSolver(double* const FL,double* const FR,const double* const QL,const double* const QR,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex){
  constexpr int numberOfVariables  = MyElasticWaveSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyElasticWaveSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyElasticWaveSolver::Order+1;
  constexpr int order              = basisSize - 1;

  kernels::idx3 idx_QLR(basisSize,basisSize,numberOfData);
  kernels::idx3 idx_FLR(basisSize,basisSize,numberOfVariables);

  //check if the face lies on a fault
  bool is_fault=true;
   for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      // TODO: Change position to static parameter
      is_fault=is_fault && std::abs(QL[idx_QLR(i,j,23)] - 5) < 1.0e-5; 
    }
  }

  double n_p[3]={0,0,0};
  double n_m[3]={0,0,0};

  double m_p[3]={0,0,0};
  double m_m[3]={0,0,0};

  double l_p[3]={0,0,0};  
  double l_m[3]={0,0,0};

  double norm_p_qr;
  double norm_m_qr;
  
  double FLn, FLm, FLl, FRn,FRm,FRl;
  double FL_n,FL_m,FL_l,FR_n,FR_m,FR_l;
  double FLx,FLy,FLz,FRx,FRy,FRz;
  double FL_x,FL_y,FL_z,FR_x,FR_y,FR_z;

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      double rho_p=QR[idx_QLR(i,j,10)];
      double c_p_p=QR[idx_QLR(i,j,11)];
      double c_s_p=QR[idx_QLR(i,j,12)];
      double mu_p=c_s_p*c_s_p*rho_p;
      double lam_p = rho_p*c_p_p*c_p_p-2*mu_p;      

      double rho_m=QL[idx_QLR(i,j,10)];
      double c_p_m=QL[idx_QLR(i,j,11)];
      double c_s_m=QL[idx_QLR(i,j,12)];
      double mu_m=c_s_m*c_s_m*rho_m;
      double lam_m = rho_m*c_p_m*c_p_m-2*mu_m;      

      get_normals(normalNonZeroIndex, norm_p_qr, n_p, QR + idx_QLR(i,j,0));
      get_normals(normalNonZeroIndex, norm_m_qr, n_m, QL + idx_QLR(i,j,0));    

      double Tx_m,Ty_m,Tz_m,Tx_p,Ty_p,Tz_p;
      double vx_m,vy_m,vz_m,vx_p,vy_p,vz_p;
      
      extract_tractions_and_particle_velocity(n_p,QR+idx_QLR(i,j,0),Tx_p,Ty_p,Tz_p,vx_p,vy_p,vz_p );
      extract_tractions_and_particle_velocity(n_m,QL+idx_QLR(i,j,0),Tx_m,Ty_m,Tz_m,vx_m,vy_m,vz_m ); 
      
      localBasis(n_p, m_p, l_p, 3);
      localBasis(n_m, m_m, l_m, 3);

      double Tn_m,Tm_m,Tl_m,vn_m,vm_m,vl_m;
      double Tn_p,Tm_p,Tl_p,vn_p,vm_p,vl_p;

      // rotate fields into l, m, n basis
      rotate_into_orthogonal_basis(n_m,m_m,l_m,Tx_m,Ty_m,Tz_m,Tn_m,Tm_m,Tl_m);
      rotate_into_orthogonal_basis(n_m,m_m,l_m,vx_m,vy_m,vz_m,vn_m,vm_m,vl_m);
      rotate_into_orthogonal_basis(n_p,m_p,l_p,Tx_p,Ty_p,Tz_p,Tn_p,Tm_p,Tl_p);
      rotate_into_orthogonal_basis(n_p,m_p,l_p,vx_p,vy_p,vz_p,vn_p,vm_p,vl_p);      
      
  
      // extract local s-wave and p-wave impedances
      double zs_p=rho_p*c_s_p;
      double zp_p=rho_p*c_p_p;      
      double zs_m=rho_m*c_s_m;
      double zp_m=rho_m*c_p_m;
      
      // impedance must be greater than zero !
      if (zp_p <= 0.0 || zp_m <= 0.0){
	std::cout<<zs_p<<" "<<zs_m<<" "<<zp_p<<" "<<zp_m<<"\n";
	std::cout<<" Impedance must be greater than zero ! "<< std::endl;
	std::exit(-1);
      }

      // generate interface data preserving the amplitude of the outgoing charactertritics
      // and satisfying interface conditions exactly.
      double vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p;        
      double vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m;


      if(is_fault){
	double x[3] = {QR[idx_QLR(i,j,23)], QR[idx_QLR(i,j,24)], QR[idx_QLR(i,j,25)]};
	double S = QR[idx_QLR(i,j,9)];
	SlipWeakeningFriction(vn_p,vn_m, Tn_p,Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p,Tn_hat_m, vm_p,vm_m,
			      Tm_p,Tm_m, zs_p,zs_m, vm_hat_p, vm_hat_m, Tm_hat_p,Tm_hat_m, vl_p,vl_m,Tl_p,Tl_m, zs_p,
			      zs_m, vl_hat_p , vl_hat_m, Tl_hat_p,Tl_hat_m, l_p, m_p, n_p, x, S);


      }else{
	if (isBoundaryFace) {
	  double r= faceIndex==1 ? 1 : 0;
	  riemannSolver_boundary(faceIndex,r,vn_m,vm_m,vl_m,Tn_m,Tm_m,Tl_m,zp_m,zs_m,vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m);
	  riemannSolver_boundary(faceIndex,r,vn_p,vm_p,vl_p,Tn_p,Tm_p,Tl_p,zp_p,zs_p,vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p);      
	}else{
	  riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
	  riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
	  riemannSolver_Nodal(vl_p,vl_m, Tl_p, Tl_m, zs_p , zs_m, vl_hat_p , vl_hat_m, Tl_hat_p, Tl_hat_m);
	}
      }
      
      //generate fluctuations in the local basis coordinates: n, m, l
      generate_fluctuations_left(zp_m,Tn_m,Tn_hat_m,vn_m,vn_hat_m,FLn);
      generate_fluctuations_left(zs_m,Tm_m,Tm_hat_m,vm_m,vm_hat_m,FLm);
      generate_fluctuations_left(zs_m,Tl_m,Tl_hat_m,vl_m,vl_hat_m,FLl);

      generate_fluctuations_right(zp_p,Tn_p,Tn_hat_p,vn_p,vn_hat_p,FRn);
      generate_fluctuations_right(zs_p,Tm_p,Tm_hat_p,vm_p,vm_hat_p,FRm);
      generate_fluctuations_right(zs_p,Tl_p,Tl_hat_p,vl_p,vl_hat_p,FRl);

      FL_n = FLn/zp_m;
      if(zs_m > 0){
	FL_m = FLm/zs_m;
	FL_l = FLl/zs_m;
      }else{
	FL_m=0;
	FL_l=0;
      }
    
      FR_n = FRn/zp_p;
      if(zs_p > 0){    
	FR_m = FRm/zs_p;
	FR_l = FRl/zs_p;
      }else{
	FR_m=0;
	FR_l=0;
      }
    
      // rotate back to the physical coordinates x, y, z
      rotate_into_physical_basis(n_m,m_m,l_m,FLn,FLm,FLl,FLx,FLy,FLz);
      rotate_into_physical_basis(n_p,m_p,l_p,FRn,FRm,FRl,FRx,FRy,FRz);
      rotate_into_physical_basis(n_m,m_m,l_m,FL_n,FL_m,FL_l,FL_x,FL_y,FL_z);
      rotate_into_physical_basis(n_p,m_p,l_p,FR_n,FR_m,FR_l,FR_x,FR_y,FR_z);
     
      // construct flux fluctuationn vectors obeying the eigen structure of the PDE
      // and choose physically motivated penalties such that we can prove
      // numerical stability.

      FR[idx_FLR(i,j, 0)] = norm_p_qr/rho_p*FRx;
      FL[idx_FLR(i,j, 0)] = norm_m_qr/rho_m*FLx;
    
      FR[idx_FLR(i,j, 1)] = norm_p_qr/rho_p*FRy;
      FL[idx_FLR(i,j, 1)] = norm_m_qr/rho_m*FLy;

      FR[idx_FLR(i,j, 2)] = norm_p_qr/rho_p*FRz;
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

      if(is_fault){
	double slip_rate= std::sqrt((vl_hat_p- vl_hat_m)*(vl_hat_p- vl_hat_m)+(vm_hat_p- vm_hat_m)*(vm_hat_p- vm_hat_m));
	FL[idx_FLR(i,j, 9)] = -slip_rate;
	FR[idx_FLR(i,j, 9)] = -slip_rate;
      }else{
	FL[idx_FLR(i,j, 9)] = 0;
	FR[idx_FLR(i,j, 9)] = 0;
      }
    }    
  }
}


//Gram Schmidt orthonormalization
void Elastic::MyElasticWaveSolver::Gram_Schmidt(double* const y, double* const z){
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
  }
  
  double norm_z = std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  
  for (int i = 0; i< 3; i++){
    z[i] =  z[i]/norm_z;
  }
}

void Elastic::MyElasticWaveSolver::localBasis(double* const n, double* const m, double* const l, int d){
  if (d == 2){
      l[0] = 0.;
      l[1] = 0.;
      l[2] = 1.0;
      
      m[0] = n[1]*l[2]-n[2]*l[1];
      m[1] = -(n[0]*l[2]-n[2]*l[0]);
      m[2] = n[0]*l[1]-n[1]*l[0];
  }else if (d == 3){
      double tol, diff_norm1, diff_norm2;
      tol = 1e-12;
      m[0] = 0.;
      m[1] = 1.;
      m[2] = 0.;
      
      diff_norm1 =  std::sqrt(pow(n[0]-m[0],2) + pow(n[1]-m[1], 2) + pow(n[2]-m[2], 2));
      diff_norm2 =  std::sqrt(pow(n[0]+m[0],2) + pow(n[1]+m[1], 2) + pow(n[2]+m[2], 2));
      
      if (diff_norm1 >= tol && diff_norm2 >= tol){
      	Gram_Schmidt(n, m);
      }else{
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

void Elastic::MyElasticWaveSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
   double p = 0.5*(z*v + sigma);
   if(z > 0){
     v_hat = (1+r)/z*p;
     sigma_hat = (1-r)*p;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }
}

void Elastic::MyElasticWaveSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
   double q = 0.5*(z*v - sigma);
   if(z > 0){
     v_hat = (1+r)/z*q;
     sigma_hat = -(1-r)*q;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }
}

void Elastic::MyElasticWaveSolver::get_normals(int normalNonZeroIndex,double& norm, double* const n,const double* const Q){

  double q_x;
  double q_y;
  double q_z;  
  double r_x;
  double r_y;
  double r_z;
  double s_x;
  double s_y;
  double s_z;  
  
  extractTransformation(Q,q_x,q_y,q_z,r_x,r_y,r_z,s_x,s_y,s_z);
  
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

void Elastic::MyElasticWaveSolver::extract_tractions_and_particle_velocity(double* const n,const double* const Q, double& Tx,double& Ty,double& Tz,double& vx,double& vy,double& vz ){
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

void Elastic::MyElasticWaveSolver::rotate_into_orthogonal_basis(double* const n,double* const m,double* const l, double Tx,double Ty,double Tz, double& Tn, double& Tm, double& Tl){
    Tn= Tx*n[0] + Ty*n[1] + Tz*n[2];
    Tm= Tx*m[0] + Ty*m[1] + Tz*m[2];
    Tl= Tx*l[0] + Ty*l[1] + Tz*l[2];
}

void Elastic::MyElasticWaveSolver::rotate_into_physical_basis(double* const n,double* const m,double* const l, double Fn,double Fm,double Fl, double& Fx, double& Fy, double& Fz){
  Fx = n[0]*Fn + m[0]*Fm + l[0]*Fl;
  Fy = n[1]*Fn + m[1]*Fm + l[1]*Fl;
  Fz = n[2]*Fn + m[2]*Fm + l[2]*Fl;
}

void Elastic::MyElasticWaveSolver::generate_fluctuations_left(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) + (T-T_hat));
}

void Elastic::MyElasticWaveSolver::generate_fluctuations_right(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) - (T-T_hat));
}

void Elastic::MyElasticWaveSolver::riemannSolver_boundary(int faceIndex,double r, double vn , double vm , double vl, double Tn , double Tm ,double Tl , double zp, double zs , double& vn_hat , double& vm_hat ,double& vl_hat , double& Tn_hat , double& Tm_hat ,double& Tl_hat)
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

// Rupture Dynamics
void Elastic::MyElasticWaveSolver::SlipWeakeningFriction(double vn_p,double vn_m, double Tn_p, double Tn_m, double zn_p , double zn_m, double& vn_hat_p , double& vn_hat_m, double& Tn_hat_p, double& Tn_hat_m, double vm_p,double vm_m, double Tm_p, double Tm_m, double zl_p , double zl_m, double& vm_hat_p , double& vm_hat_m, double& Tm_hat_p, double& Tm_hat_m, double vl_p,double vl_m, double Tl_p, double Tl_m, double zm_p , double zm_m, double& vl_hat_p , double& vl_hat_m, double& Tl_hat_p, double& Tl_hat_m, double* const l, double* const m, double* const n, double* const x,  double S)
{
  // compute characteristics
  double p_l = zl_p*vl_p + Tl_p;
  double p_m = zm_p*vm_p + Tm_p;
  double p_n = zn_p*vn_p + Tn_p;
  
  double q_l = zl_m*vl_m - Tl_m;
  double q_m = zm_m*vm_m - Tm_m;
  double q_n = zn_m*vn_m - Tn_m;

  // half of the harmonic mean of Z1_s, Z2_s
  double eta_s=(zl_p*zl_m)/(zl_p+zl_m);                                    
  double eta_n=(zn_p*zn_m)/(zn_p+zn_m); 
  
  double  phi_l= eta_s*(p_l/zl_p - q_l/zl_m);
  double  phi_m= eta_s*(p_m/zm_p - q_m/zm_m);
  double  phi_n= eta_n*(p_n/zn_p - q_n/zn_m);
    
  double T0_l=0;
  double T0_m=0;
  double T0_n=0;

  // get prestress (where normal traction is effective normal traction)
  prestress(T0_n, T0_m, T0_l, x, 0.0, l, m, n);

  vn_hat_m = (q_n - phi_n)/zn_p;   //! continuity of normal velocity
  Tn_hat_m = phi_n;                //! continuity of normal stress
  vn_hat_p = (p_n + phi_n)/zn_m;   //! continuity of normal velocity
  Tn_hat_p = phi_n;

  double tau_lock = std::sqrt(std::pow(T0_l + phi_l, 2) + std::pow(T0_m + phi_m, 2));
  double sigma_n = std::max(0.0, -(T0_n + phi_n));   // including prestress
  double tau_str;
  double Vel = 0.0;
  double Tl, Tm, vv_l, vv_m; 
  TauStrength(tau_str, sigma_n, S, x, 0.0);
  
  if (tau_lock >= tau_str){
    slip_weakening(vv_m,  vv_l, Vel, Tm, Tl, phi_m+T0_m, phi_l+T0_l, eta_s, tau_str, sigma_n);

    Tm_hat_m = Tm - T0_m;
    Tl_hat_m = Tl - T0_l;

    Tm_hat_p = Tm - T0_m;
    Tl_hat_p = Tl - T0_l;
  }else{
    Tm_hat_m = phi_m;
    Tl_hat_m = phi_l;

    Tm_hat_p = phi_m;
    Tl_hat_p = phi_l;
       
    vv_m = 0.0;
    vv_l = 0.0;
    
    Vel = 0.0;
  }
  vm_hat_p = (Tm_hat_m + p_m)/zm_m + vv_m;
  vl_hat_p = (Tl_hat_m + p_l)/zl_m + vv_l;
    
  vm_hat_m = (q_m - Tm_hat_p)/zm_p - vv_m;
  vl_hat_m = (q_l - Tl_hat_p)/zl_p - vv_l;
}

  // solve for slip-rate (vv):  
void Elastic::MyElasticWaveSolver::slip_weakening(double& v1, double& v2, double& Vel, double& tau1, double& tau2,
								double phi_1, double phi_2, double eta, double tau_str, double sigma_n){
  
  double Phi = std::sqrt(std::pow(phi_1, 2) + std::pow(phi_2, 2));   // stress-transfer functional
  Vel = (Phi - tau_str)/eta;                // slip-rate

  //compute slip velocities
  v1 = phi_1/(eta+tau_str/Vel);
  v2 = phi_2/(eta+tau_str/Vel);
  
  //compute shear stress on the fault
  tau1 = phi_1 - eta*v1;
  tau2 = phi_2 - eta*v2;
}

void Elastic::MyElasticWaveSolver::TauStrength(double& tau_str, double sigma_n, double S, double* const x, double t)
{
  double mu_s = 0.677;                          // stastic friction
  double mu_d = 0.525;                          // dynamic friction
  double sigma0 = 120.0;                        // normal stress
  double S_c = 0.40;                            // critical slip

  mu_s = mu_s + 0*1e10*(1.0-boxcar(x[1]-5,15.0)*boxcar(x[2]-5.0,15.0));
  double fric_coeff = mu_s - (mu_s-mu_d) * std::min(S,S_c)/S_c;     // friction coefficient
  tau_str = fric_coeff*sigma_n;                               // fault strength
}

double Elastic::MyElasticWaveSolver::boxcar(double x, double w) 
{
  // f(x) is boxcar of unit amplitude in (x-w,x+w)
  double f;
  double tol = 1e-5;
  
  if (-w+tol<x and x<w-tol){    // inside
    f = 1.0;
  }
  else if (std::abs(-w-x)<=tol or std::abs(x-w)<=tol){     // boundary
    f = 0.5;
  }
  else{    // outside
    f = 0.0;
  }
  return f;
}

void Elastic::MyElasticWaveSolver::prestress(double& T0_n, double& T0_m, double& T0_l, double* const x, double t, double* const l, double* const m, double* const n)
{
  T0_n = -120.0;
  T0_m = 0.0;
  T0_l = 70.0;
  //T0_l = T0_l+11.6;
  //T0_l = T0_l+8.0   *boxcar(x[2]-5,1.5)*boxcar(x[3]-2.5,1.5);
  if ((x[1] >= 4.25 and x[1] <= 5.75) and (x[2] >= 4.25 and x[2] <= 5.75)){
    T0_l = T0_l+11.6; //*boxcar(x[2]-5,1.5);
    //T0_l = T0_l-8.0   *boxcar(x[2]-5,1.5)*boxcar(x[3]-7.5,1.5);

  }
}





