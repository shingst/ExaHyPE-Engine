#include "ElasticWaveEquation.h"

#include "ElasticWaveEquation_Variables.h"

#include "../../../ExaHyPE/kernels/KernelUtils.h"

#if !defined(_USE_NETCDF)
#if !defined(_GLL)
#include "kernels/GaussLegendreBasis.h"
#else
#include "kernels/GaussLobattoBasis.h"
#endif //_GLL
#endif //_USE_NETCDF

#ifdef OPT_KERNELS
#include "kernels/ElasticWaveEquation/converter.h"
#endif  


// We define the constructor of the actual solver here in order to regenerate it easily.
ElasticWaveEquation3D::ElasticWaveEquation::ElasticWaveEquation(double maximumMeshSize,int maximumAdaptiveMeshDepth,int DMPObservables,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs):
#if defined(_CURVILINEAR)
  curvilinear_transformation(Order+1,1/9.0),
#elif defined(_USE_NETCDF)
  reader("topography.nc",Order+1),
#endif  
  AbstractElasticWaveEquation::AbstractElasticWaveEquation(maximumMeshSize,maximumAdaptiveMeshDepth,DMPObservables,timeStepping) {
  init(cmdlineargs);
}

tarch::logging::Log ElasticWaveEquation3D::ElasticWaveEquation::_log( "ElasticWaveEquation3D::ElasticWaveEquation" );

void ElasticWaveEquation3D::ElasticWaveEquation::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void ElasticWaveEquation3D::ElasticWaveEquation::adjustSolution(double *luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 36 + 19
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
  constexpr int num_nodes = ElasticWaveEquation::Order+1;
  constexpr int numberOfData = ElasticWaveEquation::NumberOfVariables + ElasticWaveEquation::NumberOfParameters;
  kernels::idx4 id_4(num_nodes,num_nodes,num_nodes,numberOfData);
  kernels::idx3 id_3(num_nodes,num_nodes,num_nodes);

  double gl_vals_x[num_nodes*num_nodes*num_nodes];
  double gl_vals_y[num_nodes*num_nodes*num_nodes];
  double gl_vals_z[num_nodes*num_nodes*num_nodes];

  double jacobian[num_nodes*num_nodes*num_nodes];
  double q_x[num_nodes*num_nodes*num_nodes];
  double q_y[num_nodes*num_nodes*num_nodes];
  double q_z[num_nodes*num_nodes*num_nodes];
  
  double r_x[num_nodes*num_nodes*num_nodes];
  double r_y[num_nodes*num_nodes*num_nodes];
  double r_z[num_nodes*num_nodes*num_nodes];

  double s_x[num_nodes*num_nodes*num_nodes];
  double s_y[num_nodes*num_nodes*num_nodes];
  double s_z[num_nodes*num_nodes*num_nodes];
  
  int num_nodes_domain = std::round(1.0/dx[0]) * Order +1;
  
  double offset_x = center[0] -dx[0]*0.5;
  double offset_y = center[1] -dx[1]*0.5;
  double offset_z = center[2] -dx[2]*0.5;    
  
  int i_e = std::round(offset_x * 1./dx[0]);
  int j_e = std::round(offset_y * 1./dx[1]);
  int k_e = std::round(offset_z * 1./dx[2]);
 
  int i_m = i_e * Order;
  int j_m = j_e * Order;
  int k_m = k_e * Order;

#if defined(_USE_NETCDF)  
  reader.metricDerivativesAndJacobian3D(i_m,j_m,k_m,dx[0],
  					gl_vals_x,gl_vals_y,gl_vals_z,
  					jacobian,
  					q_x,q_y,q_z,
  					r_x,r_y,r_z,
  					s_x,s_y,s_z);
#elif defined(_CURVILINEAR)
  curvilinear_transformation.genCoordinates(center,dx[0],
					    gl_vals_x,gl_vals_y,gl_vals_z,
					    jacobian,
					    q_x,q_y,q_z,
					    r_x,r_y,r_z,
					    s_x,s_y,s_z);
#endif  

  double dpml =  5.0/9.;
  int n = 0;
  double tol = 1e-4; 
  
   
  double xa = 5.0/9;
  double xb = 5.-5.0/9;
   
  double ya = 5.0/9;
  double yb = 5.-5.0/9;
   
  double za = 5.0/9;
  double zb = 5.-5.0/9;
   
  double d_x = 0.0;
  double d_y = 0.0;
  double d_z = 0.0;
  
  for (int k=0; k< num_nodes; k++){
    for (int j=0; j< num_nodes; j++){
      for (int i=0; i< num_nodes; i++){
#if defined(_USE_NETCDF) || defined(_CURVILINEAR)	
	double x= gl_vals_x[id_3(k,j,i)];
	double y= gl_vals_y[id_3(k,j,i)];
	double z= gl_vals_z[id_3(k,j,i)];
#else
#if !defined(_GLL)	
	double x= offset_x + kernels::legendre::nodes[Order][i]*dx[0];
	double y= offset_y + kernels::legendre::nodes[Order][j]*dx[1];
	double z= offset_z + kernels::legendre::nodes[Order][k]*dx[2];
#else
	double x= offset_x + kernels::lobatto::nodes[Order][num_nodes-1-i]*dx[0];
	double y= offset_y + kernels::lobatto::nodes[Order][num_nodes-1-j]*dx[1];
	double z= offset_z + kernels::lobatto::nodes[Order][num_nodes-1-k]*dx[2];
#endif //_GLL	
	jacobian[id_3(k,j,i)]=1.;
	q_x[id_3(k,j,i)]=1.;
	q_y[id_3(k,j,i)]=0.;
	q_z[id_3(k,j,i)]=0.;
	r_x[id_3(k,j,i)]=0.;
	r_y[id_3(k,j,i)]=1.;
	r_z[id_3(k,j,i)]=0.;
	s_x[id_3(k,j,i)]=0.;
	s_y[id_3(k,j,i)]=0.;
	s_z[id_3(k,j,i)]=1.;

#endif //_USE_NETCDF

	double x1 = 2.5;
	double y1 = 2.5;
	double z1 = 2.5;


	//	luh[id_4(k,j,i,0)]  = std::exp(-((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))/1.0);
	luh[id_4(k,j,i,0)]  = 0;
	luh[id_4(k,j,i,1)]  = 0;
	luh[id_4(k,j,i,2)]  = 0;
	luh[id_4(k,j,i,3)]  = 0;
	luh[id_4(k,j,i,4)]  = 0;
	luh[id_4(k,j,i,5)]  = 0;  
	luh[id_4(k,j,i,6)]  = 0;
	luh[id_4(k,j,i,7)]  = 0;
	luh[id_4(k,j,i,8)]  = 0;

	for(int l=9 ; l<36; l++){
	  luh[id_4(k,j,i,l)]  = 0; //pml auxiliary variables 
	}
	
#if defined(_LOH1)
	luh[id_4(k,j,i,36)] = 2.7;   //rho
	luh[id_4(k,j,i,37)] = 6.0;    //cp
	luh[id_4(k,j,i,38)] = 3.464;  //cs

	if(x < 1.0){
	  luh[id_4(k,j,i,36)] = 2.6;   //rho
	  luh[id_4(k,j,i,37)] = 4.0;    //cp
	  luh[id_4(k,j,i,38)] = 2.0;  //cs
	}
#else
    // impulse through homgeneous medium
    luh[id_4(k,j,i,36)] = 2.67;   //rho
    luh[id_4(k,j,i,37)] = 6.0;    //cp
    luh[id_4(k,j,i,38)] = 3.343;  //cs
#endif	
	
	  // luh[id_4(k,j,i,9)]   = 1.0;   //rho
	  // luh[id_4(k,j,i,10)]  = 0.0; //c(1)
	  // luh[id_4(k,j,i,11)]  = 1.484; //c(2)
	
	  // luh[id_4(k,j,i,9)]   = 2.7;   //rho
	  // luh[id_4(k,j,i,10)] = 3.343; //c(1)
	  // luh[id_4(k,j,i,11)]  = 6.0; //c(2)

	//damping
	luh[id_4(k,j,i,39)] = 0;
	luh[id_4(k,j,i,40)] = 0;
	luh[id_4(k,j,i,41)] = 0;

	//double dpmlx = 5.0-offset_x;
	//double dpmlx = 5.0-offset_x;

	double cp = luh[id_4(k,j,i,37)];

	double d0 = 0.0*(n+1.0)*cp/(2.0*dpml)* log(1.0/tol);

#if defined(_LOH1)
	d_x=0;
#else	
	if ((xa - x) >= 0 ){
	    d_x = d0*pow((xa-x)/dpml, n);
	  }
#endif	

	if ((x-xb) >= 0.0){
	    d_x = d0*pow((x-xb)/dpml, n);
	  }

	if ((ya-y) >= 0){
	     d_y = d0*pow((ya-y)/dpml, n);
	   }

	if ((y-yb) >= 0.0){
	    d_y = d0*pow((y-yb)/dpml, n);
	  }

	if ((za- z) >= 0){
	    d_z = d0*pow((za-z)/dpml, n);
	  }

	if ((z-zb) >= 0.0){
	    d_z = d0*pow((z-zb)/dpml, n);
	  }

	  luh[id_4(k,j,i,39)] = d_x;
	  luh[id_4(k,j,i,40)] = d_y;
	  luh[id_4(k,j,i,41)] = d_z;

	// if(x < 0.5){
	//   luh[id_4(k,j,i,39)] = d0;
	// }

	// if(x > 4.5){
	//   luh[id_4(k,j,i,39)] = d0;
	// }

	// if(y < 0.5){
	//   luh[id_4(k,j,i,40)] = d0;
	// }
	
	// if(y > 4.5){
	//   luh[id_4(k,j,i,40)] = d0;
	// }

	

	// if(z < 0.5){
	//   luh[id_4(k,j,i,41)] = d0;
	// }

	// if(z > 4.5){
	//   luh[id_4(k,j,i,41)] = d0;
	// }


	
	//transformation parameter

	

	luh[id_4(k,j,i,42)] = jacobian[id_3(k,j,i)];
	luh[id_4(k,j,i,43)] = q_x[id_3(k,j,i)];
	luh[id_4(k,j,i,44)] = q_y[id_3(k,j,i)];
	luh[id_4(k,j,i,45)] = q_z[id_3(k,j,i)];
	luh[id_4(k,j,i,46)] = r_x[id_3(k,j,i)];
	luh[id_4(k,j,i,47)] = r_y[id_3(k,j,i)];
	luh[id_4(k,j,i,48)] = r_z[id_3(k,j,i)];
	luh[id_4(k,j,i,49)] = s_x[id_3(k,j,i)];
	luh[id_4(k,j,i,50)] = s_y[id_3(k,j,i)];
	luh[id_4(k,j,i,51)] = s_z[id_3(k,j,i)];     

#if defined(_USE_NETCDF) // stretching
	x=x*100;
	y=y*100;
	z=z*10;
#endif
	
	luh[id_4(k,j,i,52)] = x;
	luh[id_4(k,j,i,53)] = y;
	luh[id_4(k,j,i,54)] = z;
      }
    }
  }
  }
}  

void ElasticWaveEquation3D::ElasticWaveEquation::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 36 + 19

  // @todo Please implement/augment if required

  for(int i = 0 ; i< 55; i++){
    stateOut[ i] = stateIn[i];
  }

  for(int i = 0 ; i< 36; i++){
    fluxOut[ i] = fluxIn[ i];
  }
}

exahype::solvers::Solver::RefinementControl ElasticWaveEquation3D::ElasticWaveEquation::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {

  //double dist2=(center[0]-2.5)*(center[0]-2.5)+(center[1]-2.5)*(center[1]-2.5)+(center[2]-2.5)*(center[2]-2.5);

  if(center[0]-dx[0]*0.5 < 1.0){
    return exahype::solvers::Solver::RefinementControl::Refine;
  }else if( (center[0]-dx[0]*0.5 < 3) && (level ==  getCoarsestMeshLevel() )){
    return exahype::solvers::Solver::RefinementControl::Refine;
  }
  
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void ElasticWaveEquation3D::ElasticWaveEquation::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 36 + 19
  
  // @todo Please implement/augment if required

  double rho,c_p,c_s;
  extractMaterialParameters(Q,rho,c_p,c_s);

  double jacobian;
  double q_x,q_y,q_z;
  double r_x,r_y,r_z;
  double s_x,s_y,s_z;

  extractMetricParameters(Q,jacobian,
			  q_x,q_y,q_z,
			  r_x,r_y,r_z,
			  s_x,s_y,s_z);
  
  lambda[ 0] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*c_p;
  lambda[ 1] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*c_p;
  lambda[ 2] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*c_p;
  lambda[ 3] = 0;
  lambda[ 4] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*c_s;
  lambda[ 5] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*c_s;
  lambda[ 6] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*c_s;
  lambda[ 7] = 0;
  lambda[ 8] = 0;
  for(int i=9 ; i< 36 ; i++){
    lambda[ i] = 0;
  }
}


void ElasticWaveEquation3D::ElasticWaveEquation::flux(const double* const Q,double** F) {
  // Dimensions                        = 3
  // Number of variables + parameters  = 36 + 19
  
  // @todo Please implement/augment if required

  double u,v,w;
  double sigma_xx,sigma_yy,sigma_zz;
  double sigma_xy,sigma_xz,sigma_yz;

  extractVariables(Q,
		   u,v,w,			    
		   sigma_xx,sigma_yy,sigma_zz,
		   sigma_xy,sigma_xz,sigma_yz);

  double jacobian;
  double q_x,q_y,q_z;
  double r_x,r_y,r_z;
  double s_x,s_y,s_z;

  extractMetricParameters(Q,jacobian,
			  q_x,q_y,q_z,
			  r_x,r_y,r_z,
			  s_x,s_y,s_z);


  
  F[0][0] = -jacobian*(q_x*sigma_xx+q_y*sigma_xy+q_z*sigma_xz);
  F[0][1] = -jacobian*(q_x*sigma_xy+q_y*sigma_yy+q_z*sigma_yz);
  F[0][2] = -jacobian*(q_x*sigma_xz+q_y*sigma_yz+q_z*sigma_zz);

  for(int i = 3 ; i < 9 ; i++){
    F[0][i] = 0.0;
  }
  
  F[0][9]  = -jacobian*(q_x*sigma_xx+q_y*sigma_xy+q_z*sigma_xz);
  F[0][10] = -jacobian*(q_x*sigma_xy+q_y*sigma_yy+q_z*sigma_yz);
  F[0][11] = -jacobian*(q_x*sigma_xz+q_y*sigma_yz+q_z*sigma_zz);

  for(int i = 12 ; i < 36 ; i++){
    F[0][i] = 0.0;
  }

  F[1][0] = -jacobian*(r_x*sigma_xx+r_y*sigma_xy+r_z*sigma_xz);
  F[1][1] = -jacobian*(r_x*sigma_xy+r_y*sigma_yy+r_z*sigma_yz);
  F[1][2] = -jacobian*(r_x*sigma_xz+r_y*sigma_yz+r_z*sigma_zz);
  
  for(int i = 3 ; i < 18 ; i++){
    F[1][i] = 0.0;
  }
  
  F[1][18] = -jacobian*(r_x*sigma_xx+r_y*sigma_xy+r_z*sigma_xz);
  F[1][19] = -jacobian*(r_x*sigma_xy+r_y*sigma_yy+r_z*sigma_yz);
  F[1][20] = -jacobian*(r_x*sigma_xz+r_y*sigma_yz+r_z*sigma_zz);

  for(int i = 21 ; i < 36 ; i++){
    F[1][i] = 0.0;
  }

  F[2][0] = -jacobian*(s_x*sigma_xx+s_y*sigma_xy+s_z*sigma_xz);
  F[2][1] = -jacobian*(s_x*sigma_xy+s_y*sigma_yy+s_z*sigma_yz);
  F[2][2] = -jacobian*(s_x*sigma_xz+s_y*sigma_yz+s_z*sigma_zz);

  for(int i = 3 ; i < 27 ; i++){
    F[2][i] = 0.0;
  }
  
  F[2][27] = -jacobian*(s_x*sigma_xx+s_y*sigma_xy+s_z*sigma_xz);
  F[2][28] = -jacobian*(s_x*sigma_xy+s_y*sigma_yy+s_z*sigma_yz);
  F[2][29] = -jacobian*(s_x*sigma_xz+s_y*sigma_yz+s_z*sigma_zz);
  
  for(int i = 30 ; i < 36 ; i++){
    F[2][i] = 0.0;
  }
}


//You can either implement this method or modify fusedSource
void ElasticWaveEquation3D::ElasticWaveEquation::algebraicSource(const double* const Q,double* S) {
  
  double rho,c_p,c_s;
  extractMaterialParameters(Q,rho,c_p,c_s);

  double jacobian;
  double q_x,q_y,q_z;
  double r_x,r_y,r_z;
  double s_x,s_y,s_z;
  extractMetricParameters(Q,jacobian,
			  q_x,q_y,q_z,
			  r_x,r_y,r_z,
			  s_x,s_y,s_z);
  
  double d_x,d_y,d_z;
  extractDampingParameters(Q,d_x,d_y,d_z);

  double alpha = 0.0;
  
  double alpha_x,alpha_y,alpha_z;
  alpha_x = alpha+0.0*d_x;
  alpha_y = alpha+0.0*d_y;
  alpha_z = alpha+0.0*d_z;
  
  double mu     = rho*c_s*c_s;
  double lambda = rho*c_p*c_p-2*mu;
  double rho_inv= 1.0/(rho);

  double Px[9],Py[9],Pz[9];
  extractAuxiliaryParameters(Q,Px,Py,Pz);

 

  S[0]=rho_inv * (Px[0] + Py[0] + Pz[0]);
  S[1]=rho_inv * (Px[1] + Py[1] + Pz[1]);
  S[2]=rho_inv * (Px[2] + Py[2] + Pz[2]);

  S[3]=((2*mu+lambda)*Px[3] + lambda*(Px[4] + Px[5]))
      +((2*mu+lambda)*Py[3] + lambda*(Py[4] + Py[5]))
      +((2*mu+lambda)*Pz[3] + lambda*(Pz[4] + Pz[5]));
  S[4]=((2*mu+lambda)*Px[4] + lambda*(Px[3] + Px[5]))
      +((2*mu+lambda)*Py[4] + lambda*(Py[3] + Py[5]))
      +((2*mu+lambda)*Pz[4] + lambda*(Pz[3] + Pz[5]));
  S[5]=((2*mu+lambda)*Px[5] + lambda*(Px[4] + Px[3]))
      +((2*mu+lambda)*Py[5] + lambda*(Py[4] + Py[3]))
      +((2*mu+lambda)*Pz[5] + lambda*(Pz[4] + Pz[3]));

  S[6]= mu*(Px[6] + Py[6] + Pz[6]);
  S[7]= mu*(Px[7] + Py[7] + Pz[7]);
  S[8]= mu*(Px[8] + Py[8] + Pz[8]);
    
  for (int j = 0; j < 9; j++){
    S[9+j]  = (d_x + alpha_x)*Px[j];
    S[18+j] = (d_y + alpha_y)*Py[j];
    S[27+j] = (d_z + alpha_z)*Pz[j];
  }
}

void  ElasticWaveEquation3D::ElasticWaveEquation::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // @todo Please implement/augment if required

  double u_q = gradQ[0];
  double v_q = gradQ[1];
  double w_q = gradQ[2];

  double u_r = gradQ[NumberOfVariables + 0];
  double v_r = gradQ[NumberOfVariables + 1];
  double w_r = gradQ[NumberOfVariables + 2];

  double u_s = gradQ[2*NumberOfVariables + 0];
  double v_s = gradQ[2*NumberOfVariables + 1];
  double w_s = gradQ[2*NumberOfVariables + 2];

  double jacobian;
  double q_x,q_y,q_z;
  double r_x,r_y,r_z;
  double s_x,s_y,s_z;

  extractMetricParameters(Q,jacobian,
			  q_x,q_y,q_z,
			  r_x,r_y,r_z,
			  s_x,s_y,s_z);

  for(int i = 0 ; i<3*NumberOfVariables ;i ++){
    BgradQ[i] = 0;
  }
  
  BgradQ[3] = -q_x*u_q;
  BgradQ[4] = -q_y*v_q;
  BgradQ[5] = -q_z*w_q;
  BgradQ[6] = -(q_y*u_q+q_x*v_q); //sigma_xy
  BgradQ[7] = -(q_z*u_q+q_x*w_q); //sigma_xz
  BgradQ[8] = -(q_z*v_q+q_y*w_q); //sigma_yz
  
  BgradQ[12] = BgradQ[3];
  BgradQ[13] = BgradQ[4];
  BgradQ[14] = BgradQ[5];
  BgradQ[15] = BgradQ[6];
  BgradQ[16] = BgradQ[7];
  BgradQ[17] = BgradQ[8];


  BgradQ[NumberOfVariables+3] = -r_x*u_r;
  BgradQ[NumberOfVariables+4] = -r_y*v_r;
  BgradQ[NumberOfVariables+5] = -r_z*w_r;
  BgradQ[NumberOfVariables+6] = -(r_y*u_r+r_x*v_r); //sigma_xy
  BgradQ[NumberOfVariables+7] = -(r_z*u_r+r_x*w_r); //sigma_xz
  BgradQ[NumberOfVariables+8] = -(r_z*v_r+r_y*w_r); //sigma_yz

  BgradQ[NumberOfVariables+21] = BgradQ[NumberOfVariables+3];
  BgradQ[NumberOfVariables+22] = BgradQ[NumberOfVariables+4];
  BgradQ[NumberOfVariables+23] = BgradQ[NumberOfVariables+5];
  BgradQ[NumberOfVariables+24] = BgradQ[NumberOfVariables+6];
  BgradQ[NumberOfVariables+25] = BgradQ[NumberOfVariables+7];
  BgradQ[NumberOfVariables+26] = BgradQ[NumberOfVariables+8];

  BgradQ[2*NumberOfVariables+3]  = -s_x*u_s;
  BgradQ[2*NumberOfVariables+4]  = -s_y*v_s;
  BgradQ[2*NumberOfVariables+5]  = -s_z*w_s;
  BgradQ[2*NumberOfVariables+6]  = -(s_y*u_s+s_x*v_s); //sigma_xy
  BgradQ[2*NumberOfVariables+7]  = -(s_z*u_s+s_x*w_s); //sigma_xz
  BgradQ[2*NumberOfVariables+8]  = -(s_z*v_s+s_y*w_s); //sigma_yz

  BgradQ[3*NumberOfVariables-6] = BgradQ[2*NumberOfVariables+3];
  BgradQ[3*NumberOfVariables-5] = BgradQ[2*NumberOfVariables+4];
  BgradQ[3*NumberOfVariables-4] = BgradQ[2*NumberOfVariables+5];
  BgradQ[3*NumberOfVariables-3] = BgradQ[2*NumberOfVariables+6];
  BgradQ[3*NumberOfVariables-2] = BgradQ[2*NumberOfVariables+7];
  BgradQ[3*NumberOfVariables-1] = BgradQ[2*NumberOfVariables+8];

}

void  ElasticWaveEquation3D::ElasticWaveEquation::pointSource(const double* const Q,const double* const x,const double t,const double dt, double* forceVector,int n) {
  
  static tarch::logging::Log _log("MyLinearWaveSolver::pointSource");
  double jacobian;
  double q_x,q_y,q_z;
  double r_x,r_y,r_z;
  double s_x,s_y,s_z;
  
  extractMetricParameters(Q,jacobian,
   			  q_x,q_y,q_z,
   			  r_x,r_y,r_z,
   			  s_x,s_y,s_z);

  
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  double f = 0.0;
  double M0 = 1000.0 / jacobian;

  double a_x = 0.0;
  double a_y = 0.0;
  double a_z = 0.0;
  
  double b_x = 15.0;
  double b_y = 15.0;
  double b_z = 15.0;


#if defined(_CURVILINEAR)
  double blockWidth_y = (b_y-a_y);
  double blockWidth_x = (b_x-a_x);
  double blockWidth_z = (b_z-a_z);
#else
  double blockWidth_y = 1.;
  double blockWidth_x = 1.;
  double blockWidth_z = 1.;
#endif  

  if(n == 0){
    
#if defined(_LOH1)
    t0=0.1;
    f = M0*t/(t0*t0)*std::exp(-t/t0);
    double x1 = 2.0;
    double y1 = 4.0;
    double z1 = 4.0;
#else
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    double x1 = 2.7;
    double y1 = 2.5;
    double z1 = 2.5;
#endif
    double fault_ref_x = (x1-a_x)/(blockWidth_x);
    double fault_ref_y = (y1-a_y)/(blockWidth_y);
    double fault_ref_z = (z1-a_z)/(blockWidth_z);
    
    x0[0] = fault_ref_x;
    x0[1] = fault_ref_y;
    x0[2] = fault_ref_z;

#if defined(LOH1)
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
    forceVector[4] = 0.0;
    forceVector[5] = 0.0;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 1.*f;
#else    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 1.*f;
    forceVector[4] = 1.*f;
    forceVector[5] = 1.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
#endif
    
    
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

    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.*f;
    forceVector[5] = 0.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
    
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
    
    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.*f;
    forceVector[5] = 0.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
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

    forceVector[0] = 0.0;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.*f;
    forceVector[4] = 0.*f;
    forceVector[5] = 0.*f;
    forceVector[6] = 0.0;
    forceVector[7] = 0.0;
    forceVector[8] = 0.0;
  }
}


    /**
     * @TODO LR : document
     */
void ElasticWaveEquation3D::ElasticWaveEquation::multiplyMaterialParameterMatrix(const double* const Q, double* rhs) {
  
  double rho,c_p,c_s;
  extractMaterialParameters(Q,rho,c_p,c_s);

  double jacobian;
  double q_x,q_y,q_z;
  double r_x,r_y,r_z;
  double s_x,s_y,s_z;
  extractMetricParameters(Q,jacobian,
			  q_x,q_y,q_z,
			  r_x,r_y,r_z,
			  s_x,s_y,s_z);
  
  double d_x,d_y,d_z;
  extractDampingParameters(Q,d_x,d_y,d_z);
  
  double mu     = rho*c_s*c_s;
  double lambda = rho*c_p*c_p-2*mu;
  double rho_j_inv= 1.0/(rho*jacobian);

  rhs[0]=rho_j_inv * rhs[0];
  rhs[1]=rho_j_inv * rhs[1];
  rhs[2]=rho_j_inv * rhs[2];
  
  double lam_temp = lambda * (rhs[3] + rhs[4] + rhs[5]);
  rhs[3]=(2*mu) * rhs[3] +lam_temp;
  rhs[4]=(2*mu) * rhs[4] +lam_temp;
  rhs[5]=(2*mu) * rhs[5] +lam_temp;

  rhs[6]= mu*rhs[6];
  rhs[7]= mu*rhs[7];
  rhs[8]= mu*rhs[8];


  for (int j = 0; j < 3; j++){
    rhs[9+j] = rhs[ 9+j]/jacobian;
    rhs[18+j]= rhs[18+j]/jacobian;
    rhs[27+j]= rhs[27+j]/jacobian;
  }

  for (int j = 0; j < 9; j++){
    rhs[9+j] = d_x*rhs[9+j];
    rhs[18+j]= d_y*rhs[18+j];
    rhs[27+j]= d_z*rhs[27+j];
  }

  rhs[NumberOfVariables+0]=rho_j_inv * rhs[NumberOfVariables+0];
  rhs[NumberOfVariables+1]=rho_j_inv * rhs[NumberOfVariables+1];
  rhs[NumberOfVariables+2]=rho_j_inv * rhs[NumberOfVariables+2];
  
  lam_temp = lambda * (rhs[NumberOfVariables+3] + rhs[NumberOfVariables+4] + rhs[NumberOfVariables+5]);

  rhs[NumberOfVariables+3]=(2*mu) * rhs[NumberOfVariables+3] +lam_temp;
  rhs[NumberOfVariables+4]=(2*mu) * rhs[NumberOfVariables+4] +lam_temp;
  rhs[NumberOfVariables+5]=(2*mu) * rhs[NumberOfVariables+5] +lam_temp;

  rhs[NumberOfVariables+6]= mu*rhs[NumberOfVariables+6];
  rhs[NumberOfVariables+7]= mu*rhs[NumberOfVariables+7];
  rhs[NumberOfVariables+8]= mu*rhs[NumberOfVariables+8];


  for (int j = 0; j < 3; j++){
    rhs[NumberOfVariables+9+j] = rhs[NumberOfVariables+9+j]/jacobian;
    rhs[NumberOfVariables+18+j]= rhs[NumberOfVariables+18+j]/jacobian;
    rhs[NumberOfVariables+27+j]= rhs[NumberOfVariables+27+j]/jacobian;
  }
  
  for (int j = 0; j < 9; j++){
    rhs[NumberOfVariables+9+j] = d_x*rhs[NumberOfVariables+9+j];
    rhs[NumberOfVariables+18+j]= d_y*rhs[NumberOfVariables+18+j];
    rhs[NumberOfVariables+27+j]= d_z*rhs[NumberOfVariables+27+j];
  }
  

  rhs[2*NumberOfVariables+0] = rho_j_inv * rhs[2*NumberOfVariables+0];
  rhs[2*NumberOfVariables+1] = rho_j_inv * rhs[2*NumberOfVariables+1];
  rhs[2*NumberOfVariables+2] = rho_j_inv * rhs[2*NumberOfVariables+2];
  
  lam_temp = lambda * (rhs[2*NumberOfVariables+3] + rhs[2*NumberOfVariables+4] + rhs[2*NumberOfVariables+5]);

  rhs[2*NumberOfVariables+3]=(2*mu) * rhs[2*NumberOfVariables+3] +lam_temp;
  rhs[2*NumberOfVariables+4]=(2*mu) * rhs[2*NumberOfVariables+4] +lam_temp;
  rhs[2*NumberOfVariables+5]=(2*mu) * rhs[2*NumberOfVariables+5] +lam_temp;

  rhs[2*NumberOfVariables+6]= mu*rhs[2*NumberOfVariables+6];
  rhs[2*NumberOfVariables+7]= mu*rhs[2*NumberOfVariables+7];
  rhs[2*NumberOfVariables+8]= mu*rhs[2*NumberOfVariables+8];

  for (int j = 0; j < 3; j++){
    rhs[2*NumberOfVariables+9+j] = rhs[2*NumberOfVariables+9+j] /jacobian;
    rhs[2*NumberOfVariables+18+j]= rhs[2*NumberOfVariables+18+j]/jacobian;
    rhs[2*NumberOfVariables+27+j]= rhs[2*NumberOfVariables+27+j]/jacobian;
  }
  
  for (int j = 0; j < 9; j++){
    rhs[2*NumberOfVariables+9+j] = d_x*rhs[2*NumberOfVariables+9+j];
    rhs[2*NumberOfVariables+18+j]= d_y*rhs[2*NumberOfVariables+18+j];
    rhs[2*NumberOfVariables+27+j]= d_z*rhs[2*NumberOfVariables+27+j];
  }
}

void  ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver(double* FL_,double* FR_, const double* const QL_,const double* const QR_,const double dt,const int normalNonZeroIndex, bool isBoundaryFace, int faceIndex){
  
  constexpr int NumberOfVariables  = ElasticWaveEquation::NumberOfVariables;
  constexpr int NumberOfVariables2 = NumberOfVariables*NumberOfVariables;
  constexpr int numberOfParameters = ElasticWaveEquation::NumberOfParameters;
  constexpr int numberOfData       = NumberOfVariables+numberOfParameters;
  constexpr int basisSize          = ElasticWaveEquation::Order+1;
  constexpr int order              = basisSize - 1;


#ifdef OPT_KERNELS
  double* FL = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getFFaceGenArraySize()]();
  double* FR = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getFFaceGenArraySize()]();
  double* QL = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getQFaceGenArraySize()]();
  double* QR = new double[ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::getQFaceGenArraySize()]();

  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_optimised2generic(FL_,FL);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_optimised2generic(FR_,FR);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::QFace_optimised2generic(QL_,QL);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::QFace_optimised2generic(QR_,QR);
#else
  double* FL=FL_;
  double* FR=FR_;
  const double* QL=QL_;
  const double* QR=QR_;
#endif  

  kernels::idx3 idx_QLR(basisSize,basisSize,numberOfData);
  kernels::idx3 idx_FLR(basisSize,basisSize,NumberOfVariables);

  double n[3]={0,0,0};
  n[normalNonZeroIndex]=1;

  double n_p[3]={0,0,0};
  double n_m[3]={0,0,0};

  double m_p[3]={0,0,0};
  double m_m[3]={0,0,0};

  double l_p[3]={0,0,0};  
  double l_m[3]={0,0,0};

  double norm_p_qr = 1.0;
  double norm_m_qr = 1.0;
  
  double FLn, FLm, FLl, FRn,FRm,FRl;
  double FL_n,FL_m,FL_l,FR_n,FR_m,FR_l;
  double FLx,FLy,FLz,FRx,FRy,FRz;
  double FL_x,FL_y,FL_z,FR_x,FR_y,FR_z;

  double jacobian;

  //for (int k = 0; k < 3; k++){

  // n_m[k] = n[k];
  //   n_p[k] = n[k];
  // }

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {

      double qm_x,qm_y,qm_z;
      double rm_x,rm_y,rm_z;
      double sm_x,sm_y,sm_z;
      extractMetricParameters(QL+idx_QLR(i,j,0),
			      jacobian,
			      qm_x,qm_y,qm_z,
			      rm_x,rm_y,rm_z,
			      sm_x,sm_y,sm_z);
      
      double qp_x,qp_y,qp_z;
      double rp_x,rp_y,rp_z;
      double sp_x,sp_y,sp_z;
      extractMetricParameters(QR+idx_QLR(i,j,0),
			      jacobian,
			      qp_x,qp_y,qp_z,
			      rp_x,rp_y,rp_z,
			      sp_x,sp_y,sp_z);


      double rho_m,cp_m,cs_m;
      extractMaterialParameters(QL+idx_QLR(i,j,0),rho_m,cp_m,cs_m);
      double rho_p,cp_p,cs_p;
      extractMaterialParameters(QR+idx_QLR(i,j,0),rho_p,cp_p,cs_p);

      double mu_p = cs_p*cs_p*rho_p;
      double lam_p= rho_p*cp_p*cp_p-2*mu_p;      
      
      double mu_m = cs_m*cs_m*rho_m;
      double lam_m= rho_m*cp_m*cp_m-2*mu_m;


      double dm_x,dm_y,dm_z;
      extractDampingParameters(QL+idx_QLR(i,j,0),dm_z,dm_y,dm_z);

      double dp_x,dp_y,dp_z;
      extractDampingParameters(QR+idx_QLR(i,j,0),dp_z,dp_y,dp_z);

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
      double zs_p=rho_p*cs_p;
      double zp_p=rho_p*cp_p;      
      double zs_m=rho_m*cs_m;
      double zp_m=rho_m*cp_m;
      
      // impedance must be greater than zero !
      if (zp_p <= 0.0 || zp_m <= 0.0){
	std::cout<< zp_p<<" "<<zp_m<<"\n";
	std::cout<<" Impedance must be greater than zero ! "<< std::endl;
	std::exit(-1);
      }

      // generate interface data preserving the amplitude of the outgoing charactertritics
      // and satisfying interface conditions exactly.
      double vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p;        
      double vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m;    

      if (isBoundaryFace) {
	//double r= faceIndex==10 ? 1 : 0;
#if defined(_LOH1)	
	double r= faceIndex==0 ? 1 : 0;
#else
	double r= 0;
#endif	
	riemannSolver_boundary(faceIndex,r,vn_m,vm_m,vl_m,Tn_m,Tm_m,Tl_m,zp_m,zs_m,vn_hat_m,vm_hat_m,vl_hat_m,Tn_hat_m,Tm_hat_m,Tl_hat_m);
	riemannSolver_boundary(faceIndex,r,vn_p,vm_p,vl_p,Tn_p,Tm_p,Tl_p,zp_p,zs_p,vn_hat_p,vm_hat_p,vl_hat_p,Tn_hat_p,Tm_hat_p,Tl_hat_p);      
      }else {
	riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
	riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
	riemannSolver_Nodal(vl_p,vl_m, Tl_p, Tl_m, zs_p , zs_m, vl_hat_p , vl_hat_m, Tl_hat_p, Tl_hat_m);
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
     
      // construct flux fluctuation vectors obeying the eigen structure of the PDE
      // and choose physically motivated penalties such that we can prove
      // numerical stability.

      // fluctiations for physical variables
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

      // for (int k = 9; k < NumberOfVariables; k++)
      // 	{
      // 	  FR_[idx_FLR(i,j, k)] = 0.0;
      // 	  FL_[idx_FLR(i,j, k)] = 0.0;
      // 	}


      // fluctiations for x-direction auxiliary variables
      FR[idx_FLR(i,j, 9)] = n_p[0]*dp_x*norm_p_qr*FRx;
      FL[idx_FLR(i,j, 9)] = n_m[0]*dm_x*norm_m_qr*FLx;
    
      FR[idx_FLR(i,j, 10)] = n_p[0]*dp_x*norm_p_qr*FRy;
      FL[idx_FLR(i,j, 10)] = n_m[0]*dm_x*norm_m_qr*FLy;

      FR[idx_FLR(i,j, 11)] = n_p[0]*dp_x*norm_p_qr*FRz;
      FL[idx_FLR(i,j, 11)] = n_m[0]*dm_x*norm_m_qr*FLz;

      FL[idx_FLR(i,j, 12)] = n_m[0]*dm_x*norm_m_qr*n_m[0]*FL_x;
      FL[idx_FLR(i,j, 13)] = n_m[0]*dm_x*norm_m_qr*n_m[1]*FL_y;
      FL[idx_FLR(i,j, 14)] = n_m[0]*dm_x*norm_m_qr*n_m[2]*FL_z;

      FR[idx_FLR(i,j, 12)] = -n_p[0]*dp_x*norm_p_qr*n_p[0]*FR_x;
      FR[idx_FLR(i,j, 13)] = -n_p[0]*dp_x*norm_p_qr*n_p[1]*FR_y;
      FR[idx_FLR(i,j, 14)] = -n_p[0]*dp_x*norm_p_qr*n_p[2]*FR_z;
    
      FL[idx_FLR(i,j, 15)] =  n_m[0]*dm_x*norm_m_qr*(n_m[1]*FL_x + n_m[0]*FL_y);
      FL[idx_FLR(i,j, 16)] =  n_m[0]*dm_x*norm_m_qr*(n_m[2]*FL_x + n_m[0]*FL_z);
      FL[idx_FLR(i,j, 17)] =  n_m[0]*dm_x*norm_m_qr*(n_m[2]*FL_y + n_m[1]*FL_z);

      FR[idx_FLR(i,j, 15)] = -n_p[0]*dp_x*norm_p_qr*(n_p[1]*FR_x + n_p[0]*FR_y);
      FR[idx_FLR(i,j, 16)] = -n_p[0]*dp_x*norm_p_qr*(n_p[2]*FR_x + n_p[0]*FR_z);
      FR[idx_FLR(i,j, 17)] = -n_p[0]*dp_x*norm_p_qr*(n_p[2]*FR_y + n_p[1]*FR_z);

      // fluctiations for y-direction auxiliary variables
      FR[idx_FLR(i,j, 18)] = n_p[1]*dp_y*norm_p_qr*FRx;
      FL[idx_FLR(i,j, 18)] = n_m[1]*dm_y*norm_m_qr*FLx;
    
      FR[idx_FLR(i,j, 19)] = n_p[1]*dp_y*norm_p_qr*FRy;
      FL[idx_FLR(i,j, 19)] = n_m[1]*dm_y*norm_m_qr*FLy;

      FR[idx_FLR(i,j, 20)] = n_p[1]*dp_y*norm_p_qr*FRz;
      FL[idx_FLR(i,j, 20)] = n_m[1]*dm_y*norm_m_qr*FLz;

      FL[idx_FLR(i,j, 21)] = n_m[1]*dm_y*norm_m_qr*n_m[0]*FL_x;
      FL[idx_FLR(i,j, 22)] = n_m[1]*dm_y*norm_m_qr*n_m[1]*FL_y;
      FL[idx_FLR(i,j, 23)] = n_m[1]*dm_y*norm_m_qr*n_m[2]*FL_z;

      FR[idx_FLR(i,j, 21)] = -n_p[1]*dp_y*norm_p_qr*n_p[0]*FR_x;
      FR[idx_FLR(i,j, 22)] = -n_p[1]*dp_y*norm_p_qr*n_p[1]*FR_y;
      FR[idx_FLR(i,j, 23)] = -n_p[1]*dp_y*norm_p_qr*n_p[2]*FR_z;
    
      FL[idx_FLR(i,j, 24)] =  n_m[1]*dm_y*norm_m_qr*(n_m[1]*FL_x + n_m[0]*FL_y);
      FL[idx_FLR(i,j, 25)] =  n_m[1]*dm_y*norm_m_qr*(n_m[2]*FL_x + n_m[0]*FL_z);
      FL[idx_FLR(i,j, 26)] =  n_m[1]*dm_y*norm_m_qr*(n_m[2]*FL_y + n_m[1]*FL_z);

      FR[idx_FLR(i,j, 24)] = -n_p[1]*dp_y*norm_p_qr*(n_p[1]*FR_x + n_p[0]*FR_y);
      FR[idx_FLR(i,j, 25)] = -n_p[1]*dp_y*norm_p_qr*(n_p[2]*FR_x + n_p[0]*FR_z);
      FR[idx_FLR(i,j, 26)] = -n_p[1]*dp_y*norm_p_qr*(n_p[2]*FR_y + n_p[1]*FR_z);


      // fluctiations for z-direction auxiliary variables
      FR[idx_FLR(i,j, 27)] = n_p[2]*dp_z*norm_p_qr*FRx;
      FL[idx_FLR(i,j, 27)] = n_m[2]*dm_z*norm_m_qr*FLx;
    
      FR[idx_FLR(i,j, 28)] = n_p[2]*dp_z*norm_p_qr*FRy;
      FL[idx_FLR(i,j, 28)] = n_m[2]*dm_z*norm_m_qr*FLy;

      FR[idx_FLR(i,j, 29)] = n_p[2]*dp_z*norm_p_qr*FRz;
      FL[idx_FLR(i,j, 29)] = n_m[2]*dm_z*norm_m_qr*FLz;

      FL[idx_FLR(i,j, 30)] = n_m[2]*dm_z*norm_m_qr*n_m[0]*FL_x;
      FL[idx_FLR(i,j, 31)] = n_m[2]*dm_z*norm_m_qr*n_m[1]*FL_y;
      FL[idx_FLR(i,j, 32)] = n_m[2]*dm_z*norm_m_qr*n_m[2]*FL_z;

      FR[idx_FLR(i,j, 30)] = -n_p[2]*dp_z*norm_p_qr*n_p[0]*FR_x;
      FR[idx_FLR(i,j, 31)] = -n_p[2]*dp_z*norm_p_qr*n_p[1]*FR_y;
      FR[idx_FLR(i,j, 32)] = -n_p[2]*dp_z*norm_p_qr*n_p[2]*FR_z;
    
      FL[idx_FLR(i,j, 33)] =  n_m[2]*dm_z*norm_m_qr*(n_m[1]*FL_x + n_m[0]*FL_y);
      FL[idx_FLR(i,j, 34)] =  n_m[2]*dm_z*norm_m_qr*(n_m[2]*FL_x + n_m[0]*FL_z);
      FL[idx_FLR(i,j, 35)] =  n_m[2]*dm_z*norm_m_qr*(n_m[2]*FL_y + n_m[1]*FL_z);

      FR[idx_FLR(i,j, 33)] = -n_p[2]*dp_z*norm_p_qr*(n_p[1]*FR_x + n_p[0]*FR_y);
      FR[idx_FLR(i,j, 34)] = -n_p[2]*dp_z*norm_p_qr*(n_p[2]*FR_x + n_p[0]*FR_z);
      FR[idx_FLR(i,j, 35)] = -n_p[2]*dp_z*norm_p_qr*(n_p[2]*FR_y + n_p[1]*FR_z);
    }    
  }

  #ifdef OPT_KERNELS
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_generic2optimised(FL,FL_);
  ElasticWaveEquation3D::ElasticWaveEquation_kernels::aderdg::converter::FFace_generic2optimised(FR,FR_);

  delete[] FL;  
  delete[] FR;  
  delete[] QL;  
  delete[] QR;
#endif  

}

//transformation parameter
void ElasticWaveEquation3D::ElasticWaveEquation::extractMetricParameters(const double* const Q,
									 double& jacobian,
									 double& q_x, double&q_y, double& q_z,
									 double& r_x, double&r_y, double& r_z,
									 double& s_x, double&s_y, double& s_z) {
  jacobian = Q[42];
  q_x = Q[43];
  q_y = Q[44];
  q_z = Q[45];
  r_x = Q[46];
  r_y = Q[47];
  r_z = Q[48];
  s_x = Q[49];
  s_y = Q[50];
  s_z = Q[51];
}

void ElasticWaveEquation3D::ElasticWaveEquation::extractMaterialParameters(const double* const Q,
									 double& rho, double& cp, double &cs){
  rho= Q[36];
  cp = Q[37];
  cs = Q[38];
}
  
//pml damping
void ElasticWaveEquation3D::ElasticWaveEquation::extractDampingParameters(const double* const Q,
									 double& d_x, double& d_y, double& d_z){
  d_x = Q[39];
  d_y = Q[40];
  d_z = Q[41];
}


void ElasticWaveEquation3D::ElasticWaveEquation::extractVariables(const double* const Q,
								  double& u,  double& v,  double& w,
  								  double& sigma_xx,double& sigma_yy,double& sigma_zz,
								  double& sigma_xy,double& sigma_xz,double& sigma_yz){
  u=Q[0];
  v=Q[1];
  w=Q[2];
  sigma_xx=Q[3];
  sigma_yy=Q[4];
  sigma_zz=Q[5];  
  sigma_xy=Q[6];
  sigma_xz=Q[7];
  sigma_yz=Q[8];    
}


void ElasticWaveEquation3D::ElasticWaveEquation::extractAuxiliaryParameters(const double* const Q,
									    double* Px,double* Py,double* Pz){
  for (int j = 0; j < 9; j++){
    Px[j] = Q[9+j];
    Py[j] = Q[18+j];
    Pz[j] = Q[27+j];
  }
}

//Gram Schmidt orthonormalization
void  ElasticWaveEquation3D::ElasticWaveEquation::Gram_Schmidt(double* y, double* z){
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
  }
  
  double norm_z = std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  
  for (int i = 0; i< 3; i++){
    z[i] =  z[i]/norm_z;
  }
}

void  ElasticWaveEquation3D::ElasticWaveEquation::localBasis(double* n, double * m, double* l, int d){
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



void  ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
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

void  ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
   double p = 0.5*(z*v + sigma);
   if(z > 0){
     v_hat = (1+r)/z*p;
     sigma_hat = (1-r)*p;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }
}

void  ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
   double q = 0.5*(z*v - sigma);
   if(z > 0){
     v_hat = (1+r)/z*q;
     sigma_hat = -(1-r)*q;
   }else{
     v_hat = v;
     sigma_hat = sigma;
   }
}

void  ElasticWaveEquation3D::ElasticWaveEquation::get_normals(int normalNonZeroIndex,double& norm, double* n,const double* Q){

  double q_x;
  double q_y;
  double q_z;  
  double r_x;
  double r_y;
  double r_z;
  double s_x;
  double s_y;
  double s_z;
  double jacobian;
  
  extractMetricParameters(Q,jacobian,q_x,q_y,q_z,r_x,r_y,r_z,s_x,s_y,s_z);
  
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

void  ElasticWaveEquation3D::ElasticWaveEquation::extract_tractions_and_particle_velocity(double* n,const double* Q, double& Tx,double& Ty,double& Tz,double& vx,double& vy,double& vz ){
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

void  ElasticWaveEquation3D::ElasticWaveEquation::rotate_into_orthogonal_basis(double* n,double* m,double* l, double Tx,double Ty,double Tz, double& Tn, double& Tm, double& Tl){
    Tn= Tx*n[0] + Ty*n[1] + Tz*n[2];
    Tm= Tx*m[0] + Ty*m[1] + Tz*m[2];
    Tl= Tx*l[0] + Ty*l[1] + Tz*l[2];
}

void  ElasticWaveEquation3D::ElasticWaveEquation::rotate_into_physical_basis(double* n,double* m,double* l, double Fn,double Fm,double Fl, double& Fx, double& Fy, double& Fz){
  Fx = n[0]*Fn + m[0]*Fm + l[0]*Fl;
  Fy = n[1]*Fn + m[1]*Fm + l[1]*Fl;
  Fz = n[2]*Fn + m[2]*Fm + l[2]*Fl;
}

void  ElasticWaveEquation3D::ElasticWaveEquation::generate_fluctuations_left(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) + (T-T_hat));
}

void  ElasticWaveEquation3D::ElasticWaveEquation::generate_fluctuations_right(double z,  double T,double T_hat,double v, double v_hat, double& F){
  F = 0.5*(z*(v-v_hat) - (T-T_hat));
}

void  ElasticWaveEquation3D::ElasticWaveEquation::riemannSolver_boundary(int faceIndex,double r, double vn , double vm , double vl, double Tn , double Tm ,double Tl , double zp, double zs , double& vn_hat , double& vm_hat ,double& vl_hat , double& Tn_hat , double& Tm_hat ,double& Tl_hat)
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

