#include "kernels/KernelUtils.h"
#include "CurvilinearTransformation.h"
#if defined(_GLL)
#include "kernels/GaussLobattoBasis.h"
#else
#include "kernels/GaussLegendreBasis.h"
#endif

#if defined(USE_ASAGI)
#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "reader/asagi_reader.h"

double topography_fromASAGI(double x, double z, double* topography, easi::ArraysAdapter& adapter, easi::Component* model){
  double topo;
  easi::Query query(1,3);
  query.x(0,0)=x;
  query.x(0,1)=z;
  query.x(0,2)=0;
  model->evaluate(query,adapter);
  topo = -topography[0];
  return topo;
}
#endif



// void Linear::CurvilinearTransformation::test(){
//   return;
// }

CurvilinearTransformation::CurvilinearTransformation(const int a_num_nodes,const double a_dx ):
  num_nodes(a_num_nodes),
  dx(a_dx),
  fault_position(8.),
  a_x(0.),
  a_y(0.),
  a_z(0.),
  b_x(15.),
  b_y(15.),
  b_z(15.)
  //  b_x(10.),
  //  b_y(10.),
  //  b_z(10.)
{
  int ne_x = std::round(1/dx); //number of elements in x direction on whole domain
  int ne_y = ne_x;
  int ne_z = ne_x;
  
  int nx = ne_x *(num_nodes-1) + 1; //global number of nodes in x direction considering collocation
  int ny = ne_y *(num_nodes-1) + 1; //global number of nodes in y direction considering collocation
  int nz = ne_z *(num_nodes-1) + 1; //global number of nodes in z direction considering collocation


  double block_width_x;
  double block_width_y=b_y-a_y;
  double block_width_z=b_z-a_z;

  double recWidth_x = b_x-a_x; //width in x of the rectangle

  //relative position of the fault on the unit square
  double fault_ref = (fault_position-a_x)/(recWidth_x);


  for(int n=0 ; n<2 ; ++n ){
    ne_x = std::round(1/dx);
    
    if(n == 0){
      block_width_x=(fault_position-a_x);
      ne_x = std::round((ne_x+1)*fault_ref);
      nx =  ne_x *(num_nodes-1)+1;
    }else{
      double ne_x0 = std::round((ne_x+1)*fault_ref);
      ne_x = std::round(1/dx)-ne_x0;
      nx =  ne_x *(num_nodes-1)+1;
      block_width_x= (b_x-fault_position);
    }

    left_bnd_x[n]= new double[ny*nz];
    left_bnd_y[n]= new double[ny*nz];
    left_bnd_z[n]= new double[ny*nz];  

    right_bnd_x[n]= new double[ny*nz];
    right_bnd_y[n]= new double[ny*nz];
    right_bnd_z[n]= new double[ny*nz];  

    bottom_bnd_x[n]= new double[nx*nz];
    bottom_bnd_y[n]= new double[nx*nz];
    bottom_bnd_z[n]= new double[nx*nz];  

    top_bnd_x[n]= new double[nx*nz];
    top_bnd_y[n]= new double[nx*nz];
    top_bnd_z[n]= new double[nx*nz];  

    front_bnd_x[n]= new double[nx*ny];
    front_bnd_y[n]= new double[nx*ny];
    front_bnd_z[n]= new double[nx*ny];  

    back_bnd_x[n]= new double[nx*ny];
    back_bnd_y[n]= new double[nx*ny];
    back_bnd_z[n]= new double[nx*ny];

    getBoundaryCurves3D_cutOffTopography_withFault( num_nodes,
						    nx,ny,nz,n,fault_position,
						    block_width_x,  block_width_y ,  block_width_z ,
						    left_bnd_x[n],  left_bnd_y[n],  left_bnd_z[n],
						    right_bnd_x[n],  right_bnd_y[n],  right_bnd_z[n],
						    bottom_bnd_x[n],  bottom_bnd_y[n],  bottom_bnd_z[n],
						    top_bnd_x[n],  top_bnd_y[n],  top_bnd_z[n],
						    front_bnd_x[n],  front_bnd_y[n],  front_bnd_z[n],
						    back_bnd_x[n],  back_bnd_y[n],  back_bnd_z[n]);
  }


  unif_mesh= new double[num_nodes];
  for(int i = 0; i< num_nodes ; i++){
    unif_mesh[i]=i*1.0/(num_nodes-1);
  }
  
  denominator_lagrange= new double[num_nodes];
  //initalize Lagrange Denominator
  for(int i=0; i< num_nodes ; i++){
    double temp = 1.0;

    for (int j = 0 ; j< i ; j ++){
      temp = temp* (unif_mesh[i]-unif_mesh[j]);
    }
    
    for (int j = i+1 ; j < num_nodes ; j ++){
      temp = temp*(unif_mesh[i]-unif_mesh[j]);
    }
    denominator_lagrange[i] = 1.0/temp;
  }


  
  kernels::idx2 id_xy(num_nodes,num_nodes); //nodes,polynome

  lagrange_basis_at_nodes = new double[num_nodes*num_nodes];
#if defined(_GLL)
  kernels::initGaussLobattoNodesAndWeights(std::set<int>()); //empty set as it is not used in function
#else  
  kernels::initGaussLegendreNodesAndWeights(std::set<int>()); //empty set as it is not used in function
#endif  

  for(int i=0; i< num_nodes ; i++){
    for (int j = 0 ; j< num_nodes ; j ++){
#if defined(_GLL)
      lagrange_basis_at_nodes[id_xy(j,i)]
	=lagrangeBasis(kernels::lobatto::nodes[num_nodes-1][num_nodes-1-j],i,num_nodes);
#else
      lagrange_basis_at_nodes[id_xy(j,i)]
	=lagrangeBasis(kernels::legendre::nodes[num_nodes-1][j],i,num_nodes);
#endif
    }
  }




}

void CurvilinearTransformation::genCoordinates(const tarch::la::Vector<DIMENSIONS,double>& center,
					       const tarch::la::Vector<DIMENSIONS,double>& dx,    
					       double* gl_vals_x,double* gl_vals_y,double* gl_vals_z,
					       double* jacobian,
					       double* q_x,double* q_y,double* q_z,
					       double* r_x,double* r_y,double* r_z,
					       double* s_x,double* s_y,double* s_z){

  int ne_x = std::round(1.0/this->dx); //number of elements in x direction on whole domain
  int ne_y = std::round(1.0/this->dx); //number of elements in y direction on whole domain
  int ne_z = std::round(1.0/this->dx); //number of elements in z direction on whole domain       

  int nx = ne_x *(num_nodes-1) + 1; //global number of nodes in x direction considering collocation
  int ny = ne_y *(num_nodes-1) + 1; //global number of nodes in y direction considering collocation
  int nz = ne_z *(num_nodes-1) + 1; //global number of nodes in z direction considering collocation
  
  double block_width_x;
  double block_width_y=b_y-a_y;
  double block_width_z=b_z-a_z;

  double recWidth_x = b_x-a_x; //width in x of the rectangle

  //relative position of the fault on the unit square
  double fault_ref = (fault_position-a_x)/(recWidth_x);

  double offset_x=center[0]-0.5*dx[0];
  double offset_y=center[1]-0.5*dx[1];
  double offset_z=center[2]-0.5*dx[2];

  //defines which block we are in
  int n = getBlock(center,dx);
  
  double width_x=dx[0];
  double width_y=dx[1];
  double width_z=dx[2];


  //indeces of the first node within the element 
  int i_m; 
  int j_m; 
  int k_m;  


  if(n == 0){
    block_width_x=(fault_position-a_x);
    ne_x = std::round((ne_x+1)*fault_ref);
    nx =  ne_x *(num_nodes-1)+1;

    i_m =  std::round((offset_x)/width_x) *(num_nodes-1);
    j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    k_m =  std::round((offset_z)/width_z) *(num_nodes-1);    
  }else{
    double ne_x0 = std::round((ne_x+1)*fault_ref);
    ne_x = std::round(1/dx[0])-ne_x0;
    nx =  ne_x *(num_nodes-1)+1;
    block_width_x= (b_x-fault_position);

    i_m =  std::floor((offset_x-fault_ref)/width_x) *(num_nodes-1);
    j_m =  std::round((offset_y)/width_y) *(num_nodes-1);
    k_m =  std::round((offset_z)/width_z) *(num_nodes-1);    
  }

  
  // getBoundaryCurves3D( num_nodes,
  //            offset_x,  offset_y,  offset_z,
  //            width_x,  width_y ,  width_z ,
  //            left_bnd_x,  left_bnd_y,  left_bnd_z,
  //            right_bnd_x,  right_bnd_y,  right_bnd_z,
  //            bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
  //            top_bnd_x,  top_bnd_y,  top_bnd_z,
  //            front_bnd_x,  front_bnd_y,  front_bnd_z,
  //            back_bnd_x,  back_bnd_y,  back_bnd_z);

  

  
  // getBoundaryCurves3D_fixedTopFace_forBlock( num_nodes,
  //                nx,ny,nz,n  ,       
  //                block_width_x,  block_width_y ,  block_width_z ,
  //                left_bnd_x,  left_bnd_y,  left_bnd_z,
  //                right_bnd_x,  right_bnd_y,  right_bnd_z,
  //                bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
  //                top_bnd_x,  top_bnd_y,  top_bnd_z,
  //                front_bnd_x,  front_bnd_y,  front_bnd_z,
  //                back_bnd_x,  back_bnd_y,  back_bnd_z);

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
              left_bnd_x[n],
              right_bnd_x[n],
              bottom_bnd_x[n],
              top_bnd_x[n],
              front_bnd_x[n],
              back_bnd_x[n],
              curvilinear_x
              );

  transFiniteInterpolation3D( nx,  ny,  nz,
              k_m,  k_p ,
              j_m,  j_p ,
              i_m,  i_p ,
              num_nodes,
            width_x,width_y,width_z,
              left_bnd_y[n],
              right_bnd_y[n],
              bottom_bnd_y[n],
              top_bnd_y[n],
              front_bnd_y[n],
              back_bnd_y[n],
              curvilinear_y
              );

  transFiniteInterpolation3D( nx,  ny,  nz,
              k_m,  k_p ,
              j_m,  j_p ,
              i_m,  i_p ,
              num_nodes,
            width_x,width_y,width_z,            
              left_bnd_z[n],
              right_bnd_z[n],
              bottom_bnd_z[n],
              top_bnd_z[n],
              front_bnd_z[n],
              back_bnd_z[n],
              curvilinear_z
              );

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

double CurvilinearTransformation::fault(double y, double z,double a_y, double b_y, double a_z, double b_z){
  double pi = 3.14159265359;
  double angle = 60.0/360.0 * 2 * pi;
  double depth_y = 1.0; //b_y-a_y;

  double Ly = b_y-a_y;
  double Lz
    = b_z-a_z;

  double fault_surface;
  
  fault_surface=0.25*(std::sin(2*pi*y/Ly)*std::cos(2*pi*y/Ly))*std::sin(2*pi*z/Lz)*std::cos(2*pi*z/Lz);
  //  fault_surface=0;
  return  fault_surface;
   //return y*(1-y);
}


double CurvilinearTransformation::topography(double x, double z,double a_x, double b_x,double a_z, double b_z, double depth){
  double pi = 3.14159265359;
  double angle = 60.0/360.0 * 2 * pi;
  double Lx = b_x-a_x;
  double Lz = b_z-a_z;
  double topo;

  topo = 1.0*(0.1*(x + z) + 0.25*depth*(std::sin(4*pi*x/Lx+3.34)*std::cos(4*pi*x/Lx)
					* std::sin(4*pi*z/Lz+3.34)*std::cos(4*pi*z/Lz)));
  //  topo = 0.;
  return topo;
}  


void CurvilinearTransformation::getBoundaryCurves3D(int num_points,
			 double offset_x, double offset_y, double offset_z,
			 double width_x, double width_y , double width_z ,
			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
  double pi = 3.14159265359;

  int num_elt_x= std::ceil(1/width_x);
  int num_elt_y= std::ceil(1/width_y);
  int num_elt_z= std::ceil(1/width_z);

  int nx = num_elt_x*(num_points-1) + 1;
  int ny = num_elt_y*(num_points-1) + 1;
  int nz = num_elt_z*(num_points-1) + 1;  
  
  double dx= width_x/(num_points-1);
  double dy= width_y/(num_points-1);
  double dz= width_z/(num_points-1);

  kernels::idx2 id_xy(ny,nx);// back front
  kernels::idx2 id_xz(nz,nx);// bottom top
  kernels::idx2 id_yz(nz,ny);//left right


  //top bottom
  for(int k = 0 ; k < nz ; k ++){
    for(int i = 0 ; i < nx ; i ++){
      bottom_bnd_x[id_xz(k,i)] = dx*i;
      bottom_bnd_z[id_xz(k,i)] = dz*k;
      bottom_bnd_y[id_xz(k,i)] = 1;	    

      top_bnd_x[id_xz(k,i)] = dx*i;
      top_bnd_z[id_xz(k,i)] = dz*k;
      top_bnd_y[id_xz(k,i)] = 0+0.25*std::sin(2*pi*top_bnd_x[id_xz(k,i)])*std::sin(2*pi*top_bnd_z[id_xz(k,i)]);

	}
      }
    
  
  //left right
  for(int k = 0 ; k < nz ; k ++){
    for(int j = 0 ; j < ny ; j ++){
      left_bnd_x[id_xz(k,j)]  = 0;	    
      left_bnd_y[id_xz(k,j)]  = dy*j;
      left_bnd_z[id_xz(k,j)]  = dz*k;
      
      right_bnd_x[id_xz(k,j)] = 1.0;
      right_bnd_y[id_xz(k,j)] = dy*j;
      right_bnd_z[id_xz(k,j)] = dz*k;
      
    }
  }
  
  //front back
  for(int j = 0 ; j < ny ; j ++){
    for(int i = 0 ; i < nx ; i ++){
      
      front_bnd_x[id_xy(j,i)]  = dx*i;
      front_bnd_y[id_xy(j,i)]  = dy*j;
      front_bnd_z[id_xy(j,i)]  = 0;
      
      back_bnd_x[id_xy(j,i)] = dx*i;
      back_bnd_y[id_xy(j,i)] = dy*j;
      back_bnd_z[id_xy(j,i)] = 1;
    }
  }
    
}

void CurvilinearTransformation::getBoundaryCurves3D_fixedTopFace(int num_points,
				      double offset_x, double offset_y, double offset_z,
				      double width_x, double width_y , double width_z ,
				      double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
				      double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
				      double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
				      double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
				      double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
				      double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){
  
  double pi = 3.14159265359;


  int nx = std::round(1/width_x)*(num_points-1) + 1;
  int ny = std::round(1/width_y)*(num_points-1) + 1;
  int nz = std::round(1/width_z)*(num_points-1) + 1;  

  double dx= 1.0/(nx-1);
  double dz= 1.0/(nz-1);
  double dy= 1.0/(ny-1);

  double depth=1.0;

  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); //left right
  
  //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){
      top_bnd_x[id_xz(k,i)] = 0+dx*i;      
      top_bnd_z[id_xz(k,i)] = 0+dz*k;
      top_bnd_y[id_xz(k,i)] = 0- 0*(0.1*(top_bnd_x[id_xz(k,i)] + top_bnd_z[id_xz(k,i)])
			       + 0.25*(std::sin(4*pi*top_bnd_x[id_xz(k,i)]+3.34)*std::cos(4*pi*top_bnd_x[id_xz(k,i)])
				       * std::sin(4*pi*top_bnd_z[id_xz(k,i)]+3.34)*std::cos(4*pi*top_bnd_z[id_xz(k,i)])));

    }
  }
  
  
  for(int i = 0 ; i< nx; i++){
    for(int k = 0 ; k< nz; k++){
      bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
      bottom_bnd_y[id_xz(k,i)] = depth;      
    }
  }


  getInterpolatedFace_fromBottomAndTop(nx,ny,nz,0,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       left_bnd_x,left_bnd_y,left_bnd_z);
    
  getInterpolatedFace_fromBottomAndTop(nx,ny,nz,1,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       right_bnd_x,right_bnd_y,right_bnd_z);
  
  getInterpolatedFace_fromBottomAndTop(nz,ny,nx,2,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       front_bnd_x,front_bnd_y,front_bnd_z);

  getInterpolatedFace_fromBottomAndTop(nz,ny,nx,3,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       back_bnd_x,back_bnd_y,back_bnd_z);
}


//TODO: not working with new width definition !
void CurvilinearTransformation::getBoundaryCurves3D_fixedTopFace_forBlock(int num_points,
					       int nx, int ny, int nz, int n,
					       double width_x, double width_y , double width_z,	      
					       double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
					       double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
					       double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					       double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					       double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
					       double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
  double pi = 3.14159265359;


  double dx= 10*width_x/(nx-1);
  double dz= 10*width_y/(nz-1);
  double dy= 10*width_z/(ny-1);

  double depth=b_y;

  double x;
  double z;
  double y;

  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); // left right

  
  //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){

      top_bnd_y[id_xz(k,i)] = 0;      
      top_bnd_z[id_xz(k,i)] = 0+dz*k;

     
      
      if(n == 0){
      top_bnd_x[id_xz(k,i)] = 0+dx*i;
      
      }else{
	top_bnd_x[id_xz(k,i)] = 5 + dx*i;
	
      }

      x = top_bnd_x[id_xz(k,i)];
      z = top_bnd_z[id_xz(k,i)];

      top_bnd_y[id_xz(k,i)] -= topography(x, z, a_x, b_x, a_z, b_z,  depth);
    }
  }
  
  for(int k = 0 ; k< nz; k++){
    for(int i = 0 ; i< nx; i++){
      bottom_bnd_y[id_xz(k,i)] = depth;      
      bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
      bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      
    }
  }


  int n_left_right=ny;

  
  int n_block;
  int n_top_bottom;

  //left face
  n_block = nx;
  n_top_bottom = nz;

                                     //n_block, n_left_right, n_top_bottom,
  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,0,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       left_bnd_x,left_bnd_y,left_bnd_z);

  //right face
  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,1,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       right_bnd_x,right_bnd_y,right_bnd_z);


  //front face
  n_block = nz;
  n_top_bottom = nx;
  
  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,2,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       front_bnd_x,front_bnd_y,front_bnd_z);

  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,3,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       back_bnd_x,back_bnd_y,back_bnd_z);

  
}




void CurvilinearTransformation::getBoundaryCurves3D_cutOffTopography_withFault(int num_points,
				      //				      double offset_x, double offset_y, double offset_z,
				      //				      double width_x, double width_y , double width_z ,
						    int nx, int ny, int nz, int n,double fault_position,
						    double width_x, double width_y , double width_z,	      
						    double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
						    double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
						    double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
						    double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
						    double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
						    double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

  
  double pi = 3.14159265359;


  double dx= width_x/(nx-1);
  double dz= width_z/(nz-1);
  double dy= width_y/(ny-1);

  double depth= b_y;

  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); // left right


  double* X1 = new double[ny*nz];
  double* Y1 = new double[ny*nz];
  double* Z1 = new double[ny*nz];


  double x;
  double z;
  double y;

#if defined(USE_ASAGI)
  constexpr int chunksize = 1;
  AsagiReader asagiReader("");
  easi::YAMLParser parser(3,&asagiReader);
  easi::Component* model = parser.parse("topography.yaml");
  static double topography[chunksize];
  static easi::ArraysAdapter adapter;
  //Easi binding point for topography
  adapter.addBindingPoint("z",topography);
#endif

  
  //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){
      
      top_bnd_y[id_xz(k,i)] = a_y;      
      top_bnd_z[id_xz(k,i)] = a_z+dz*k;
         
      if(n == 0){
	
	top_bnd_x[id_xz(k,i)] = a_x+dx*i;
      
      }else{
	
	top_bnd_x[id_xz(k,i)] = fault_position + dx*i;
	
      }

      x = top_bnd_x[id_xz(k,i)];
      z = top_bnd_z[id_xz(k,i)];
      
#if !defined(USE_ASAGI)      
      top_bnd_y[id_xz(k,i)] -= topography(x, z, a_x, b_x, a_z, b_z,  depth);
#else
      top_bnd_y[id_xz(k,i)] -= topography_fromASAGI(x,z,topography, adapter, model);
#endif      
    }
  }
  
  //fault(double y, double z,double depth_y,double a_z, double b_z)
  //given bottom surface
  for(int k = 0 ; k< nz; k++){
    for(int i = 0 ; i< nx; i++){
      
      bottom_bnd_y[id_xz(k,i)] = depth;      
      bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
      bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      // if(n == 0) {
      // 	std::cout<<bottom_bnd_x[id_xz(k, i)]<<std::endl;
      // }
      
    }
  }  
  
  //and the fault surface
  for(int k = 0 ; k< nz; k++){
    for(int j = 0 ; j< ny; j++){
      
      left_bnd_y[id_yz(k,j)] = a_y + dy*j;
      left_bnd_z[id_yz(k,j)] = a_z +dz*k;

      if(n==0){
  	left_bnd_x[id_yz(k,j)] = a_x;
      }else{
	
	y = left_bnd_y[id_yz(k,j)];
	z = left_bnd_z[id_yz(k,j)];
	
  	left_bnd_x[id_yz(k,j)] = fault_position;

	// synthetic fault not compatible with computational geometry
	Y1[id_yz(k,j)] = a_y-(b_y-a_y)*0.3 + 1.3*dy*j;
	Z1[id_yz(k,j)] = a_z + dz*k;

	X1[id_yz(k,j)] = fault_position - fault(Y1[id_yz(k,j)], Z1[id_yz(k,j)], a_y, b_y, a_z, b_z);
	
      }
      
    }
  }

  //and the fault surface
  for(int k = 0 ; k< nz; k++){
    for(int j = 0 ; j< ny; j++){
      right_bnd_y[id_yz(k,j)] = a_y + dy*j;
      right_bnd_z[id_yz(k,j)] = a_z + dz*k;

      if(n==0){
	
	y = right_bnd_y[id_yz(k,j)];
	z = right_bnd_z[id_yz(k,j)];
	
	right_bnd_x[id_yz(k,j)] = fault_position;

	// synthetic fault not compatible with computational geometry
	Y1[id_yz(k,j)] = a_y-(b_y-a_y)*0.3 + 1.3*dy*j;
	Z1[id_yz(k,j)] = a_z + dz*k;
	
	X1[id_yz(k,j)] = fault_position - fault(Y1[id_yz(k,j)], Z1[id_yz(k,j)], a_y, b_y, a_z, b_z);
      }else{
	right_bnd_x[id_yz(k,j)] = b_x;
      }
    }
  }

  double g;

  for(int k = 0 ; k< nz; k++){
    for(int j = 0 ; j< ny; j++){
      if (n == 0){

	y = right_bnd_y[id_yz(k,j)];
	z = right_bnd_z[id_yz(k,j)];
	
	g =interpolate_fault_surface(X1, Y1, Z1, y, z, ny, nz);

	right_bnd_x[id_yz(k,j)] = g;

	
      }

      if (n == 1){

	y = left_bnd_y[id_yz(k,j)];
	z = left_bnd_z[id_yz(k,j)];
	
	g =interpolate_fault_surface(X1, Y1, Z1, y, z, ny, nz);

	left_bnd_x[id_yz(k,j)] = g;
	
      }

    }

  }

  {
    double* top_edge_x = new double[nz];
    double* bottom_edge_x = new double[nz];
    
    double* left_edge_x = new double[nx];
    double* right_edge_x = new double[nx];

    if (n == 0) {
      double distance_x_left=(right_bnd_x[id_yz(0,ny-1)]-a_x)/(nx-1);
      double distance_x_right=(right_bnd_x[id_yz(nz-1,ny-1)]-a_x)/(nx-1);
      
      for (int k = 0; k < nz; k++){
	
	top_edge_x[k] = right_bnd_x[id_yz(k,ny-1)];
	bottom_edge_x[k] = a_x;
	
      }
      
      for (int i = 0; i < nx; i++){
	
	left_edge_x[i] = a_x + i*distance_x_left;
	right_edge_x[i] = a_x + i*distance_x_right;
	
	
      }
      
      transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, bottom_bnd_x);
      
      
      distance_x_left=(right_bnd_x[id_yz(0,0)]-a_x)/(nx-1);
      distance_x_right=(right_bnd_x[id_yz(nz-1,0)]-a_x)/(nx-1);
      
      for (int k = 0; k < nz; k++){
	
	top_edge_x[k] = right_bnd_x[id_yz(k,0)];
	bottom_edge_x[k] = a_x;
	
      }
      
      for (int i = 0; i < nx; i++){
	
	left_edge_x[i] = a_x + i*distance_x_left;
	right_edge_x[i] = a_x + i*distance_x_right;
	
	
      }
      
      transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, top_bnd_x);
      
      
    }
    
    if (n == 1) {
      double distance_x_left=(b_x - left_bnd_x[id_yz(0,ny-1)])/(nx-1);
      double distance_x_right=(b_x - left_bnd_x[id_yz(nz-1,ny-1)])/(nx-1);
      
      for (int k = 0; k < nz; k++){
	
	top_edge_x[k] = b_x;
	bottom_edge_x[k] = left_bnd_x[id_yz(k,ny-1)];
	
      }
      
      //
      
      for (int i = 0; i < nx; i++){
	
	left_edge_x[i] = left_bnd_x[id_yz(0,ny-1)] + i*distance_x_left;
	right_edge_x[i] =  left_bnd_x[id_yz(nz-1,ny-1)] + i*distance_x_right;
	
	
	
      }
      
      transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, bottom_bnd_x);
      
      
      distance_x_left=(b_x - left_bnd_x[id_yz(0,0)])/(nx-1);
      distance_x_right=(b_x - left_bnd_x[id_yz(nz-1,0)])/(nx-1);
      
      for (int k = 0; k < nz; k++){
	
	top_edge_x[k] = b_x;
	bottom_edge_x[k] = left_bnd_x[id_yz(k,0)];   
	
      }
      
      
      
      for (int i = 0; i < nx; i++){
	
	left_edge_x[i] = left_bnd_x[id_yz(0,0)] + i*distance_x_left;
	right_edge_x[i] =  left_bnd_x[id_yz(nz-1,0)] + i*distance_x_right;
	
      }
      
      transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, top_bnd_x);
      
    }
  }
   //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){
      
      x = top_bnd_x[id_xz(k,i)];
      z = top_bnd_z[id_xz(k,i)];
#if !defined(USE_ASAGI)      
      top_bnd_y[id_xz(k,i)] -= topography(x, z, a_x, b_x, a_z, b_z,  depth);
#else
      top_bnd_y[id_xz(k,i)] -= topography_fromASAGI(x,z,topography, adapter, model);
#endif      
    }
  }
  

  {
  double* top_edge_y = new double[nz];
  double* bottom_edge_y = new double[nz];

  double* left_edge_y = new double[ny];
  double* right_edge_y = new double[ny];


  if ( n == 0){
    double distance_y_left=(bottom_bnd_y[id_xz(0,nx-1)]-top_bnd_y[id_xz(0,nx-1)])/(ny-1);
    double distance_y_right=(bottom_bnd_y[id_xz(nz-1,nx-1)]-top_bnd_y[id_xz(nz-1,nx-1)])/(ny-1);
    

    for (int k = 0; k < nz; k++){

      top_edge_y[k] = top_bnd_y[id_xz(k,nx-1)];
      bottom_edge_y[k] = bottom_bnd_y[id_xz(k,nx-1)];
     
   }

    for (int j = 0; j < ny; j++){

      left_edge_y[j]  = top_bnd_y[id_xz(0,nx-1)]    + j*distance_y_left;
      right_edge_y[j] = top_bnd_y[id_xz(nz-1,nx-1)] + j*distance_y_right;
     
   }

    transFiniteInterpolation_singleCoordinate(ny, nz, top_edge_y, bottom_edge_y, left_edge_y, right_edge_y,  right_bnd_y);


    distance_y_left=(bottom_bnd_y[id_xz(0,0)]-top_bnd_y[id_xz(0,0)])/(ny-1);
    distance_y_right=(bottom_bnd_y[id_xz(nz-1,0)]-top_bnd_y[id_xz(nz-1,0)])/(ny-1);


    for (int k = 0; k < nz; k++){

      top_edge_y[k] = top_bnd_y[id_xz(k,0)];
      bottom_edge_y[k] = bottom_bnd_y[id_xz(k,0)];
     
   }

    for (int j = 0; j < ny; j++){

      left_edge_y[j]  = top_bnd_y[id_xz(0,0)]     + j*distance_y_left;
      right_edge_y[j] = top_bnd_y[id_xz(nz-1,0)]  + j*distance_y_right;
     
   }

    transFiniteInterpolation_singleCoordinate(ny, nz, top_edge_y, bottom_edge_y, left_edge_y, right_edge_y,  left_bnd_y);
  }
  

  if ( n == 1){
    double distance_y_left=(bottom_bnd_y[id_xz(0,0)]-top_bnd_y[id_xz(0,0)])/(ny-1);
    double distance_y_right=(bottom_bnd_y[id_xz(nz-1,0)]-top_bnd_y[id_xz(nz-1,0)])/(ny-1);
    
    for (int k = 0; k < nz; k++){

      top_edge_y[k] = top_bnd_y[id_xz(k,0)];
      bottom_edge_y[k] = bottom_bnd_y[id_xz(k,0)];
     
   }

    for (int j = 0; j < ny; j++){

      left_edge_y[j] = top_bnd_y[id_xz(0,0)] + j*distance_y_left;
      right_edge_y[j] = top_bnd_y[id_xz(nz-1,0)] + j*distance_y_right;
     
   }

    transFiniteInterpolation_singleCoordinate(ny, nz, top_edge_y, bottom_edge_y, left_edge_y, right_edge_y,  left_bnd_y);

    distance_y_left=(bottom_bnd_y[id_xz(0,nx-1)]-top_bnd_y[id_xz(0,nx-1)])/(ny-1);
    distance_y_right=(bottom_bnd_y[id_xz(nz-1,nx-1)]-top_bnd_y[id_xz(nz-1,nx-1)])/(ny-1);
    
    for (int k = 0; k < nz; k++){

      top_edge_y[k] = top_bnd_y[id_xz(k,nx-1)];
      bottom_edge_y[k] = bottom_bnd_y[id_xz(k,nx-1)];
     
   }

    for (int j = 0; j < ny; j++){

      left_edge_y[j] = top_bnd_y[id_xz(0,nx-1)] + j*distance_y_left;
      right_edge_y[j] = top_bnd_y[id_xz(nz-1,nx-1)] + j*distance_y_right;
     
   }

    transFiniteInterpolation_singleCoordinate(ny, nz, top_edge_y, bottom_edge_y, left_edge_y, right_edge_y,  right_bnd_y);
  }
  
 
  }


  double* top_edge_y = new double[nx];
  double* bottom_edge_y = new double[nx];

  double* left_edge_y = new double[ny];
  double* right_edge_y = new double[ny];

 double distance_y_left=(bottom_bnd_y[id_xz(0,0)]-top_bnd_y[id_xz(0,0)])/(ny-1);
 double distance_y_right=(bottom_bnd_y[id_xz(0,nx-1)]-top_bnd_y[id_xz(0,nx-1)])/(ny-1);
    
    for (int i = 0; i < nx; i++){

      top_edge_y[i] = top_bnd_y[id_xz(0,i)];
      bottom_edge_y[i] = bottom_bnd_y[id_xz(0,i)];

      
     
   }

    

    for (int j = 0; j < ny; j++){

      left_edge_y[j]  = top_bnd_y[id_xz(0,0)] + j*distance_y_left;
      right_edge_y[j] = top_bnd_y[id_xz(0,nx-1)] + j*distance_y_right;

      
     
   }

    

    transFiniteInterpolation_singleCoordinate(nx, ny, left_edge_y, right_edge_y, top_edge_y, bottom_edge_y, front_bnd_y);


    distance_y_left=(bottom_bnd_y[id_xz(nz-1,0)]-top_bnd_y[id_xz(nz-1,0)])/(ny-1);
    distance_y_right=(bottom_bnd_y[id_xz(nz-1,nx-1)]-top_bnd_y[id_xz(nz-1,nx-1)])/(ny-1);

    
    for (int i = 0; i < nx; i++){
      
      top_edge_y[i] = top_bnd_y[id_xz(nz-1,i)];
      bottom_edge_y[i] = bottom_bnd_y[id_xz(nz-1,i)];
     
   }

    for (int j = 0; j < ny; j++){

      left_edge_y[j] = top_bnd_y[id_xz(nz-1,0)] + j*distance_y_left;
      right_edge_y[j] = top_bnd_y[id_xz(nz-1,nx-1)] + j*distance_y_right;
      
   }


    transFiniteInterpolation_singleCoordinate(nx, ny, left_edge_y, right_edge_y,  top_edge_y, bottom_edge_y,   back_bnd_y);


    double* top_edge_x = new double[ny];
    double* bottom_edge_x = new double[ny];
    
    double* left_edge_x = new double[nx];
    double* right_edge_x = new double[nx];
    if (n == 0)
      {
    double distance_x_left=(right_bnd_x[id_yz(nz-1,0)]-0)/(nx-1);
    double distance_x_right=(right_bnd_x[id_yz(nz-1,ny-1)]-0)/(nx-1);



    for (int j = 0; j < ny; j++){

      bottom_edge_x[j] = a_x;
      top_edge_x[j] = right_bnd_x[id_yz(nz-1,j)];
     
   }


    for (int i = 0; i < nx; i++){


       left_edge_x[i] = a_x + i*distance_x_left;
      right_edge_x[i] = a_x + i*distance_x_right;

     
   }

    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x,  top_edge_x, left_edge_x, right_edge_x,  back_bnd_x);

    

    distance_x_left=(right_bnd_x[id_yz(0,0)]-0)/(nx-1);
    distance_x_right=(right_bnd_x[id_yz(0,ny-1)]-0)/(nx-1);


    for (int j = 0; j < ny; j++){

      top_edge_x[j] = right_bnd_x[id_yz(0,j)];
      bottom_edge_x[j] = a_x;

     
   }


    for (int i = 0; i < nx; i++){

      left_edge_x[i] = a_x + i*distance_x_left;
      right_edge_x[i] = a_x + i*distance_x_right;
 
     
   }
    
    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, front_bnd_x);
      }


    if (n == 1)
      {
    double distance_x_left=(b_x-left_bnd_x[id_yz(nz-1,0)])/(nx-1);
    double distance_x_right=(b_x-left_bnd_x[id_yz(nz-1,ny-1)])/(nx-1);



    for (int j = 0; j < ny; j++){

      bottom_edge_x[j] = left_bnd_x[id_yz(nz-1,j)];
      top_edge_x[j] = b_x;
     
   }


    for (int i = 0; i < nx; i++){


       left_edge_x[i] = left_bnd_x[id_yz(nz-1,0)] + i*distance_x_left;
      right_edge_x[i] = left_bnd_x[id_yz(nz-1,ny-1)] + i*distance_x_right;


     
   }

    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x,  top_edge_x, left_edge_x, right_edge_x,  back_bnd_x);


    distance_x_left=(b_x-left_bnd_x[id_yz(0,0)])/(nx-1);
    distance_x_right=(b_x-left_bnd_x[id_yz(0,ny-1)])/(nx-1);


    for (int j = 0; j < ny; j++){

      top_edge_x[j] = b_x;
      bottom_edge_x[j] = left_bnd_x[id_yz(0,j)];
     
   }


    for (int i = 0; i < nx; i++){

      left_edge_x[i] = left_bnd_x[id_yz(0,0)] + i*distance_x_left;
      right_edge_x[i] = left_bnd_x[id_yz(0,ny-1)] + i*distance_x_right;

     
   }

    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, front_bnd_x);
      }

    
    for(int j = 0 ; j< ny; j++){
      for(int i = 0 ; i< nx; i++){
	
	front_bnd_z[id_xy(j,i)] = a_z;
	back_bnd_z[id_xy(j,i)] = b_z;
      }
    }
  
}





void CurvilinearTransformation::getInterpolatedFace_fromBottomAndTop( int n_block, int n_left_right, int n_top_bottom, int face,
					   double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					   double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					   double* face_x, double* face_y , double* face_z){

 
  double* top_edge_x= new double[n_top_bottom];
  double* top_edge_y= new double[n_top_bottom];
  double* top_edge_z= new double[n_top_bottom];

  double* bottom_edge_x= new double[n_top_bottom];
  double* bottom_edge_y= new double[n_top_bottom];
  double* bottom_edge_z= new double[n_top_bottom];

  double* left_edge_x= new double[n_left_right];
  double* left_edge_y= new double[n_left_right];
  double* left_edge_z= new double[n_left_right];

  double* right_edge_x= new double[n_left_right];
  double* right_edge_y= new double[n_left_right];
  double* right_edge_z= new double[n_left_right];

  // face 0:left 1:right 2:front 3:back
  bool left_right_face;
  bool first_face;
  switch(face){
  case(0): left_right_face=true; first_face=true ; break; //left
  case(1): left_right_face=true; first_face=false; break; //right
  case(2): left_right_face=false;first_face=true;  break; //bottom
  case(3): left_right_face=false;first_face=false; break; //top
  }


  int nx_on_top_bottom_face;
  int nz_on_top_bottom_face;
  if(left_right_face){
    nz_on_top_bottom_face=n_top_bottom;    
    nx_on_top_bottom_face=n_block;

  }else{
    nz_on_top_bottom_face=n_block;    
    nx_on_top_bottom_face=n_top_bottom;
  }
  
  kernels::idx2 id_xz(nz_on_top_bottom_face,nx_on_top_bottom_face);

  //generate top and bottom edge

  //  index = id_xz(nz-1-0,nx-1);
  // index = id_xz(0,nx-1);  
  // if(left_right & !first){
  //   std::cout <<  "top" << std::endl;
  //   std::cout <<  top_bnd_x[index] << std::endl;
  //   std::cout <<  top_bnd_y[index] << std::endl;;
  //   std::cout <<  top_bnd_z[index] << std::endl;;
  //   std::cout <<  "bottom" << std::endl;
  //   std::cout <<  bottom_bnd_x[index] << std::endl;
  //   std::cout <<  bottom_bnd_y[index] << std::endl;;
  //   std::cout <<  bottom_bnd_z[index] << std::endl;;     
  //   std::exit(-1);
  // } 


  // for(int i=0 ; i< nz; i++){
  //   for(int j=0 ; j< nx; j++){

  //     std::cout <<  "top" << std::endl;
  //     std::cout <<  top_bnd_x[id_xz(i,j)] << std::endl;
  //     std::cout <<  top_bnd_y[id_xz(i,j)] << std::endl;;
  //     std::cout <<  top_bnd_z[id_xz(i,j)] << std::endl;;
  //     std::cout <<  std::endl;
  //   }
  // }

  //   for(int i=0 ; i< nz; i++){
  //   for(int j=0 ; j< nx; j++){

  //     std::cout <<  "bottom" << std::endl;
  //     std::cout <<  bottom_bnd_x[id_xz(i,j)] << std::endl;
  //     std::cout <<  bottom_bnd_y[id_xz(i,j)] << std::endl;;
  //     std::cout <<  bottom_bnd_z[id_xz(i,j)] << std::endl;;
  //     std::cout <<  std::endl;
     
  //   }
  // } 
  
  //  std::cout << "top bottom edge" << std::endl;
  int index;
  for(int i=0; i < n_top_bottom ; i++){
    
    // left right front back
    if(left_right_face){
      if(first_face){ //left
	index = id_xz(i,0);
      }else{ //right
	index = id_xz(i,n_block-1);			
      }
    }else{//front or back
      if(first_face){ //front
	index=id_xz(0,i);
      }else{ //back
	//	index=id_xz(nz-1,nx-1-i);
	index=id_xz(n_block-1,i);
	// std::cout<< std::endl;	
	// std::cout << i << std::endl;
	// std::cout << index << std::endl;
      }
    }

    //    index = left_right ? (first ? id_xz(i,0) : id_xz(nz-1-i,nx-1)) : (first ? id_xz(0,i) :  id_xz(nz-1,nx-1-i));

    
    top_edge_x[i]=top_bnd_x[index];
    top_edge_y[i]=top_bnd_y[index];
    top_edge_z[i]=top_bnd_z[index];

    bottom_edge_x[i]=bottom_bnd_x[index];
    bottom_edge_y[i]=bottom_bnd_y[index];
    bottom_edge_z[i]=bottom_bnd_z[index];


    // std::cout << top_edge_x[i] <<std::endl;
    // std::cout << top_edge_y[i] <<std::endl;
    // std::cout << top_edge_z[i] <<std::endl;
                    
    // std::cout << bottom_edge_x[i] <<std::endl;
    // std::cout << bottom_edge_y[i] <<std::endl;
    // std::cout << bottom_edge_z[i] <<std::endl;
  }


  //
  double h_x=(bottom_edge_x[0]-top_edge_x[0])/(n_left_right-1);
  double h_y=(bottom_edge_y[0]-top_edge_y[0])/(n_left_right-1);
  double h_z=(bottom_edge_z[0]-top_edge_z[0])/(n_left_right-1);  

  

  // std::cout << "front edge" << std::endl;
  // std::cout << "hx: " << h_x<< "hy: " << h_y << "hz: " << h_z<<std::endl;
  for(int i=0; i < n_left_right ; i++){
    left_edge_x[i]=top_edge_x[0]+i*h_x;
    left_edge_y[i]=top_edge_y[0]+i*h_y;
    left_edge_z[i]=top_edge_z[0]+i*h_z;

    // std::cout << front_edge_x[i] <<std::endl;
    // std::cout << front_edge_y[i] <<std::endl;
    // std::cout << front_edge_z[i] <<std::endl;
  }

  h_x=(bottom_edge_x[n_top_bottom-1]-top_edge_x[n_top_bottom-1])/(n_left_right-1);
  h_y=(bottom_edge_y[n_top_bottom-1]-top_edge_y[n_top_bottom-1])/(n_left_right-1);
  h_z=(bottom_edge_z[n_top_bottom-1]-top_edge_z[n_top_bottom-1])/(n_left_right-1);  

  //  std::cout << "back edge" << std::endl;
  //  std::cout << "hx: " << h_x<< "hy: " << h_y << "hz: " << h_z<<std::endl;
  for(int i=0; i < n_left_right ; i++){
    right_edge_x[i]=top_edge_x[n_top_bottom-1]+i*h_x;
    right_edge_y[i]=top_edge_y[n_top_bottom-1]+i*h_y;
    right_edge_z[i]=top_edge_z[n_top_bottom-1]+i*h_z;

    // std::cout << back_edge_x[i] <<std::endl;
    // std::cout << back_edge_y[i] <<std::endl;
    // std::cout << back_edge_z[i] <<std::endl;

  }

  transFiniteInterpolation_singleCoordinate( n_top_bottom,  n_left_right,  left_edge_x,  right_edge_x,  top_edge_x,  bottom_edge_x,  face_x );
  transFiniteInterpolation_singleCoordinate( n_top_bottom,  n_left_right,  left_edge_y,  right_edge_y,  top_edge_y,  bottom_edge_y,  face_y );
  transFiniteInterpolation_singleCoordinate( n_top_bottom,  n_left_right,  left_edge_z,  right_edge_z,  top_edge_z,  bottom_edge_z,  face_z );


   // if(!left_right & !first){
   //  std::cout <<  "front" << std::endl;
   //  for(int i=0; i < ny ; i++){
   //  std::cout <<  front_edge_x[i] << std::endl;
   //  std::cout <<  front_edge_y[i] << std::endl;;
   //  std::cout <<  front_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;          
   //  }
   //  std::cout <<  "back" << std::endl;
   //  for(int i=0; i < ny ; i++){
   //    std::cout <<  back_edge_x[i] << std::endl;
   //    std::cout <<  back_edge_y[i] << std::endl;;
   //    std::cout <<  back_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;            
   //  }

   //  std::cout <<  "top" << std::endl;
   //  for(int i=0; i < nz ; i++){
   //    std::cout <<  top_edge_x[i] << std::endl;
   //    std::cout <<  top_edge_y[i] << std::endl;;
   //    std::cout <<  top_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;            
   //  }

   //  std::cout <<  "bottom" << std::endl;
   //  for(int i=0; i < nz ; i++){
   //    std::cout <<  bottom_edge_x[i] << std::endl;
   //    std::cout <<  bottom_edge_y[i] << std::endl;;
   //    std::cout <<  bottom_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;      
   //  }

    

    //    std::exit(-1);
  //} 

  

  // kernels::idx2 id_zy(ny,nz);

  // std::cout << "face" <<std::endl;
  // std::cout << face <<std::endl;  


  // std::cout << "x" <<std::endl;    
  // for(int i=0; i < nz ; i++){
  //   for(int j=0; j < ny ; j++){
  //     std::cout << face_x[id_zy(j,i)] <<std::endl;
  //   }
  // }
   
  // std::cout << "y" <<std::endl;    
  //   for(int i=0; i < nz ; i++){
  //   for(int j=0; j < ny ; j++){
  //     std::cout << face_y[id_zy(j,i)] <<std::endl;
  //   }
  // }
    
  // std::cout << "z" <<std::endl;    
  //     for(int i=0; i < nz ; i++){
  //   for(int j=0; j < ny ; j++){
  //     std::cout << face_z[id_zy(j,i)] <<std::endl;
  //   }
  // }

}  


void CurvilinearTransformation::getBoundaryCurves(int num_points,double offset_x, double offset_y,double width_x, double width_y ,double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y){

  double pi = 3.14159265359;

  int nx = 1/width_x*num_points;
  int ny = 1/width_y*num_points;

  double dx= 1.0/(nx-1);
  double dy= 1.0/(ny-1);

  // int i_0 =  (offset_x-0.0)/width_x;
  // int j_0 =  (offset_y-0.0)/width_y;
  

  //std::cout<<i_0<<std::endl;
  //std::cout<<j_0<<std::endl;
  //std::cout<<std::endl;

  
  // for(int i = 0 ; i< num_points; i++){
  //   left_bnd_x[i] =  offset_x;
  //   right_bnd_x[i] = width_x+offset_x;
  //   bottom_bnd_x[i] = width_x*dx*i + offset_x;
  //   top_bnd_x[i] = width_x*dx*i + offset_x;

  //   left_bnd_y[i] = width_y*dy*i + offset_y;
  //   right_bnd_y[i] = width_y*dy*i + offset_y;
  //   bottom_bnd_y[i] =offset_y;
  //   top_bnd_y[i] = width_y+offset_y + 0.0*std::sin(2*pi*top_bnd_x[i]);
  // }


  for(int i = 0 ; i< nx; i++){
    
    bottom_bnd_x[i] = 0 + dx*i;
    top_bnd_x[i]    = 0 + dx*i;
    
    bottom_bnd_y[i] = 0.0;
    top_bnd_y[i]    = 1 + 0.0*std::sin(2*pi*top_bnd_x[i]);

    //std::cout<<top_bnd_x[i]<<std::endl;
    //std::cout<<i<<std::endl;
    //std::cout<<std::endl;
  }

  //std::exit(-1);
  
  for(int i = 0 ; i< ny; i++){
    
    left_bnd_x[i]  = 0;
    right_bnd_x[i] = 1;

    left_bnd_y[i]  = dy*i;
    right_bnd_y[i] = dy*i;
  }


  double h0y = (top_bnd_y[0] - bottom_bnd_y[0])/(ny-1);
  double hny = (top_bnd_y[nx-1] - bottom_bnd_y[nx-1])/(ny-1);

  for(int i = 0 ; i< ny; i++){

    left_bnd_y[i]  = bottom_bnd_y[0] + h0y*i;
    right_bnd_y[i] = bottom_bnd_y[nx-1] + hny*i;
    
  }
  
}


void CurvilinearTransformation::transFiniteInterpolation3D(int mx, int my, int mz,
				int k_m, int k_p ,
				int j_m, int j_p ,
				int i_m, int i_p ,
				int num_points,
				double width_x, double width_y, double width_z,
				double* left_bnd,
				double* right_bnd,
				double* bottom_bnd,
				double* top_bnd,
				double* front_bnd,
				double* back_bnd,
				double* curvilinear
				){
  kernels::idx3 id_xyz(num_points,num_points,num_points);
  kernels::idx2 id_xy(my,mx);// back front
  kernels::idx2 id_xz(mz,mx);// bottom top
  kernels::idx2 id_yz(mz,my);//left right

   // double mesh_size_x = width_x/(num_points-1);
   // double mesh_size_y = width_y/(num_points-1);
   // double mesh_size_z = width_z/(num_points-1);


  //double mesh_size_x = 1.0/(num_points-1);
  //double mesh_size_y = 1.0/(num_points-1);
  //double mesh_size_z = 1.0/(num_points-1);

   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);
   double mesh_size_z = 1.0/(mz-1);

  // std::cout << mesh_size_x << " "<< mesh_size_y << " "<< mesh_size_z <<std::endl;
  // std::cout << mx << " "<< my << " "<< mz <<std::endl;

  int i_0;
  int j_0;
  int k_0;    

  double u,v,w;
  double uv,vw,uw,uvw1,uvw2;

  double q,r,s;

  for(int i=i_m ; i< i_p ; i++){
    i_0=i-i_m;
    q=mesh_size_x * i; //(i-std::floor((i+1)/num_points));
    //q=mesh_size_x *i_0;
    for(int j=j_m ; j< j_p ; j++){
      j_0=j-j_m;
      //r=mesh_size_y * (j-std::floor((j+1)/num_points));
      r=mesh_size_y *j;
      for(int k=k_m ; k< k_p ; k++){
	k_0=k-k_m;
	//s=mesh_size_z * (k-std::floor((k+1)/num_points));
	s=mesh_size_z * k;

	//std ::cout << i_0 << " "<< j_0 << " "<< k_0 <<std::endl;
	
	u=(1-q)*left_bnd[id_yz(k,j)]+q*right_bnd[id_yz(k,j)];

	v=(1-r)*top_bnd[id_xz(k,i)]+r*bottom_bnd[id_xz(k,i)];

	w=(1-s)*front_bnd[id_xy(j,i)]+s*back_bnd[id_xy(j,i)];

	uw=(1-q)* ((1 -r)*left_bnd[id_yz(k,0)] 
		   +r *left_bnd[id_yz(k,my-1)])
	  +q*((1-r)*right_bnd[id_yz(k,0)] 
	      +r *right_bnd[id_yz(k,my-1)]);

	uv=(1-r)*((1-s)*top_bnd[id_xz(0,i)]
		  + s*top_bnd[id_xz(mz-1,i)])
	  + r*((1-s)*bottom_bnd[id_xz(0,i)]
	       + s*bottom_bnd[id_xz(mz-1,i)]);

	vw = (1-s)*((1-q)*front_bnd[id_xy(j,0)] 
		    + q*front_bnd[id_xy(j,mx-1)])
	       + s*((1-q)*back_bnd[id_xy(j,0)] 
	              + q*back_bnd[id_xy(j,mx-1)]);
	
	uvw1=(1-q)*((1-r)*((1-s)*top_bnd[id_xz(0,0)] 
		              +s*top_bnd[id_xz(mz-1,0)])
		       +r*((1-s)*bottom_bnd[id_xz(0,0)]
			   +s*bottom_bnd[id_xz(mz-1,0)]));

	uvw2=q*((1-r)*((1-s)*top_bnd[id_xz(0,mx-1)] 
			 + s*top_bnd[id_xz(mz-1,mx-1)])
		+ r*((1-s)*bottom_bnd[id_xz(0,mx-1)]
		       + s*bottom_bnd[id_xz(mz-1,mx-1)]));

	curvilinear[id_xyz(k_0,j_0,i_0)] = u + v + w - uv - uw - vw + uvw1 + uvw2;


	// 	uw=(1-q)* ((1 -r)*left_bnd[id_yz(k,j_m)] 
	// 	   +r *left_bnd[id_yz(k,j_p-1)])
	//   +q*((1-r)*right_bnd[id_yz(k,j_m)] 
	//       +r *right_bnd[id_yz(k,j_p-1)]);

	// uv=(1-r)*((1-s)*top_bnd[id_xz(k_m,i)]
	// 	  + s*top_bnd[id_xz(k_p-1,i)])
	//   + r*((1-s)*bottom_bnd[id_xz(k_m,i)]
	//        + s*bottom_bnd[id_xz(k_p-1,i)]);

	// vw = (1-s)*((1-q)*front_bnd[id_xy(j,i_m)] 
	// 	    + q*front_bnd[id_xy(j,i_p-1)])
	//        + s*((1-q)*back_bnd[id_xy(j,i_m)] 
	//               + q*back_bnd[id_xy(j,i_p-1)]);
	
	// uvw1=(1-q)*((1-r)*((1-s)*top_bnd[id_xz(k_m,i_m)] 
	// 	              +s*top_bnd[id_xz(k_p-1,i_m)])
	// 	       +r*((1-s)*bottom_bnd[id_xz(k_m,i_m)]
	// 		   +s*bottom_bnd[id_xz(k_p-1,i_m)]));

	// uvw2=q*((1-r)*((1-s)*top_bnd[id_xz(k_m,i_p-1)] 
	// 		 + s*top_bnd[id_xz(k_p-1,i_p-1)])
	// 	+ r*((1-s)*bottom_bnd[id_xz(k_m,i_p-1)]
	// 	       + s*bottom_bnd[id_xz(k_p-1,i_p-1)]));

	// curvilinear[id_xyz(k_0,j_0,i_0)] = u + v + w - uv - uw - vw + uvw1 + uvw2;

      }
    }
  }
}


void CurvilinearTransformation::transFiniteInterpolation_singleCoordinate(int mx, int my, double* left_bnd, double* right_bnd, double* bottom_bnd, double* top_bnd, double* curvilinear ){

   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);

   
   // double mesh_size_x = width_x/(num_points-1);
   // double mesh_size_y = width_y/(num_points-1);
   // double mesh_size_z = width_z/(num_points-1);
  

   double r;
   double q;

   
   // std::cout << mesh_size_x << std::endl;
   // std::cout << mesh_size_y << std::endl;

   // std::cout << my << std::endl;
   // std::cout << mx << std::endl;

   
   //   std::exit(-1);

   kernels::idx2 id_xy(my,mx);
   
   for(int j =0 ; j < my ; j++) {
     for(int i = 0 ; i < mx ; i++) {
       
       q = (i)*mesh_size_x;
       r = (j)*mesh_size_y;
       
       curvilinear[id_xy(j,i)] = (1-q)*left_bnd[j]+q*right_bnd[j]+(1-r)*bottom_bnd[i]+r*top_bnd[i]-
	 (1-q)*(1-r)*left_bnd[0]-q*(1-r)*right_bnd[0]-r*(1-q)*top_bnd[0]-
	 (r*q)*top_bnd[mx-1];
     }
   }
}


void CurvilinearTransformation::transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y ){
      

   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);

   double r;
   double q;
   kernels::idx2 id_xy(num_points,num_points);
   
   int j_0 = 0;
  
   for(int j =j_m ; j < j_p ; j++) {
     int i_0 = 0;
     for(int i =i_m ; i < i_p ; i++) {
       q = (i)*mesh_size_x;
       r = (j)*mesh_size_y;
       
       curvilinear_x[id_xy(j_0,i_0)] = (1-q)*left_bnd_x[j]+q*right_bnd_x[j]+(1-r)*bottom_bnd_x[i]+r*top_bnd_x[i]-
	 (1-q)*(1-r)*left_bnd_x[0]-q*(1-r)*right_bnd_x[0]-r*(1-q)*top_bnd_x[0]-
	 (r*q)*top_bnd_x[mx-1];
       
       curvilinear_y[id_xy(j_0,i_0)] = (1-q)*left_bnd_y[j]+q*right_bnd_y[j]+(1-r)*bottom_bnd_y[i]+r*top_bnd_y[i]-
	 (1-q)*(1-r)*left_bnd_y[0]-q*(1-r)*right_bnd_y[0]-r*(1-q)*top_bnd_y[0]-
	 (r*q)*top_bnd_y[mx-1];
       
       i_0 = i_0+1; 
     }
     j_0 = j_0+1;
   }
}

// void CurvilinearTransformation::transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y ){
  
//   transFiniteInterpolation_singleCoordinate(mx, my, left_bnd_x,right_bnd_x,bottom_bnd_x,top_bnd_x,curvilinear_x );
//   transFiniteInterpolation_singleCoordinate(mx, my, left_bnd_y,right_bnd_y,bottom_bnd_y,top_bnd_y,curvilinear_y );
// }


double CurvilinearTransformation::lagrangeBasis(double x,int i,int num_points){
  double result=1;

  for (int j = 0 ; j< i ; j ++){
    // unifMesh
    result *= (x-unif_mesh[j]);
  }
  
  for (int j = i+1 ; j < num_points ; j ++){
    // 
    result *= (x-unif_mesh[j]);
  }

  result*= denominator_lagrange[i];
  
  return result;
}


void CurvilinearTransformation::interpolate3D(int x, int y ,int z, double* dest_mesh, int num_nodes,double& result){

  double a_x=0;
  double a_y=0;
  double a_z=0;  
  
  result=0;
  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);
  kernels::idx2 id_xy(num_nodes,num_nodes); //nodes,polynome

  for (int k = 0 ; k< num_nodes ; k++){    
    for (int j = 0 ; j< num_nodes ; j ++){
      for (int i = 0 ; i< num_nodes ; i ++){
	a_x=lagrange_basis_at_nodes[id_xy(x,i)];
	a_y=lagrange_basis_at_nodes[id_xy(y,j)];
	a_z=lagrange_basis_at_nodes[id_xy(z,k)];
	
	result += dest_mesh[id_xyz(k,j,i)] * a_x*a_y*a_z;
      }
    }
  }
}



void CurvilinearTransformation::interpolate(double x, double y, double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes,double& result){

  double a_x=0;
  double a_y=0;
  
  result=0;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);

  
  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
      a_x=lagrangeBasis(x,i,num_nodes);
      a_y=lagrangeBasis(y,j,num_nodes);
      result += dest_mesh[id_xy(j,i)] * a_x*a_y;
    }
  }
}



void CurvilinearTransformation::getValuesAtQuadNodes3D(double* dest_mesh, int num_nodes, double* results){
  // unifMesh unifMesh unifMesh

  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);


  for (int k = 0 ; k< num_nodes ; k ++){
    for (int j = 0 ; j< num_nodes ; j ++){
      for (int i = 0 ; i< num_nodes ; i ++){
	interpolate3D(i,j,k,dest_mesh,num_nodes,results[id_xyz(k,j,i)]);
      }
    }
  }
}



void CurvilinearTransformation::getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results){

  kernels::idx2 id_xy(num_nodes,num_nodes);

  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
#if defined(_GLL)
      interpolate(kernels::lobatto::nodes[num_nodes-1][num_nodes-1-i],kernels::lobatto::nodes[num_nodes-1][num_nodes-1-j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
      #else
      interpolate(kernels::legendre::nodes[num_nodes-1][i],kernels::legendre::nodes[num_nodes-1][j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
#endif
    }
  }
}

void CurvilinearTransformation::computeDerivatives_x(int i, int j , double* values , int num_nodes, double& der_x, double dx){
   kernels::idx2 id_xy(num_nodes,num_nodes);

   der_x = 0.0;

   for (int n = 0 ; n< num_nodes ; n ++){
     der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xy(j,n)]/dx;
   }
}

void CurvilinearTransformation::computeDerivatives_y (int i, int j , double* values , int num_nodes, double& der_y, double dy){
  kernels::idx2 id_xy(num_nodes,num_nodes); 

  der_y = 0.0;
  for (int n = 0 ; n< num_nodes ; n ++){
    der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xy(n,i)]/dy;
  }
}


void CurvilinearTransformation::computeDerivatives_x_3D(int i, int j , int k, double* values , int num_nodes, double& der_x, double dx){
  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);

  der_x = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xyz(k,j,n)]/dx;
  }

}

void CurvilinearTransformation::computeDerivatives_y_3D (int i, int j,int k, double* values , int num_nodes, double& der_y, double dy){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes); 

  der_y = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xyz(k,n,i)]/dy;
  }
}

void CurvilinearTransformation::computeDerivatives_z_3D (int i, int j , int k ,double* values , int num_nodes, double& der_z, double dz){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes); 

  der_z = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_z += kernels::dudx[num_nodes-1][k][n] * values[id_xyz(n,j,i)]/dz;
  }

}

void CurvilinearTransformation::metricDerivativesAndJacobian(int num_nodes, double* curvilinear_x, double* curvilinear_y,double* gl_vals_x,double* gl_vals_y,double* q_x,double* q_y,double* r_x,double* r_y,double* jacobian, double dx, double dy){


  getValuesAtQuadNodes(unif_mesh,
		       unif_mesh,
		       curvilinear_x,
		       num_nodes,
		       gl_vals_x);

  getValuesAtQuadNodes(unif_mesh,
		       unif_mesh,
		       curvilinear_y,
		       num_nodes,
		       gl_vals_y);

  double x_der_x; 
  double x_der_y;
  double y_der_x; 
  double y_der_y;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);
  
  for(int j = 0 ; j < num_nodes ; j ++){
    for(int i = 0 ; i< num_nodes ; i ++){

      computeDerivatives_x(j,i,gl_vals_x,num_nodes,x_der_x, dx);
      computeDerivatives_y(j,i,gl_vals_x,num_nodes,x_der_y, dy);

      computeDerivatives_x(j,i,gl_vals_y,num_nodes,y_der_x, dx);
      computeDerivatives_y(j,i,gl_vals_y,num_nodes,y_der_y, dy);

// J = x_q.*y_r - y_q.*x_r; % Jacobian (determinant of metric)
// q_x = y_r./J;   % 1/J y_r
// r_x = - y_q./J; % -1/J y_q
// q_y = - x_r./J; % -1/J x_r
// r_y = x_q./J;   %  1/J x_q
		     
     jacobian[id_xy(j,i)]=x_der_x*y_der_y-x_der_y*y_der_x;
     
     q_x[id_xy(j,i)]=y_der_y/jacobian[id_xy(j,i)];
     q_y[id_xy(j,i)]=-x_der_y/jacobian[id_xy(j,i)];
     r_x[id_xy(j,i)]=-y_der_x/jacobian[id_xy(j,i)];
     r_y[id_xy(j,i)]=x_der_x/jacobian[id_xy(j,i)];

     //      printf("%f\n",der_x);
     //      printf("%f\n",der_y);
    }
 }

  


}


void CurvilinearTransformation::metricDerivativesAndJacobian3D(int num_nodes,
				  double* curvilinear_x, double* curvilinear_y, double* curvilinear_z,
				  double* gl_vals_x, double* gl_vals_y, double* gl_vals_z,
				  double* q_x, double* q_y, double* q_z,
				  double* r_x, double* r_y, double* r_z,
				  double* s_x, double* s_y, double* s_z,				  
				  double* jacobian,
				  double dx, double dy, double dz
				  ){

  getValuesAtQuadNodes3D(curvilinear_x,num_nodes,gl_vals_x);
  getValuesAtQuadNodes3D(curvilinear_y,num_nodes,gl_vals_y);
  getValuesAtQuadNodes3D(curvilinear_z,num_nodes,gl_vals_z);  

  double x_der_x; 
  double x_der_y;
  double x_der_z;  
  
  double y_der_x; 
  double y_der_y;
  double y_der_z;

  double z_der_x; 
  double z_der_y;
  double z_der_z;  

  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);
  
  for(int k = 0 ; k < num_nodes ; k ++){
    for(int j = 0 ; j < num_nodes ; j ++){
      for(int i = 0 ; i< num_nodes ; i ++){

  	computeDerivatives_x_3D(i,j,k,gl_vals_x,num_nodes,x_der_x, dx);
  	computeDerivatives_y_3D(i,j,k,gl_vals_x,num_nodes,x_der_y, dy);
  	computeDerivatives_z_3D(i,j,k,gl_vals_x,num_nodes,x_der_z, dz);

	computeDerivatives_x_3D(i,j,k,gl_vals_y,num_nodes,y_der_x, dx);
	computeDerivatives_y_3D(i,j,k,gl_vals_y,num_nodes,y_der_y, dy);
	computeDerivatives_z_3D(i,j,k,gl_vals_y,num_nodes,y_der_z, dz);

	computeDerivatives_x_3D(i,j,k,gl_vals_z,num_nodes,z_der_x, dx);
	computeDerivatives_y_3D(i,j,k,gl_vals_z,num_nodes,z_der_y, dy);
	computeDerivatives_z_3D(i,j,k,gl_vals_z,num_nodes,z_der_z, dz);
	
	jacobian[id_xyz(k,j,i)]=x_der_x*(y_der_y*z_der_z-y_der_z*z_der_y)
	  -x_der_y*(y_der_x*z_der_z-y_der_z*z_der_x)
	  +x_der_z*(y_der_x*z_der_y-y_der_y*z_der_x);
	
	q_x[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(y_der_y*z_der_z - z_der_y*y_der_z);
	r_x[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(z_der_x*y_der_z - y_der_x*z_der_z);
	s_x[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(y_der_x*z_der_y - z_der_x*y_der_y);
	
	q_y[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(z_der_y*x_der_z - x_der_y*z_der_z);
	r_y[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(x_der_x*z_der_z - z_der_x*x_der_z);
	s_y[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(z_der_x*x_der_y - x_der_x*z_der_y);
	
	q_z[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(x_der_y*y_der_z - y_der_y*x_der_z);
	r_z[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(y_der_x*x_der_z - x_der_x*y_der_z);
	s_z[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(x_der_x*y_der_y - y_der_x*x_der_y);

      }
    }
  }

  


}


double CurvilinearTransformation::interpolate_fault_surface(double* X1, double* Y1, double* Z1, double y, double z, int my, int mz)
{
//real(kind = wp), dimension(:,:),intent(in) :: X1, Y1, Z1       ! data
//    real(kind = wp),intent(in) :: y, z
//    real(kind = wp),intent(inout) :: gg
//    integer,intent(in) :: my, mz
  
  double dz, dy, r, r0;
  
  int jj, kk, nmy, npy, nmz, npz, stenc, ndp;
  
  kernels::idx2 id_yz(mz,my);
  
  dy = Y1[id_yz(0,1)] - Y1[id_yz(0,0)];
  dz = Z1[id_yz(1,0)] - Z1[id_yz(0,0)];
  
  r0 = std::sqrt(pow(dy,2) + pow(dz, 2));
    
  stenc = 2;
  ndp = 5;

  double* yy = new double [my];
  double* zz = new double [mz];

  double gg;

 
                
  for (int k = 0; k<mz; k++){
    for (int j = 0; j<my; j++){

      r = std::sqrt(pow(y- Y1[id_yz(k,j)], 2) + pow(z-Z1[id_yz(k,j)], 2));

      if (r < 2*r0) {
             
	if (std::sqrt(pow(y-Y1[id_yz(k,j)], 2)) <= dy) {
                
	  nmy = j-stenc;
	  npy = j+stenc;
            
	  if (j == 0) {
                   
	    nmy = j;
	    npy = j+ndp-1;
                   
	  }
	  
	  if (j == 1) {
                   
	    nmy = j-1;
	    npy = j+ndp-2;
            
	  }

	  if (j == 2) {
                   
	    nmy = j-2;
	    npy = j+ndp-3;
            
	  }

	   if (j == my-3){
                   
	    nmy = j-(ndp-3);
	    npy = j+2;
                   
	  }
                
	  if (j == my-2){
                   
	    nmy = j-(ndp-2);
	    npy = j+1;
                   
	  }

	  if (j == my-1){
                   
	    nmy = j-(ndp-1);
	    npy = j;
            
	  }
                
                
	  if  (std::sqrt(pow(z-Z1[id_yz(k,j)], 2)) <= dz) {
	    
	    nmz = k-stenc;
	    npz = k+stenc;
                   
	    if (k == 0){
	      
	      nmz = k;
	      npz = k+ndp-1;
	    }        
       

	    if (k == 1) {
                      
	      nmz = k-1;
	      npz = k+ndp-2;
              
	    }

	     if (k == 2) {
                      
	      nmz = k-2;
	      npz = k+ndp-3;
              
	    }
                   
	    if (k == mz-1){
                      
	      nmz = k-(ndp-1);
	      npz = k;
                      
	    }

	    if (k == mz-2) {
                      
	      nmz = k-(ndp-2);
	      npz = k+1;
                      
	    }

	    if (k == mz-3) {
                      
	      nmz = k-(ndp-3);
	      npz = k+2;
                      
	    }

	    //std::cout<< "y " << std::endl;
	    for (int j0 = 0; j0<my; j0++){
	      yy[j0] = Y1[id_yz(k,j0)];
	      //std::cout<<j0<<" "<<yy[j0] << std::endl;
	    }

	    
	    //std::cout<< "z " << std::endl;
	    for (int k0 = 0; k0<mz; k0++){
	      zz[k0] = Z1[id_yz(k0,j)];
	      
	      //std::cout<<k0<<" "<<zz[k0] << std::endl;
	    }

	    //std::exit(-1);
	    
	     gg = interpol2d_dG(my, mz, nmy, npy, nmz, npz, y, z, yy, zz, X1);
                                    
                
	  }
	  //==============================================
	}
      }
    }
  }

  return gg;
  
}

double  CurvilinearTransformation::interpol2d_dG(int my, int mz, int nmy, int npy, int nmz, int npz, double y, double z, double* yj, double* zk, double* f)
{
  //!contruct the langrange interpolation g(x,y) of  f(x,y)
  
  //! g(x,y): lagrange interpolant
  //! nmx, nmy: lower bound
  //! npx, npy: upper bound
  //! xi,yj: data points
  //! x,y : interpolation point
  //  integer, intent(in) :: nmx, nmy, npx, npy             ! number of rows and column
  //  real(kind = wp), intent(out) :: g                     ! lagrange interpolant 
  //  real(kind = wp), intent(in) :: x, y                   ! interpolation point
  //  real(kind = wp), dimension(:),intent(in) :: xi, yj    ! data points
  //  real(kind = wp), dimension(:,:),intent(in) :: f       ! data
  
  double a_j, a_k, g;

 
    
  
  g = 0.0;
  
  kernels::idx2 id_yz(mz,my);
    
  for (int k = nmz; k< npz+1; k++){
    for (int j = nmy; j <npy+1; j++){
      
      a_j = lagrange(nmy, npy, j, y, yj);
      a_k = lagrange(nmz, npz, k, z, zk);
      g = g + a_j*a_k*f[id_yz(k,j)];
      
    }
  }

  return g;
  
}


int CurvilinearTransformation::getBlock(const tarch::la::Vector<DIMENSIONS,double>& center,
					const tarch::la::Vector<DIMENSIONS,double>& dx){

  double recWidth_x = b_x-a_x; //width in x of the rectangle

  //relative position of the fault on the unit square
  double fault_ref = (fault_position-a_x)/(recWidth_x);

  //defines which block we are in
  int n = center[0] > fault_ref ? 1 : 0 ;

  return n;
}

double CurvilinearTransformation::lagrange(int m, int p, int i, double x, double *xi){
  
  //!Function to calculate  Lagrange polynomial for order N and polynomial
  //![nm, np] at location x.
  
  //! nm: lower bound
  //! np: upper bound
  //! xi: data points
  //! x : interpolation point
  
  //integer, intent(in) :: m, p, i
  //real(kind = wp), intent(in) :: x, xi(:)
  
  double h, num, den;
  
  h = 1.0;
  
  for (int j = m; j <p+1; j++){
    if (j != i) {
      
      num = x - xi[j];
      den = xi[i] - xi[j];
      h = h*num/den;
      
      //std::cout<< i << " " << j << " "<< xi[i]-xi[j]<< std::endl;
      //std::cout<< x << " " << den << " "<< h << std::endl;
    }
  }

 
  return h;
}

  

