#include "kernels/KernelUtils.h"
#include "kernels/DGMatrices.h"
#include "CurvilinearTransformation.h"
#if defined(_GLL)
#include "kernels/GaussLobattoQuadrature.h"
#else
#include "kernels/GaussLegendreQuadrature.h"
#endif

#if defined(USE_ASAGI)
#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "reader/asagi_reader.h"
#endif

#define PI 3.14159265359;

CurvilinearTransformation::CurvilinearTransformation(const int num_nodes, const double dx ,const double fault_position,
						     const double* const left_vertex,
						     const double* const right_vetrex):
  _num_nodes(num_nodes),
  _dx(dx),
  _fault_position(fault_position)
{
  _left_vertex[0] =left_vertex[0];
  _left_vertex[1] =left_vertex[1];
  _left_vertex[2] =left_vertex[2];

  _right_vertex[0]=right_vertex[0];
  _right_vertex[1]=right_vertex[1];
  _right_vertex[2]=right_vertex[2];
  
  _rect_width[0]=_right_vertex[0]-_left_vertex[0];
  _rect_width[1]=_right_vertex[1]-_left_vertex[1];
  _rect_width[2]=_right_vertex[2]-_left_vertex[2];
  
  _ne_x=std::round(_rect_width[0],_dx); //number of elements in x direction on whole domain
  _ne_y=std::round(_rect_width[1],_dx); //number of elements in y direction on whole domain
  _ne_z=std::round(_rect_width[2],_dx); //number of elements in z direction on whole domain
  
  _nx=_ne_x*(_num_nodes-1)+1; //global number of nodes in x direction considering collocation
  _ny=_ne_y*(_num_nodes-1)+1; //global number of nodes in y direction considering collocation
  _nz=_ne_z*(_num_nodes-1)+1; //global number of nodes in z direction considering collocation

  double fault_ref=(fault_position-left_vertex[0])/(recWidth[0]);   //relative position of the fault
  
  for(int n=0 ; n<2 ; ++n){
    if(n == 0){
      _block_width_x[0]=fault_position -left_vertex[0];
      _ne_x_block[0] = std::round(_ne_x*fault_ref);
    }else{
      _block_width_x[1]=right_vertex[0]-fault_position;
      _ne_x_block[1] = _ne_x-std::round(_ne_x*fault_ref);
    }

    _ne_x_bock[n]=min(_ne_x-1,max(1,_ne_x_block[n])); // at least one element for each block
    _nx_block [n]=ne_x_block *(_num_nodes-1)+1;
    _boundary_x[n]=new Boundary_single_coordinate(nx_block,ny,nz);
    _boundary_y[n]=new Boundary_single_coordinate(nx_block,ny,nz);
    _boundary_z[n]=new Boundary_single_coordinate(nx_block,ny,nz);

    getBoundaryCurves3D_cutOffTopography_withFault(n,boundary_x,boundary_y,boundary_z);
  }

  //init uniform mesh 
  unif_mesh= new double[num_nodes];
  for(int i = 0; i< num_nodes ; i++){
    unif_mesh[i]=i*1.0/(num_nodes-1);
  }

  //precompute denominator for lagrange interpolation on uniform mesh for all vertices
  denominator_lagrange= new double[num_nodes];
  for(int i=0; i< num_nodes ; i++){
    double temp = 1.0;
    for (int j = 0 ; j< i ; j ++){
      temp = temp*(unif_mesh[i]-unif_mesh[j]);
    }
    for (int j = i+1 ; j < num_nodes ; j ++){
      temp = temp*(unif_mesh[i]-unif_mesh[j]);
    }
    denominator_lagrange[i] = 1.0/temp;
  }
  
  kernels::idx2 id_xy(num_nodes,num_nodes); //nodes,polynome
  lagrange_basis_at_nodes = new double[num_nodes*num_nodes];
  kernels::initGaussLegendreNodesAndWeights(std::set<int>()); //empty set as it is not used

  for(int i=0; i< num_nodes ; i++){
    for (int j = 0 ; j< num_nodes ; j ++){
#if defined(_GLL)
      lagrange_basis_at_nodes[id_xy(j,i)]
	=lagrangeBasis_uniform(kernels::gaussLobattoNodes[num_nodes-1][num_nodes-j],i,num_nodes);
#else
      lagrange_basis_at_nodes[id_xy(j,i)]
	=lagrangeBasis_uniform(kernels::gaussLegendreNodes[num_nodes-1][j],i,num_nodes);
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
  int n = getBlock(center,dx); //find current block

  //index of the first node within the block
  int i_m;
  //within the block = the domain
  int j_m = std::round((offset_y)/dx[1])*(_num_nodes-1);
  int k_m = std::round((offset_z)/dx[2])*(_num_nodes-1);    
  
  if(n == 0){
    i_m = std::round( offset_x                /width_x)*(_num_nodes-1);
  }else{
    i_m = std::floor((offset_x-fault_position)/width_x)*(_num_nodes-1); //subtract fault position from offset
  }
  // index of the last node within the block
  int i_p = i_m + num_nodes;
  int j_p = j_m + num_nodes;
  int k_p = k_m + num_nodes;

  double* curvilinear_x = new double[num_nodes*num_nodes*num_nodes];
  double* curvilinear_y = new double[num_nodes*num_nodes*num_nodes];
  double* curvilinear_z = new double[num_nodes*num_nodes*num_nodes];  
  //perform transfinite interpolation along boundary curves
  transFiniteInterpolation3D(k_m, k_p,
			     j_m, j_p,
			     i_m, i_p,
			     dx,
			     boundary_x,
			     curvilinear_x);
  
  transFiniteInterpolation3D(k_m, k_p,
			     j_m, j_p,
			     i_m, i_p,
			     dx,
			     boundary_y,
			     curvilinear_y);

  transFiniteInterpolation3D(k_m, k_p,
			     j_m, j_p,
			     i_m, i_p,
			     dx,
			     boundary_z,
			     curvilinear_z);

  metricDerivativesAndJacobian3D(curvilinear_x,
				 curvilinear_y,
				 curvilinear_z,
				 gl_vals_x, gl_vals_y, gl_vals_z,
				 q_x, q_y, q_z,
				 r_x, r_y, r_z,
				 s_x, s_y, s_z,         
				 jacobian,
				 dx);
}

double CurvilinearTransformation::fault(double y, double z){
  double pi = 3.14159265359;
  double angle = 60.0/360.0 * 2 * pi;
  double Ly = _right_vertex[1]-_left_vertex[1];
  double Lz = _right_vertex[2]-_left_vertex[2];
  double fault_surface;
  
  fault_surface=0.25*(std::sin(2*pi*y/Ly)*std::cos(2*pi*y/Ly))*std::sin(2*pi*z/Lz)*std::cos(2*pi*z/Lz);
  //  fault_surface=0;
  return  fault_surface;

}


double CurvilinearTransformation::topography(double x, double z, double depth){
  double pi = 3.14159265359;
  double angle = 60.0/360.0 * 2 * pi;
  double Lx = b_x-a_x;
  double Lz = b_z-a_z;
  double topo;

  topo = 1.0*(0.1*(x + z) + 0.25*depth*(std::sin(4*pi*x/Lx+3.34)*std::cos(4*pi*x/Lx)
				 * std::sin(4*pi*z/Lz+3.34)*std::cos(4*pi*z/Lz)));
  return topo;
}  

#ifdef(USE_ASAGI)
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

void CurvilinearTransformation::getBoundaryCurves3D_cutOffTopography_withFault(n, boundary_x, boundary_y, boundary_z){
  int nx=_nx_block[n];
  int ny=_ny;
  int nz=_nz;

  //get dx for block
  double dx= _block_width_x[n]/(nx-1); 
  double dy= _rect_width[1]   /(ny-1);
  double dz= _rect_width[2]   /(nz-1);

  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); // left right

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
  
  double x,y,z;

  //gen top surface and bottom face
  //z      fault_postion
  //^-----------------------o right vertex(x,z)
  //|           |           |
  //|           |           |
  //|           |           |
  //|           |           |
  //|     0     |     1     |
  //|           |           |
  //|           |           |
  //|           |           |
  //|           |           |
  //o-----------------------> x
  //left vertex(x,z)           top: left_vertex(y)+topography bottom: right_vertex(y)
  
  double top_left_vertice_x= (n==0) ? left_vertice[0] : fault_position;
  double depth= right_vertice[1];
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){
      top_bnd_x[id_xz(k,i)] = top_left_vertice_x+dx*i;
      top_bnd_y[id_xz(k,i)] = left_vertice[1];      
      top_bnd_z[id_xz(k,i)] = left_vertice[2]   +dz*k;

      x = top_bnd_x[id_xz(k,i)];
      z = top_bnd_z[id_xz(k,i)];
    }
  }

  for(int k = 0 ; k< nz; k++){
    for(int i = 0 ; i< nx; i++){
      bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      bottom_bnd_y[id_xz(k,i)] = right_vertice[1];      
      bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
    }
  }  

  //gen left right surface
  //y
  //^-----------------------o right vertex(y,z)
  //|                       |
  //|                       |
  //|                       |
  //|                       |
  //|                       |
  //|                       |
  //|                       |
  //|                       |
  //|                       |
  //o-----------------------> z
  //left vertex(y,z)           n==0 left=left vertex(x) right=fault_position
  //                           n==1 left=fault_position right=right_vertex(x)

  double left_boundary_left_vertice_x = n==0 ? left_vertice[0] : fault_position;
  for(int k = 0 ; k< nz; k++){
    for(int j = 0 ; j< ny; j++){
      left_bnd_x[id_yz(k,j)]=left_boundary_left_vertice_x;
      left_bnd_y[id_yz(k,j)]=left_vertice[1]+dy*j;
      left_bnd_z[id_yz(k,j)]=left_vertice[2]+dz*k;
    }
  }
  
  double right_boundary_right_vertice_x = n==0 ? fault_position : right_vertice[0];
  for(int k = 0 ; k< nz; k++){
    for(int j = 0 ; j< ny; j++){
      right_bnd_x[id_yz(k,j)]=right_boundary_right_vertice_x;
      right_bnd_y[id_yz(k,j)]=left_vertice[1]+dy*j;
      right_bnd_z[id_yz(k,j)]=left_vertice[2]+dz*k;
    }
  }


  //generate syntetic fault for interpolation
  //----------------------------------------//
  //LR TODO: Why do we need this ?)
  //IMO: This is equivalent to just evaluating fault at bnd_x and bnd_y)
  double* x_synt_fault = new double[ny*nz];
  double* y_synt_fault = new double[ny*nz];
  
  for(int k = 0 ; k< nz; k++){
    for(int j = 0 ; j< ny; j++){
      if(n==0){
	y = right_bnd_y[id_yz(k,j)];
	z = right_bnd_z[id_yz(k,j)];
      }else if(n==1){
	y = left_bnd_y[id_yz(k,j)];
	z = left_bnd_z[id_yz(k,j)];
      }
      y_synt_fault[id_yz(k,j)] = _left_vertice[1]-_rect_width[1]*0.3 + 1.3*dy*j; //stretch domain in y by 30% to cover intersection with topography
      x_synt_fault[id_yz(k,j)] = _fault_position - fault(y_synt_fault[id_yz(k,j)],left_bnd_z[id_yz(k,j)]);
    }
  }

  for(int k = 0 ; k< nz; k++){ 
    for(int j = 0 ; j< ny; j++){
      if (n == 0){
	y = right_bnd_y[id_yz(k,j)];
	z = right_bnd_z[id_yz(k,j)];
	interpolate_fault_surface(x_synt_fault, y_synt_fault, right_bnd_z, y, z, right_bnd_x[id_yz(k,j)]); //get interpolated x on syntetic fault
      }else if (n == 1){
	y = left_bnd_y[id_yz(k,j)];
	z = left_bnd_z[id_yz(k,j)];
	interpolate_fault_surface(x_synt_fault, y_synt_fault, left_bnd_z, y, z, left_bnd_x[id_yz(k,j)]); //get interpolated x on syntetic fault
      }
    }
  }

  delete[] x_synt_fault;
  delete[] y_synt_fault;
  //----------------------------------------//

  
  double* top_edge_x = new double[nz];
  double* bottom_edge_x = new double[nz];
  double* left_edge_x = new double[nx];
  double* right_edge_x = new double[nx];
  double distance_x_left;  
  double distance_x_right;
  if (n == 0) {

    /** bottom_bnd_x:
     *                         right_bnd_x(0,ny-1)	   
     *                   x-----x                  
     *          ^        |     )         
     *          |	 |      )        
     *    left_vertice_x |       (       right_bnd_x(:,ny-1)
     *          |	 |        )      
     *          v	 |         )      
     *          	 x---------x
     *           	        right_bnd_x(nz-1,ny-1)
     **/
    
    distance_x_left=(right_bnd_x[id_yz(0,ny-1)]-left_vertice[0])/(nx-1);
    distance_x_right=(right_bnd_x[id_yz(nz-1,ny-1)]-left_vertice[0])/(nx-1);
    
    for (int k = 0; k < nz; k++){
      top_edge_x[k] = right_bnd_x[id_yz(k,ny-1)];
      bottom_edge_x[k] = left_vertice[0];
    }
    for (int i = 0; i < nx; i++){
      left_edge_x[i] = a_x + i*distance_x_left;
      right_edge_x[i] = a_x + i*distance_x_right;
    }
      
    transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, bottom_bnd_x);

    /** top_bnd_x:
     *                         right_bnd_x(0,0)	   
     *                   x-----x                  
     *          ^        |     )         
     *          |	 |      )        
     *    left_vertice_x |       (       right_bnd_x(:,0)
     *          |	 |        )      
     *          v	 |         )      
     *          	 x---------x
     *           	        right_bnd_x(nz-1,0)
     **/    
      
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
  }else if (n == 1) {

    /** bottom_bnd_x:
     *
     *      left_bnd_x(0,ny-1)
     *                  x---------x         
     *                  )         |         ^
     *                   )        |         |
     * left_bnd_x(:,ny-1) (       |  right_vertice_x
     *                     )      |         |
     *                     )      |         v
     *                      x-----x         
     *       left_bnd_x(nz-1,ny-1)          
     **/
    
    distance_x_left=(right_vertice[0] - left_bnd_x[id_yz(0,ny-1)])/(nx-1);
    distance_x_right=(right_vertice[0] - left_bnd_x[id_yz(nz-1,ny-1)])/(nx-1);
    
    for (int k = 0; k < nz; k++){
      top_edge_x[k] = right_vertice[0];
      bottom_edge_x[k] = left_bnd_x[id_yz(k,ny-1)];
    }
    for (int i = 0; i < nx; i++){
      left_edge_x[i] = left_bnd_x[id_yz(0,ny-1)] + i*distance_x_left;
      right_edge_x[i] =  left_bnd_x[id_yz(nz-1,ny-1)] + i*distance_x_right;
    }
      
    transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, bottom_bnd_x);

  /** top_bnd_x:
   *
   *      left_bnd_x(0,0)        
   *                x---------x         
   *                )         |         ^
   *                 )        |         |
   * left_bnd_x(:,0)  (       |  right_vertice_x
   *                   )      |         |
   *                   )      |         v
   *                    x-----x         
   *       left_bnd_x(nz-1,0)          
   **/
      
    distance_x_left=(right_vertice[0] - left_bnd_x[id_yz(0,0)])/(nx-1);
    distance_x_right=(right_vertice[0] - left_bnd_x[id_yz(nz-1,0)])/(nx-1);
    for (int k = 0; k < nz; k++){
      top_edge_x[k] = right_vertice[0];
      bottom_edge_x[k] = left_bnd_x[id_yz(k,0)];   
    }
    for (int i = 0; i < nx; i++){
      left_edge_x[i] = left_bnd_x[id_yz(0,0)] + i*distance_x_left;
      right_edge_x[i] =  left_bnd_x[id_yz(nz-1,0)] + i*distance_x_right;
    }
    // left right bottom top
    transFiniteInterpolation_singleCoordinate(nx, nz, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, top_bnd_x);
  }

  //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){
      x = top_bnd_x[id_xz(k,i)];
      z = top_bnd_z[id_xz(k,i)];
#if !defined(USE_ASAGI)      
      top_bnd_y[id_xz(k,i)] -= topography(x, z, depth);
#else
      top_bnd_y[id_xz(k,i)] -= topography_fromASAGI(x,z,topography, adapter, model);
#endif      
    }
  }

  double* top_edge_y = new double[nz];
  double* bottom_edge_y = new double[nz];
  double* left_edge_y = new double[ny];
  double* right_edge_y = new double[ny];
  double distance_y_left;        
  double distance_y_right;

  /** left_bnd_y:
   *
   *                     top_bnd_y(:,0)
   *                           /\  /x top_bnd_y(nz-1,0)
   *        top_bnd_y(0,0) x/\/  \/ |         ^
   *              ^        |        |         |
   *              |        |        |         |
   *              v        |        |         v
   *     bottom_bnd_y(0,0) x--------x bottom_bnd_y(nz-1,0)
   *                    bottom_bnd_y(:,0)
   **/
  
  distance_y_left=(bottom_bnd_y[id_xz(0,0)]-top_bnd_y[id_xz(0,0)])/(ny-1);
  distance_y_right=(bottom_bnd_y[id_xz(nz-1,0)]-top_bnd_y[id_xz(nz-1,0)])/(ny-1);
  
  for (int k = 0; k < nz; k++){
    top_edge_y[k] = top_bnd_y[id_xz(k,0)];
    bottom_edge_y[k] = bottom_bnd_y[id_xz(k,0)];
  }
  for (int j = 0; j < ny; j++){
    left_edge_y[j] = top_bnd_y[id_xz(0,0)] + j*distance_y_left;
    right_edge_y[j] = top_bnd_y[id_xz(nz-1,0)] + j*distance_y_right;
  }
  // left right bottom top
  transFiniteInterpolation_singleCoordinate(ny, nz, top_edge_y, bottom_edge_y, left_edge_y, right_edge_y, left_bnd_y);
  
  /** right_bnd_y:
   *
   *                    top_bnd_y(:,nx-1)
   *                           /\  /x top_bnd_y(nz-1,nx-1)
   *     top_bnd_y(0,nx-1) x/\/  \/ |         ^
   *              ^        |        |         |
   *              |        |        |         |
   *              v        |        |         v
   *  bottom_bnd_y(0,nx-1) x--------x bottom_bnd_y(nz-1,nx-1)
   *                   bottom_bnd_y(:,nx-1)
   **/
  
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
  // left right bottom top
  transFiniteInterpolation_singleCoordinate(ny, nz, top_edge_y, bottom_edge_y, left_edge_y, right_edge_y, right_bnd_y);

  double* top_edge_y = new double[nx];
  double* bottom_edge_y = new double[nx];
  
  double* left_edge_y = new double[ny];
  double* right_edge_y = new double[ny];
  
  /** front_bnd_y:
   *
   *                    top_bnd_y(0,:)
   *                           /\  /x top_bnd_y(0,nx-1)
   *        top_bnd_y(0,0) x/\/  \/ |         ^
   *              ^        |        |         |
   *              |        |        |         |
   *              v        |        |         v
   *     bottom_bnd_y(0,0) x--------x bottom_bnd_y(0,nx-1)
   *                   bottom_bnd_y(0,:)
   **/

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
  // left right bottom top
  transFiniteInterpolation_singleCoordinate(nx, ny, left_edge_y, right_edge_y, top_edge_y, bottom_edge_y,front_bnd_y);

  /** back_bnd_y:
   *
   *                    top_bnd_y(nz-1,:)
   *                           /\  /x top_bnd_y(nz-1,nx-1)
   *     top_bnd_y(nz-1,0) x/\/  \/ |         ^
   *              ^        |        |         |
   *              |        |        |         |
   *              v        |        |         v
   *  bottom_bnd_y(nz-1,0) x--------x bottom_bnd_y(nz-1,nx-1)
   *                   bottom_bnd_y(nz-1,:)
   **/
  
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
  // left right bottom top
  transFiniteInterpolation_singleCoordinate(nx, ny, left_edge_y, right_edge_y,  top_edge_y, bottom_edge_y,back_bnd_y);
  
  //Generate back and front boundary
  double* top_edge_x = new double[ny];
  double* bottom_edge_x = new double[ny];
  double* left_edge_x = new double[nx];
  double* right_edge_x = new double[nx];

  double distance_x_left;
  double distance_x_right;
 
  if (n == 0){
    /** back_bnd_x:
     *
     *                      right_bnd_x(nz-1,ny-1)         
     *       ^	       x---x      
     *       |         |    )       
     *       |	       |     (    
     *  left_vertice_x |      ) right_bnd_x(nz-1,:)
     *       |	       |       (  
     *       |	       |        ) 
     *       v	       x--------x 
     *                       right_bnd_x(nz-1,0)         	    
     **/
    
    distance_x_left=(right_bnd_x[id_yz(nz-1,0)]-0)/(nx-1);
    distance_x_right=(right_bnd_x[id_yz(nz-1,ny-1)]-0)/(nx-1);

    for (int j = 0; j < ny; j++){
      bottom_edge_x[j] = left_vertex[0]; //bottom is constant in x
      top_edge_x[j] = right_bnd_x[id_yz(nz-1,j)]; 
    }

    for (int i = 0; i < nx; i++){
      left_edge_x[i] = a_x + i*distance_x_left;
      right_edge_x[i] = a_x + i*distance_x_right;
    }
    // left right bottom top
    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x,  top_edge_x, left_edge_x, right_edge_x,  back_bnd_x); //interpolate back_bnd_x

    /** front_bnd_x:
     *
     *                      right_bnd_x(0,ny-1)         
     *       ^	       x---x      
     *       |         |    )       
     *       |	       |     (    
     *  left_vertice_x |      ) right_bnd_x(0,:)
     *       |	       |       (  
     *       |	       |        ) 
     *       v	       x--------x 
     *                       right_bnd_x(0,0)         	    
     **/
    
    distance_x_left=(right_bnd_x[id_yz(0,0)]-left_vertice[0])/(nx-1);
    distance_x_right=(right_bnd_x[id_yz(0,ny-1)]-left_vertice[0])/(nx-1);


    for (int j = 0; j < ny; j++){
      top_edge_x[j] = right_bnd_x[id_yz(0,j)];
      bottom_edge_x[j] = a_x;
    }

    for (int i = 0; i < nx; i++){
      left_edge_x[i] = a_x + i*distance_x_left;
      right_edge_x[i] = a_x + i*distance_x_right;
    }   
    // left right bottom top
    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, front_bnd_x); //interpolate front_bnd_y
  }else if (n == 1){
    /** back_bnd_x:
     *
     *    left_bnd_x(nz-1,ny-1)         
     *                  x--------x       ^
     *                  )        |       |
     *                   (       |       |
     * left_bnd_x(nz-1,:) )      | right_vertice_x
     *                     (     |       |
     *                      )    |       |
     *                       x---x       v
     *       left_bnd_x(nz-1,0)         
     **/    
    
    distance_x_left=(b_x-left_bnd_x[id_yz(nz-1,0)])/(nx-1);
    distance_x_right=(b_x-left_bnd_x[id_yz(nz-1,ny-1)])/(nx-1);

    for (int j = 0; j < ny; j++){
      bottom_edge_x[j] = left_bnd_x[id_yz(nz-1,j)];
      top_edge_x[j] = b_x;
    }
    for (int i = 0; i < nx; i++){
      left_edge_x[i] = left_bnd_x[id_yz(nz-1,0)] + i*distance_x_left;
      right_edge_x[i] = left_bnd_x[id_yz(nz-1,ny-1)] + i*distance_x_right;
    }
    // left right bottom top
    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x,  top_edge_x, left_edge_x, right_edge_x,  back_bnd_x); //interpolate back_bnd_x

    /** front_bnd_x:
     *
     *   left_bnd_x(0,ny-1)         
     *               x--------x       ^
     *               )        |       |
     *                (       |       |
     * left_bnd_x(0,:) )      | right_vertice_x
     *                  (     |       |
     *                   )    |       |
     *                    x---x       v
     *      left_bnd_x(0,0)         
     **/    

    distance_x_left=(right_vertice[0]-left_bnd_x[id_yz(0,0)])/(nx-1);
    distance_x_right=(right_vertice[0]-left_bnd_x[id_yz(0,ny-1)])/(nx-1);

    for (int j = 0; j < ny; j++){
      top_edge_x[j]   = right_vertice[0];
      bottom_edge_x[j]= left_bnd_x[id_yz(0,j)];
    }
    for (int i = 0; i < nx; i++){
      left_edge_x[i] = left_bnd_x[id_yz(0,0)] + i*distance_x_left;
      right_edge_x[i] = left_bnd_x[id_yz(0,ny-1)] + i*distance_x_right;
    }
    // left right bottom top
    transFiniteInterpolation_singleCoordinate(nx, ny, bottom_edge_x, top_edge_x, left_edge_x, right_edge_x, front_bnd_x);
  }

  for(int j = 0 ; j< ny; j++){
    for(int i = 0 ; i< nx; i++){
      front_bnd_z[id_xy(j,i)]= left_vertice[2];  //constant in z
      back_bnd_z[id_xy(j,i)] = right_vertice[2]; //constant in z
    }
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

  double mesh_size_x = 1.0/(mx-1);
  double mesh_size_y = 1.0/(my-1);
  double mesh_size_z = 1.0/(mz-1);

  int i_0;
  int j_0;
  int k_0;    

  double u,v,w;
  double uv,vw,uw,uvw1,uvw2;

  double q,r,s;

  for(int i=i_m ; i< i_p ; i++){
    i_0=i-i_m;
    q=mesh_size_x * i;
    for(int j=j_m ; j< j_p ; j++){
      j_0=j-j_m;
      r=mesh_size_y *j;
      for(int k=k_m ; k< k_p ; k++){
	k_0=k-k_m;
	s=mesh_size_z * k;
	
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

      }
    }
  }
}



/**Transfinite interpolation
 * Returns the interpolated coordiantes on a uniform mesh of size mx*my
 * in a planar spanned by boundary curves left_bnd right_bnd bottom_bnd top_bnd
 *
 *               top_bnd
 *  top_bnd(0)             top_bnd(mx-1)          
 *      (0,my)o~~o~~o~~o~~o(mx,my)
 *            (  :  :  :  (
 *            o--o--o--o--o
 *            )  :  :  :  )
 *   left_bnd o--o--o--o--o right_bnd
 *            (  :  :  :  (
 *            o--o--o--o--o
 *            )  :  :  :  )
 *       (0,0)o~~o~~o~~o~~o(mx,0)
 * left_bnd(0)             right_bnd(0)
 *             bottom_bnd 
**/
void CurvilinearTransformation::transFiniteInterpolation_singleCoordinate(int mx, int my, double* left_bnd, double* right_bnd, double* bottom_bnd, double* top_bnd, double* curvilinear ){

   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);
   
   double r;
   double q;
 
   kernels::idx2 id_xy(my,mx);
   
   for(int j =0 ; j < my ; j++) {
     for(int i = 0 ; i < mx ; i++) {
       q = (i)*mesh_size_x;
       r = (j)*mesh_size_y;

//From https://en.wikipedia.org/wiki/Transfinite_interpolation
       curvilinear[id_xy(j,i)]=
	 (1-q)*left_bnd[j]+q*right_bnd[j]
	 +(1-r)*bottom_bnd[i]+r*top_bnd[i]
	 -(1-q)*(1-r)*left_bnd[0] //intersection left bottom
	 -q*(1-r)*right_bnd[0]    //intersection right bottom
	 -r*(1-q)*top_bnd[0]      //intersection left top
	 -(r*q)*top_bnd[mx-1];    //intersection right top
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


/**
 * Evaluates the ith 1d lagrange polynomial on a uniform mesh of degree _num_nodes at x 
 **/
double CurvilinearTransformation::lagrangeBasis_uniform(double x,int i){
  double result=1;

  for (int j = 0 ; j< i ; j ++){
    // unifMesh
    result *= (x-unif_mesh[j]);
  }
  
  for (int j = i+1 ; j < _num_points ; j ++){
    // 
    result *= (x-unif_mesh[j]);
  }

  result*= denominator_lagrange[i];
  
  return result;
}

/**
 * Evaluates the ith 1d lagrange polynomial on the set of nodes of length n at x 
 **/
double CurvilinearTransformation::lagrangeBasis(double x,int i,int n, double* nodes){
  double result=1;
  double denominator=1;

  for (int j = 0 ; j< i ; j ++){
    // unifMesh
    result *= (x-nodes[j]);
    denominator *= (nodes[i]-nodes[j]);;
  }
  for (int j = i+1 ; j < n ; j ++){
    result *= (x-nodes[j]);
    denominator *= (nodes[i]-nodes[j]);
  }
  result/= denominator;
  return result;
}

double  CurvilinearTransformation::interpolate2D(double const x, double const y,
						 int const number_nodes_x, int const number_nodes_y,
						 double* const nodes_x, double* const nodes_y,
						 const double* const values){
  double result=0.;
  kernels::idx2 id_yz(number_nodes_y,number_nodes_x);

  double a_x,a_y;
    
  for (int k = 0; k< number_nodes_y; k++){
    for (int j = 0; j <number_nodes_x; j++){
      a_x = lagrangeBasis(x, j, number_nodes_x, nodes_x);
      a_y = lagrangeBasis(y, k, number_nodes_y, nodes_y);
      result += a_x*a_y*values[id_yz(k,j)];
    }
  }

  return result;
}  


void CurvilinearTransformation::interpolate3D_uniform(int x, int y ,int z, double* values, double& result){
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
	
	result += a_x*a_y*a_z*values[id_xyz(k,j,i)];
      }
    }
  }
}



void CurvilinearTransformation::interpolate2D_uniform(double x, double y,  double* values, double& result){
  double a_x=0;
  double a_y=0;
  result=0;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);
  
  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
      a_x=lagrangeBasis_uniform(x,i,num_nodes);
      a_y=lagrangeBasis_uniform(y,j,num_nodes);
      result += values[id_xy(j,i)] * a_x*a_y;
    }
  }
}



void CurvilinearTransformation::getValuesAtQuadNodes3D(double* dest_mesh, int num_nodes, double* results){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);
  
  for (int k = 0 ; k< num_nodes ; k ++){
    for (int j = 0 ; j< num_nodes ; j ++){
      for (int i = 0 ; i< num_nodes ; i ++){
	interpolate3D_uniform(i,j,k,dest_mesh,num_nodes,results[id_xyz(k,j,i)]);
      }
    }
  }
}



void CurvilinearTransformation::getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results){
  kernels::idx2 id_xy(num_nodes,num_nodes);

  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
#if defined(_GLL)
      interpolate(kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-i],kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
      #else
      interpolate(kernels::gaussLegendreNodes[num_nodes-1][i],kernels::gaussLegendreNodes[num_nodes-1][j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
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
      
      jacobian[id_xy(j,i)]=x_der_x*y_der_y-x_der_y*y_der_x;
      
      q_x[id_xy(j,i)]=y_der_y/jacobian[id_xy(j,i)];
      q_y[id_xy(j,i)]=-x_der_y/jacobian[id_xy(j,i)];
      r_x[id_xy(j,i)]=-y_der_x/jacobian[id_xy(j,i)];
      r_y[id_xy(j,i)]=x_der_x/jacobian[id_xy(j,i)];
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


const double CurvilinearTransformation::interpolate_fault_surface( const double* const fault_x, const double* const fault_y, const double* const fault_z,
								   const double y, const double z, double& x){ //TODO: I don't like the siganture

  kernels::idx2 id_yz(_nz,_ny);

  //mesh width
  double dy=fault_y[id_yz(0,1)]-fault_y[id_yz(0,0)];
  double dz=fault_z[id_yz(1,0)]-fault_z[id_yz(0,0)];
  
  dr = std::sqrt(pow(dy,2) + pow(dz, 2));
    
  stencil_size=2;
  ndp  =5;

  int k,j

  //find first point within mesh witdh near to (y,z)
  for (int i_k = 0; i_k<_nz; i_k++){
    for (int i_j = 0; i_j<ny; i_j++){
      if (std::abs(y-fault_y[id_yz(i_k,i_j)]) <= dy) { //TODO: Check if radius really war unnecessary
	if (std::abs(z-fault_z[id_yz(i_k,i_j)]) <= dz) {
	  k=i_k; j=i_j; break;
	}
      }
    }
  }


  //select stencil around that point
  //if point is at boundary select stencil along boundary
  if(j < stencil){
    nmy = 0;
    npy = stencil_size*2+1
  }
  else if(j > my-stencil){
    nmy = my-stencil_size*2-1;
    npy = my; 
  }else{
    nmy = j-stencil_size;
    npy = j+stencil_size+1;
  }

  if (k <= 3){
    nmz = 0;
    npz = stencil_size*2+1;
  }        
  else if (k => mz-3){
    nmz = mz-stencil_size*2-1;
    npz = mz;
  }else{
    nmz = k-stencil_size;
    npz = k+stencil_size+1; 
  }


  //extract sctencil nodes
  double nodes_y[stencil_size*2+1];
  double nodes_z[stencil_size*2+1];
  for (int i_j = nmy; i_j<npy; i_j++){
    nodes_y[i_j] = fault_y[id_yz(k,i_j)];
  }
  for (int i_k = 0; i_k<mz; i_k++){
    nodes_z[i_k] = fault_z[id_yz(i_k,j)];
  }

  // 2d interpolation along stencil
  x=interpolate_2d(stencil_size*2+1, stencil_size*2+1, nodes_y, nodes_z, fault_x+id_xy(nmz,nmx));
}


int CurvilinearTransformation::getBlock(const tarch::la::Vector<DIMENSIONS,double>& center,
					const tarch::la::Vector<DIMENSIONS,double>& dx){

  //check four outer elements
  if(center[0] < _dx){
    return 0;
  }else if(center[0] > _right_vertice[0]-_dx){
    return 1;
  }
  
  //check for fault position
  return center[0] > _fault_position ? 1 : 0;
}


  


// void CurvilinearTransformation::getBoundaryCurves3D(int num_points,
// 			 double offset_x, double offset_y, double offset_z,
// 			 double width_x, double width_y , double width_z ,
// 			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
// 			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
// 			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
// 			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
// 			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
// 			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
//   double pi = 3.14159265359;

//   int num_elt_x= std::ceil(1/width_x);
//   int num_elt_y= std::ceil(1/width_y);
//   int num_elt_z= std::ceil(1/width_z);

//   int nx = num_elt_x*(num_points-1) + 1;
//   int ny = num_elt_y*(num_points-1) + 1;
//   int nz = num_elt_z*(num_points-1) + 1;  
  
//   double dx= width_x/(num_points-1);
//   double dy= width_y/(num_points-1);
//   double dz= width_z/(num_points-1);

//   kernels::idx2 id_xy(ny,nx);// back front
//   kernels::idx2 id_xz(nz,nx);// bottom top
//   kernels::idx2 id_yz(nz,ny);//left right


//   //top bottom
//   for(int k = 0 ; k < nz ; k ++){
//     for(int i = 0 ; i < nx ; i ++){
//       bottom_bnd_x[id_xz(k,i)] = dx*i;
//       bottom_bnd_z[id_xz(k,i)] = dz*k;
//       bottom_bnd_y[id_xz(k,i)] = 1;	    

//       top_bnd_x[id_xz(k,i)] = dx*i;
//       top_bnd_z[id_xz(k,i)] = dz*k;
//       top_bnd_y[id_xz(k,i)] = 0+0.25*std::sin(2*pi*top_bnd_x[id_xz(k,i)])*std::sin(2*pi*top_bnd_z[id_xz(k,i)]);

// 	}
//       }
    
  
//   //left right
//   for(int k = 0 ; k < nz ; k ++){
//     for(int j = 0 ; j < ny ; j ++){
//       left_bnd_x[id_xz(k,j)]  = 0;	    
//       left_bnd_y[id_xz(k,j)]  = dy*j;
//       left_bnd_z[id_xz(k,j)]  = dz*k;
      
//       right_bnd_x[id_xz(k,j)] = 1.0;
//       right_bnd_y[id_xz(k,j)] = dy*j;
//       right_bnd_z[id_xz(k,j)] = dz*k;
      
//     }
//   }
  
//   //front back
//   for(int j = 0 ; j < ny ; j ++){
//     for(int i = 0 ; i < nx ; i ++){
      
//       front_bnd_x[id_xy(j,i)]  = dx*i;
//       front_bnd_y[id_xy(j,i)]  = dy*j;
//       front_bnd_z[id_xy(j,i)]  = 0;
      
//       back_bnd_x[id_xy(j,i)] = dx*i;
//       back_bnd_y[id_xy(j,i)] = dy*j;
//       back_bnd_z[id_xy(j,i)] = 1;
//     }
//   }
    
// }

// void CurvilinearTransformation::getBoundaryCurves3D_fixedTopFace(int num_points,
// 				      double offset_x, double offset_y, double offset_z,
// 				      double width_x, double width_y , double width_z ,
// 				      double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
// 				      double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
// 				      double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
// 				      double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
// 				      double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
// 				      double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){
  
//   double pi = 3.14159265359;


//   int nx = std::round(1/width_x)*(num_points-1) + 1;
//   int ny = std::round(1/width_y)*(num_points-1) + 1;
//   int nz = std::round(1/width_z)*(num_points-1) + 1;  

//   double dx= 1.0/(nx-1);
//   double dz= 1.0/(nz-1);
//   double dy= 1.0/(ny-1);

//   double depth=1.0;

//   kernels::idx2 id_xy(ny,nx); // back front
//   kernels::idx2 id_xz(nz,nx); // botton top
//   kernels::idx2 id_yz(nz,ny); //left right
  
//   //given top surface
//   for(int k = 0 ; k< nz; k++){  
//     for(int i = 0 ; i< nx; i++){
//       top_bnd_x[id_xz(k,i)] = 0+dx*i;      
//       top_bnd_z[id_xz(k,i)] = 0+dz*k;
//       top_bnd_y[id_xz(k,i)] = 0- 0*(0.1*(top_bnd_x[id_xz(k,i)] + top_bnd_z[id_xz(k,i)])
// 			       + 0.25*(std::sin(4*pi*top_bnd_x[id_xz(k,i)]+3.34)*std::cos(4*pi*top_bnd_x[id_xz(k,i)])
// 				       * std::sin(4*pi*top_bnd_z[id_xz(k,i)]+3.34)*std::cos(4*pi*top_bnd_z[id_xz(k,i)])));

//     }
//   }
  
  
//   for(int i = 0 ; i< nx; i++){
//     for(int k = 0 ; k< nz; k++){
//       bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
//       bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
//       bottom_bnd_y[id_xz(k,i)] = depth;      
//     }
//   }


//   getInterpolatedFace_fromBottomAndTop(nx,ny,nz,0,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       left_bnd_x,left_bnd_y,left_bnd_z);
    
//   getInterpolatedFace_fromBottomAndTop(nx,ny,nz,1,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       right_bnd_x,right_bnd_y,right_bnd_z);
  
//   getInterpolatedFace_fromBottomAndTop(nz,ny,nx,2,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       front_bnd_x,front_bnd_y,front_bnd_z);

//   getInterpolatedFace_fromBottomAndTop(nz,ny,nx,3,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       back_bnd_x,back_bnd_y,back_bnd_z);
// }


// //TODO: not working with new width definition !
// void CurvilinearTransformation::getBoundaryCurves3D_fixedTopFace_forBlock(int num_points,
// 					       int nx, int ny, int nz, int n,
// 					       double width_x, double width_y , double width_z,	      
// 					       double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
// 					       double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
// 					       double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
// 					       double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
// 					       double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
// 					       double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
//   double pi = 3.14159265359;


//   double dx= 10*width_x/(nx-1);
//   double dz= 10*width_y/(nz-1);
//   double dy= 10*width_z/(ny-1);

//   double depth=b_y;

//   double x;
//   double z;
//   double y;

//   kernels::idx2 id_xy(ny,nx); // back front
//   kernels::idx2 id_xz(nz,nx); // botton top
//   kernels::idx2 id_yz(nz,ny); // left right

  
//   //given top surface
//   for(int k = 0 ; k< nz; k++){  
//     for(int i = 0 ; i< nx; i++){

//       top_bnd_y[id_xz(k,i)] = 0;      
//       top_bnd_z[id_xz(k,i)] = 0+dz*k;

     
      
//       if(n == 0){
//       top_bnd_x[id_xz(k,i)] = 0+dx*i;
      
//       }else{
// 	top_bnd_x[id_xz(k,i)] = 5 + dx*i;
	
//       }

//       x = top_bnd_x[id_xz(k,i)];
//       z = top_bnd_z[id_xz(k,i)];

//       top_bnd_y[id_xz(k,i)] -= topography(x, z, a_x, b_x, a_z, b_z,  depth);
//     }
//   }
  
//   for(int k = 0 ; k< nz; k++){
//     for(int i = 0 ; i< nx; i++){
//       bottom_bnd_y[id_xz(k,i)] = depth;      
//       bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
//       bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      
//     }
//   }


//   int n_left_right=ny;

  
//   int n_block;
//   int n_top_bottom;

//   //left face
//   n_block = nx;
//   n_top_bottom = nz;

//                                      //n_block, n_left_right, n_top_bottom,
//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,0,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       left_bnd_x,left_bnd_y,left_bnd_z);

//   //right face
//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,1,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       right_bnd_x,right_bnd_y,right_bnd_z);


//   //front face
//   n_block = nz;
//   n_top_bottom = nx;
  
//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,2,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       front_bnd_x,front_bnd_y,front_bnd_z);

//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,3,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       back_bnd_x,back_bnd_y,back_bnd_z);

  
// }

// void CurvilinearTransformation::getBoundaryCurves(int num_points,double offset_x, double offset_y,double width_x, double width_y ,double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y){

//   double pi = 3.14159265359;

//   int nx = 1/width_x*num_points;
//   int ny = 1/width_y*num_points;

//   double dx= 1.0/(nx-1);
//   double dy= 1.0/(ny-1);

//   // int i_0 =  (offset_x-0.0)/width_x;
//   // int j_0 =  (offset_y-0.0)/width_y;
  

//   //std::cout<<i_0<<std::endl;
//   //std::cout<<j_0<<std::endl;
//   //std::cout<<std::endl;

  
//   // for(int i = 0 ; i< num_points; i++){
//   //   left_bnd_x[i] =  offset_x;
//   //   right_bnd_x[i] = width_x+offset_x;
//   //   bottom_bnd_x[i] = width_x*dx*i + offset_x;
//   //   top_bnd_x[i] = width_x*dx*i + offset_x;

//   //   left_bnd_y[i] = width_y*dy*i + offset_y;
//   //   right_bnd_y[i] = width_y*dy*i + offset_y;
//   //   bottom_bnd_y[i] =offset_y;
//   //   top_bnd_y[i] = width_y+offset_y + 0.0*std::sin(2*pi*top_bnd_x[i]);
//   // }


//   for(int i = 0 ; i< nx; i++){
    
//     bottom_bnd_x[i] = 0 + dx*i;
//     top_bnd_x[i]    = 0 + dx*i;
    
//     bottom_bnd_y[i] = 0.0;
//     top_bnd_y[i]    = 1 + 0.0*std::sin(2*pi*top_bnd_x[i]);

//     //std::cout<<top_bnd_x[i]<<std::endl;
//     //std::cout<<i<<std::endl;
//     //std::cout<<std::endl;
//   }

//   //std::exit(-1);
  
//   for(int i = 0 ; i< ny; i++){
    
//     left_bnd_x[i]  = 0;
//     right_bnd_x[i] = 1;

//     left_bnd_y[i]  = dy*i;
//     right_bnd_y[i] = dy*i;
//   }


//   double h0y = (top_bnd_y[0] - bottom_bnd_y[0])/(ny-1);
//   double hny = (top_bnd_y[nx-1] - bottom_bnd_y[nx-1])/(ny-1);

//   for(int i = 0 ; i< ny; i++){

//     left_bnd_y[i]  = bottom_bnd_y[0] + h0y*i;
//     right_bnd_y[i] = bottom_bnd_y[nx-1] + hny*i;
    
//   }
  
// }
