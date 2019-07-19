#include "kernels/KernelUtils.h"

#if defined(_GLL)
#include "kernels/GaussLobattoBasis.h"
#else
#include "kernels/GaussLegendreBasis.h"
#endif
#include <cmath>
#include <string>
#include "netcdf.h" 

class MeshReader{

 public:
 MeshReader(std::string filename, int num_nodes):
  m_num_nodes(num_nodes), m_filename(filename){
    
    int retval;
    if ((retval = nc_open(m_filename.c_str(), NC_NOWRITE, &m_ncid))){
      std::cout<< "Severe: Couldn't open "<< m_filename << " " << retval << std::endl;
    }

    if ((retval = nc_get_att_int(m_ncid,NC_GLOBAL,"lod_min",&m_min_lod))){
      std::cout<< "Severe: Couldn't read attribute lod_min "<< retval << std::endl;
    }

    if ((retval = nc_get_att_int(m_ncid,NC_GLOBAL,"lod_max",&m_max_lod))){
      std::cout<< "Severe: Couldn't read attribute lod_max "<< retval << std::endl;
    }

    int order;

    if ((retval = nc_get_att_int(m_ncid,NC_GLOBAL,"order",&order))){
      std::cout<< "Severe: Couldn't read attribute order "<< retval << std::endl;
    }

    if(num_nodes-1 != order){
      std::cout << "File is made for order " << order << " not " << num_nodes-1 << std::endl;
    }
    


    x_mesh_id = new int[m_max_lod-m_min_lod +1];
    y_mesh_id = new int[m_max_lod-m_min_lod +1];
    z_mesh_id = new int[m_max_lod-m_min_lod +1];
    
    for(int lod = m_min_lod; lod< m_max_lod+1; lod++){
      std::string x_mesh_str = "curvilinear_x_";
      x_mesh_str += std::to_string(lod);
      std::string y_mesh_str = "curvilinear_y_";
      y_mesh_str += std::to_string(lod);
      std::string z_mesh_str = "curvilinear_z_";
      z_mesh_str += std::to_string(lod);
      
      if ((retval = nc_inq_varid(m_ncid, x_mesh_str.c_str(), &x_mesh_id[lod-m_min_lod]))){
	std::cout<< "Severe: Couldn find Variable "<< x_mesh_str << " " << retval << std::endl;
      }
      if ((retval = nc_inq_varid(m_ncid, y_mesh_str.c_str(), &y_mesh_id[lod-m_min_lod]))){
	std::cout<< "Severe: Couldn find Variable "<< y_mesh_str << " " << retval << std::endl;
      }

      if ((retval = nc_inq_varid(m_ncid, z_mesh_str.c_str(), &z_mesh_id[lod-m_min_lod]))){
	std::cout<< "Severe: Couldn find Variable "<< z_mesh_str << " " << retval << std::endl;
      }
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
	  =lagrangeBasis(kernels::lobatto::nodes[num_nodes-1][num_nodes-j-1],i,num_nodes);
#else
	lagrange_basis_at_nodes[id_xy(j,i)]
	  =lagrangeBasis(kernels::legendre::nodes[num_nodes-1][j],i,num_nodes);
#endif
      }
    }
  }

  ~MeshReader(){
    nc_close(m_ncid);
    delete[] unif_mesh;
    delete[] denominator_lagrange;
    delete[] lagrange_basis_at_nodes;
  }


  void metricDerivativesAndJacobian3D(int i_m,int j_m,int k_m,double dx,
						  double* gl_vals_x, double* gl_vals_y, double* gl_vals_z,
						  double* q_x, double* q_y, double* q_z,
						  double* r_x, double* r_y, double* r_z,
						  double* s_x, double* s_y, double* s_z,
						  double* jacobian
						  ){

  double curvilinear_x[m_num_nodes*m_num_nodes*m_num_nodes];
  double curvilinear_y[m_num_nodes*m_num_nodes*m_num_nodes];
  double curvilinear_z[m_num_nodes*m_num_nodes*m_num_nodes];
  int lod   = std::round(std::log(1./dx)/std::log(3.));

  // number of colocal points
  int num_nodes_domain = std::round(1.0/dx) * (m_num_nodes-1) +1;
  
  size_t start[3]= {k_m,j_m,i_m};
  size_t count[3]= {m_num_nodes,m_num_nodes,m_num_nodes};

  int retval;
  if((retval=nc_get_vara_double(m_ncid, x_mesh_id[lod-m_min_lod], start ,count, curvilinear_x))){
    std::cout<< "Severe: Couldn't read x mesh " << retval << std::endl;
  }
  if((retval=nc_get_vara_double(m_ncid, y_mesh_id[lod-m_min_lod], start ,count, curvilinear_y))){
    std::cout<< "Severe: Couldn't read y mesh " << retval << std::endl;
  }
  std::cout << y_mesh_id[lod-m_min_lod] << std::endl;
  if((retval=nc_get_vara_double(m_ncid, z_mesh_id[lod-m_min_lod], start ,count, curvilinear_z))){
    std::cout<< "Severe: Couldn't read z mesh " << retval << std::endl;
  }


  /* for(int i = 0 ; i< m_num_nodes*m_num_nodes*m_num_nodes;i++){ */
  /*   std::cout << curvilinear_x[i] << std::endl; */
  /* } */
  
  getValuesAtQuadNodes3D(curvilinear_x,m_num_nodes,gl_vals_x);
  getValuesAtQuadNodes3D(curvilinear_y,m_num_nodes,gl_vals_y);
  getValuesAtQuadNodes3D(curvilinear_z,m_num_nodes,gl_vals_z);

  /* for(int i = 0 ; i< m_num_nodes*m_num_nodes*m_num_nodes;i++){ */
  /*   std::cout << gl_vals_x[i] << std::endl; */
  /* } */

  double x_der_x; 
  double x_der_y;
  double x_der_z;  
  
  double y_der_x; 
  double y_der_y;
  double y_der_z;

  double z_der_x; 
  double z_der_y;
  double z_der_z;  

  
  kernels::idx3 id_xyz(m_num_nodes,m_num_nodes,m_num_nodes);
  
  for(int k = 0 ; k < m_num_nodes ; k ++){
    for(int j = 0 ; j < m_num_nodes ; j ++){
      for(int i = 0 ; i< m_num_nodes ; i ++){

  	computeDerivatives_x_3D(i,j,k,gl_vals_x,m_num_nodes,x_der_x, dx);
  	computeDerivatives_y_3D(i,j,k,gl_vals_x,m_num_nodes,x_der_y, dx);
  	computeDerivatives_z_3D(i,j,k,gl_vals_x,m_num_nodes,x_der_z, dx);

	computeDerivatives_x_3D(i,j,k,gl_vals_y,m_num_nodes,y_der_x, dx);
	computeDerivatives_y_3D(i,j,k,gl_vals_y,m_num_nodes,y_der_y, dx);
	computeDerivatives_z_3D(i,j,k,gl_vals_y,m_num_nodes,y_der_z, dx);

	computeDerivatives_x_3D(i,j,k,gl_vals_z,m_num_nodes,z_der_x, dx);
	computeDerivatives_y_3D(i,j,k,gl_vals_z,m_num_nodes,z_der_y, dx);
	computeDerivatives_z_3D(i,j,k,gl_vals_z,m_num_nodes,z_der_z, dx);
	
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
};


  
 private:

 std::string m_filename="topography.nc";

 int m_min_lod;
 int m_max_lod;
 int m_num_nodes;
 int m_ncid;

 /*Variable ids*/
 int* lat_var_id;
 int* lon_var_id;
 int* x_mesh_id;
 int* y_mesh_id;
 int* z_mesh_id;

 
 double* lagrange_basis_at_nodes;
 double* unif_mesh;
 double* denominator_lagrange;


 void getValuesAtQuadNodes3D(double* dest_mesh, int num_nodes, double* results){
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

 void interpolate3D(int x, int y ,int z, double* dest_mesh, int num_nodes,double& result){

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

 void computeDerivatives_x_3D(int i, int j , int k, double* values , int num_nodes, double& der_x, double dx){
  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);

  der_x = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xyz(k,j,n)]/dx;
  }

}

void computeDerivatives_y_3D (int i, int j,int k, double* values , int num_nodes, double& der_y, double dy){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes); 

  der_y = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xyz(k,n,i)]/dy;
  }
}

void computeDerivatives_z_3D (int i, int j , int k ,double* values , int num_nodes, double& der_z, double dz){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes); 

  der_z = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_z += kernels::dudx[num_nodes-1][k][n] * values[id_xyz(n,j,i)]/dz;
  }

}

double lagrangeBasis(double x,int i,int num_points){
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


};
