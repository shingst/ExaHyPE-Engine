#ifndef __CurvilinearTransformation_CLASS_HEADER__
#define __CurvilinearTransformation_CLASS_HEADER__
#include <array>

#include "tarch/la/Vector.h"
#include "kernels/aderdg/generic/Kernels.h"
struct Boundary_single_coordinate{
  std::vector<double> left;
  std::vector<double> right;
  std::vector<double> bottom;
  std::vector<double> top;
  std::vector<double> front;
  std::vector<double> back;
  Boundary_single_coordinate(int nx, int ny, int nz){
    left  .resize(ny*nz);
    right .resize(ny*nz);
    bottom.resize(nx*nz);
    top   .resize(nx*nz);
    front .resize(nx*ny);
    back  .resize(nx*ny);
  }
  ~Boundary_single_coordinate(){
  }

  friend std::ostream& operator<< (std::ostream& stream, const Boundary_single_coordinate& boundary){

    stream << std::endl;
    stream << "Left" << std::endl;
    for(int i = 0 ; i<boundary.left.size(); i++){
      stream << boundary.left[i] << std::endl;
    }

    stream << "Right" << std::endl;
    for(int i = 0 ; i<boundary.right.size(); i++){
      stream << boundary.right[i] << std::endl;
    }

    stream << "Bottom" << std::endl;
    for(int i = 0 ; i<boundary.bottom.size(); i++){
      stream << boundary.bottom[i] << std::endl;
    }

    stream << "Top" << std::endl;
    for(int i = 0 ; i<boundary.top.size(); i++){
      stream << boundary.top[i] << std::endl;
    }

    stream << "Front" << std::endl;
    for(int i = 0 ; i<boundary.front.size(); i++){
      stream << boundary.front[i] << std::endl;
    }

    stream << "Back" << std::endl;
    for(int i = 0 ; i<boundary.back.size()-1; i++){
      stream << boundary.back[i] << std::endl;
    }

    return stream << boundary.back[boundary.back.size()-1] << std::endl;
  }


};

class CurvilinearTransformation{

 public:
 //Default Constructor
 CurvilinearTransformation():_num_nodes(0), _dx(0),
    _fault_position(0),
    _right_vertex{0.},
    _left_vertex{0.}{};

    CurvilinearTransformation(const int num_nodes, const int mesh_level ,
			      const double fault_position,
			      double* domain_offset,
			      double* domain_size);

    
   
 //Destructor
 ~CurvilinearTransformation(){
   delete[] _boundary_x;
   delete[] _boundary_y;
   delete[] _boundary_z;
 };


  void genCoordinates(const tarch::la::Vector<DIMENSIONS,double>& center,
		      const tarch::la::Vector<DIMENSIONS,double>& dx,    
		      double* gl_vals_x,double* gl_vals_y,double* gl_vals_z,
		      double* jacobian,
		      double* q_x,double* q_y,double* q_z,
		      double* r_x,double* r_y,double* r_z,
		      double* s_x,double* s_y,double* s_z);

  void metricDerivativesAndJacobian3D(const double* const curvilinear_x,
				      const double* const curvilinear_y,
				      const double* const curvilinear_z,
				      const double* const dx,
				      double* gl_vals_x, double* gl_vals_y, double* gl_vals_z,
				      double* q_x, double* q_y, double* q_z,
				      double* r_x, double* r_y, double* r_z,
				      double* s_x, double* s_y, double* s_z,
				      double* jacobian);




 private:

  int _num_nodes;
  double _dx;
  double _block_width_x[2];
  double _fault_position;
  std::array<double,DIMENSIONS> _left_vertex;
  std::array<double,DIMENSIONS> _right_vertex;
  std::array<double,DIMENSIONS> _rect_width;

  //number of elements over whole domain
  int _ne_x_block[2];
  int _ne_x;
  int _ne_y;
  int _ne_z;

  //number of nodes over whole domain
  int _nx_block[2];
  int _nx;
  int _ny;
  int _nz;

  //boundaries
  Boundary_single_coordinate* _boundary_x[2];
  Boundary_single_coordinate* _boundary_y[2];
  Boundary_single_coordinate* _boundary_z[2];

  double* lagrange_basis_at_nodes;
  double* denominator_lagrange;
  double* unif_mesh;

  double fault(double y, double z);
  double topography(double x, double z, double depth);
  
/* #if defined(USE_ASAGI) */
/*   double topography_fromASAGI(double x, double z, double* topography, easi::ArraysAdapter& adapter, easi::Component* model); */
/* #endif   */

  void getInterpolatedFace_fromBottomAndTop( int n_block, int n_left_right, int n_top_bottom, int face,
					   double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					   double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					   double* face_x, double* face_y , double* face_z);


  double interpolate2D(double x, double y,
		       int number_nodes_x, int number_nodes_y,
		       const double* const nodes_x, const double* const nodes_y,
		       const double* const values);

  double interpolate2D_uniform(double x, double y,  double* values);
  double interpolate3D_uniform(int x, int y ,int z, const double* const values);

  double lagrangeBasis(double x, int i, int n, const double* const nodes);
  double lagrangeBasis_uniform(double x,int i);

  void getValuesAtQuadNodes2D(double* dest_mesh,  double* results);
  void getValuesAtQuadNodes3D(const double* const dest_mesh, double* results);

  void computeDerivatives_x(int i, int j, double* values , int num_nodes, double& der_x, double dx);
  void computeDerivatives_y(int i, int j, double* values , int num_nodes, double& der_y, double dy);

  void computeDerivatives_x_3D(int i, int j, int k, double* values , int num_nodes, double& der_x, double dx);
  void computeDerivatives_y_3D(int i, int j, int k, double* values , int num_nodes, double& der_y, double dy);
  void computeDerivatives_z_3D(int i, int j, int k ,double* values , int num_nodes, double& der_z, double dz);


  void metricDerivativesAndJacobian(int num_nodes, double* curvilinear_x, double* curvilinear_y,double* gl_vals_x,double* gl_vals_y,double* q_x,double* q_y,double* r_x,double* r_y,double* jacobian, double dx, double dy);

  double interpolate_fault_surface(const double* const fault_x,
				   const double* const fault_y,
				   const double* const fault_z,
				   const double y, const double z);
  
  double interpol2d_dG(int my, int mz, int nmy, int npy, int nmz, int npz, double y, double z, double* yj, double* zk, double* f);
  double lagrange(int m, int p, int i, double x, double *xi);

  void getBoundaryCurves3D_cutOffTopography_withFault(int n,
						      Boundary_single_coordinate* boundary_x,
						      Boundary_single_coordinate* boundary_y,
     						      Boundary_single_coordinate* boundary_z);


  void transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y );



  void transFiniteInterpolation_singleCoordinate(int mx, int my,
						 double* left_bnd, double* right_bnd,
						 double* bottom_bnd, double* top_bnd,
						 double* curvilinear);

  void transFiniteInterpolation_singleCoordinate(int mx, int my,
						 std::vector<double> left_bnd,
						 std::vector<double> right_bnd,
						 std::vector<double> bottom_bnd,
						 std::vector<double> top_bnd,
						 double* curvilinear){
    transFiniteInterpolation_singleCoordinate(mx,  my,
					  left_bnd.data(),  right_bnd.data(),
					  bottom_bnd.data(),  top_bnd.data(),
					      curvilinear);
      };

  void transFiniteInterpolation3D(int n,
				  int k_m, int k_p ,
				  int j_m, int j_p ,
				  int i_m, int i_p ,
				  Boundary_single_coordinate* boundary,
				  double* curvilinear);
  
  int getBlock(const tarch::la::Vector<DIMENSIONS,double>& center,
	       const tarch::la::Vector<DIMENSIONS,double>& dx);

};
#endif
