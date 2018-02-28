#ifndef __CurvilinearTransformation_CLASS_HEADER__
#define __CurvilinearTransformation_CLASS_HEADER__

struct Boundary_single_coordiante{
  double* left;
  double* right;
  double* bottom;
  double* top;
  double* front;
  double* back;
  Boundary_single_coordiante(int nx, int ny, int nz){
    left  =new double[ny*nz];
    right =new double[ny*nz];
    bottom=new double[nx*nz];
    top   =new double[nx*nz];
    front =new double[nx*ny];
    back  =new double[nx*ny];
  }
}

class CurvilinearTransformation{

 public:
  //CurvilinearTransformation(){};

  CurvilinearTransformation(const int a_num_nodes, const double a_dx);

  void genCoordinates(const tarch::la::Vector<DIMENSIONS,double>& center,
		      const tarch::la::Vector<DIMENSIONS,double>& dx,    
		      double* gl_vals_x,double* gl_vals_y,double* gl_vals_z,
		      double* jacobian,
		      double* q_x,double* q_y,double* q_z,
		      double* r_x,double* r_y,double* r_z,
		      double* s_x,double* s_y,double* s_z);
  

  
  void getBoundaryCurves(int num_points,double offset_x, double offset_y,
			 double width_x, double width_y ,
			 double* left_bnd_x, double* left_bnd_y,
			 double* right_bnd_x, double* right_bnd_y,
			 double* bottom_bnd_x, double* bottom_bnd_y,
			 double* top_bnd_x, double* top_bnd_y);


  void getBoundaryCurves3D(int num_points,
			 double offset_x, double offset_y, double offset_z,
			 double width_x, double width_y , double width_z ,
			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z);

  void getBoundaryCurves3D_fixedTopFace(int num_points,
			double offset_x, double offset_y, double offset_z,
			double width_x, double width_y , double width_z ,
			double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
			double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
			double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
			double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
			double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
			double* back_bnd_x, double* back_bnd_y, double* back_bnd_z);


  void getBoundaryCurves3D_fixedTopFace_forBlock(int num_points,
			 int nx, int ny, int nz, int n,
			 double width_x, double width_y , double width_z,
			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z);


  void getBoundaryCurves3D_cutOffTopography_withFault(int num_points,
		      int nx, int ny, int nz, int n,double fault_position,
		      double width_x, double width_y , double width_z,
		      double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
		      double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
		      double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
		      double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
		      double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
		      double* back_bnd_x, double* back_bnd_y, double* back_bnd_z);

//void transFiniteInterpolation( int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y );

  void transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p,
				int num_points,
				double* left_bnd_x, double* left_bnd_y,
				double* right_bnd_x, double* right_bnd_y,
				double* bottom_bnd_x, double* bottom_bnd_y,
				double* top_bnd_x, double* top_bnd_y,
				double* curvilinear_x , double* curvilinear_y );

  void transFiniteInterpolation_singleCoordinate(int mx, int my,
  				double* left_bnd, double* right_bnd,
  				double* bottom_bnd, double* top_bnd,
  				double* curvilinear);

  void transFiniteInterpolation3D(int mx, int my, int mz,
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
  				  );

  int getBlock(const tarch::la::Vector<DIMENSIONS,double>& center,
	       const tarch::la::Vector<DIMENSIONS,double>& dx);

   private:

  const int num_nodes;
  const double dx;
  const double fault_position;

  double* lagrange_basis_at_nodes;
  double* denominator_lagrange;
  double* unif_mesh;

  Boundary_single_coordinate* boundary_x[2];
  Boundary_single_coordinate* boundary_y[2];
  Boundary_single_coordinate* boundary_z[2];

  double a_x;
  double a_y;
  double a_z;
  double b_x;
  double b_y;
  double b_z;

  double fault(double y, double z,double a_y, double b_y, double a_z, double b_z);
  double topography(double x, double z,double a_x, double b_x,double a_z, double b_z, double depth);
  
/* #if defined(USE_ASAGI) */
/*   double topography_fromASAGI(double x, double z, double* topography, easi::ArraysAdapter& adapter, easi::Component* model); */
/* #endif   */

  void getInterpolatedFace_fromBottomAndTop( int n_block, int n_left_right, int n_top_bottom, int face,
					   double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					   double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					   double* face_x, double* face_y , double* face_z);

  void interpolate(double x, double y, double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double& result);
  
  void interpolate3D(int x, int y, int z, double* dest_mesh, int num_nodes,double& result);

  double lagrangeBasis(double x,int i,int num_points);

  void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results);

  void getValuesAtQuadNodes3D(double* dest_mesh, int num_nodes, double* results);

  void computeDerivatives_x(int i, int j, double* values , int num_nodes, double& der_x, double dx);

  void computeDerivatives_y(int i, int j, double* values , int num_nodes, double& der_y, double dy);

  void computeDerivatives_x_3D(int i, int j, int k, double* values , int num_nodes, double& der_x, double dx);
  void computeDerivatives_y_3D(int i, int j, int k, double* values , int num_nodes, double& der_y, double dy);
  void computeDerivatives_z_3D(int i, int j, int k ,double* values , int num_nodes, double& der_z, double dz);


  void metricDerivativesAndJacobian(int num_nodes, double* curvilinear_x, double* curvilinear_y,double* gl_vals_x,double* gl_vals_y,double* q_x,double* q_y,double* r_x,double* r_y,double* jacobian, double dx, double dy);


  void metricDerivativesAndJacobian3D(int num_nodes,
				  double* curvilinear_x, double* curvilinear_y, double* curvilinear_z,
				  double* gl_vals_x, double* gl_vals_y, double* gl_vals_z,
				  double* q_x, double* q_y, double* q_z,
				  double* r_x, double* r_y, double* r_z,
				  double* s_x, double* s_y, double* s_z,
				  double* jacobian,
				  double dx, double dy, double dz
				    );


  double interpolate_fault_surface(double* X1, double* Y1, double* Z1, double y, double z, int my, int mz);
  double interpol2d_dG(int my, int mz, int nmy, int npy, int nmz, int npz, double y, double z, double* yj, double* zk, double* f);
  double lagrange(int m, int p, int i, double x, double *xi);

};
#endif
