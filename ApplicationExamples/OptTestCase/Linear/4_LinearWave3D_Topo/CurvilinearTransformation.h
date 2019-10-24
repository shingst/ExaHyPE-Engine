#include "kernels/KernelUtils.h"
#include "kernels/GaussLegendreBasis.h"



//namespace Linear{
//  class CurvilinearTransformation;
//}


//class Linear::CurvilinearTransformation {
 // private:
  
// public:
//  CurvilinearTransformation();

double fault(double y, double z,double a_y, double b_y, double a_z, double b_z);


void getBoundaryCurves(int num_points,double offset_x, double offset_y,double width_x, double width_y , double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y);


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
						    double a_x, double a_y , double a_z,
						    double b_x, double b_y , double b_z,
						    double width_x, double width_y , double width_z,
						    double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
						    double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
						    double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
						    double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
						    double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
						    double* back_bnd_x, double* back_bnd_y, double* back_bnd_z);



//void transFiniteInterpolation( int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y );

void transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y );

void transFiniteInterpolation_singleCoordinate(int mx, int my, double* left_bnd, double* right_bnd, double* bottom_bnd, double* top_bnd, double* curvilinear);

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


void getInterpolatedFace_fromBottomAndTop( int n_block, int n_left_right, int n_top_bottom, int face,
					   double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					   double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					   double* face_x, double* face_y , double* face_z);

void interpolate(double x, double y, double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double& result);

void interpolate3D(double x, double y, double z, double* orig_mesh_x , double* orig_mesh_y, double* orig_mesh_z, double* dest_mesh, int num_nodes,double& result);

double lagrangeBasis(double x,double* points,int i,int num_points);

void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results);

void getValuesAtQuadNodes3D(double* orig_mesh_x , double* orig_mesh_y, double* orig_mesh_z, double* dest_mesh, int num_nodes, double* results);

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

