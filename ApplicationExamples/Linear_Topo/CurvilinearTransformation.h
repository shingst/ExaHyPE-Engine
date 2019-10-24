#include "kernels/KernelUtils.h"
#include "kernels/GaussLegendreBasis.h"



//namespace Linear{
//  class CurvilinearTransformation;
//}


//class Linear::CurvilinearTransformation {
 // private:
  
// public:
//  CurvilinearTransformation();

void getBoundaryCurves(int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y);


void transFiniteInterpolation( int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y );

void interpolate(double x, double y, double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double& result);

double lagrangeBasis(double x,double* points,int i,int num_points);

void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results);

void computeDerivatives_x(int i, int j, double* values , int num_nodes, double& der_x);

void computeDerivatives_y(int i, int j, double* values , int num_nodes, double& der_y);

void metricDerivativesAndJacobian(int num_nodes, double* curvilinear_x, double* curvilinear_y,double* gl_vals_x,double* gl_vals_y,double* q_x,double* q_y,double* r_x,double* r_y,double* jacobian);

//  void test();

//};
