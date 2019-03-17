#include "tov.h"

#include <cassert>
#include <sstream>


using namespace TOV;

void TOV::TOVSolver::TOV_C_Allocate()
{
  TOV_Surface = new double[TOV_Num_TOVs];
  TOV_R_Surface = new double[TOV_Num_TOVs];
  TOV_RProp_Surface = new double[TOV_Num_TOVs];

  //TODO find better way than storing these as
  //TOV_Num_Radial = 40000000
  TOV_rbar_1d = new double[TOV_Num_Radial * TOV_Num_TOVs];
  TOV_press_1d = new double[TOV_Num_Radial * TOV_Num_TOVs];
  TOV_phi_1d = new double[TOV_Num_Radial * TOV_Num_TOVs];

  //surface helper variables
  double factor;
  int TOV_Surface_Index;
  double Surface_Mass;
}

TOV::TOVSolver::~TOVSolver ()
{
  delete TOV_Surface;
  delete TOV_R_Surface;
  delete TOV_RProp_Surface;

  delete TOV_rbar_1d;
  delete TOV_press_1d;
  delete TOV_phi_1d;
}

/* - utility routine
   - fills an real-array 'var' of size 'i' with value 'r' */
void TOV::TOV_C_fill(double *var, int i, double r)
{
  for (i-- ;i >= 0; i--)
    var[i]=r;
}

/*void TOV::TOV_Copy(int size, double *var_p, double *var)
{
    for(int i=0; i<size; i++)
        var_p[i] = var[i];
}*/

#define SS(x) (static_cast<std::ostringstream&>(( std::ostringstream() << x )).str())
#define INDENT "    "
#define NL SS(", " << std::endl << INDENT)
#define DESC(var) SS(#var << " = " << var)
#define MEANS(var, num, text) SS((var == num ? SS(DESC(var) << " = " text) : ""))
#define VEC3(var) SS(#var << " = " << "[" << var[0] << ","	<< var[1] << "," << var[2] << "]")
#define VEC4(var) SS(#var << " = " << "[" << var[0] << ","	<< var[1] << "," << var[2] << "," << var[3] << "]")
#define MAT(var) SS(#var << " = " << "[" << NL << VEC3(var[0]) << NL << VEC3(var[1]) << NL << VEC3(var[2]) << NL <<  "]")

std::string TOV::idvars::toString() const {
	std::stringstream out;
	out << "idvars(" << NL;
	
	out << DESC(rho) << NL;
	out << VEC3(vel) << NL;
	out << DESC(press) << NL;
	out << DESC(eps) << NL;

	out << DESC(psi) << NL;
	out << DESC(alp) << NL;
	out << VEC3(beta) << NL;
	out << MAT(gam) << NL;
	
	out << ")"; // closing bracket of InitialData()

	return out.str();
}

