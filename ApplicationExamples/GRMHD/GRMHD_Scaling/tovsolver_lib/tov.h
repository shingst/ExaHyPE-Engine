#ifndef __TOVSOLVER_STANDALONE__
#define __TOVSOLVER_STANDALONE__

#include "param.h"
#include <cmath>

namespace TOV {

  struct idvars {
    double rho;
    double vel[3];
    double press;
    double w_lorentz;

    double eps;

    double gam[3][3], alp, beta[3];

    double psi; // only set if some conformal treatment is turned on

    idvars() {
      rho = NAN;
      for(int i=0; i<3; i++) vel[i] = NAN;
      press = NAN;
      w_lorentz = NAN;
      eps = NAN;
      for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        beta[i] = NAN;
        gam[i][j] = NAN;
      }
      alp = NAN;
      psi = NAN;
    }

    std::string toString() const;
  };


  struct TOVSolver : public Parameters {
    TOVSolver() {}
    virtual ~TOVSolver();

    double * TOV_Surface=0;
    double * TOV_R_Surface=0;
    double * TOV_RProp_Surface=0;

    double * TOV_rbar_1d=0;
    double * TOV_press_1d=0;
    double * TOV_phi_1d=0;

    void TOV_C_Allocate();

    //surface helper variables
    double factor;
    int TOV_Surface_Index;
    double Surface_Mass;

    // Main Function to compute TOV solution
    void TOV_C_Integrate_RHS();

    void Setup() {
        TOV_C_Allocate();
        TOV_C_Integrate_RHS(); 
    }


    void TOV_C_interp_tov_isotropic(
        int  star,
        double *TOV_press_1d_local,
        double *TOV_phi_1d_local,
        double *TOV_rbar_1d_local,
        double *r,
        double surface,
        double *press_point,
        double *phi_point,
        double *r_point);


    void Interpolate(const double pos[3], idvars& id);
  };

  // Functions which do not depend on the solution or parameters
  void TOV_C_fill(double *var, int i, double r);
  //void TOV_Copy(int size, double *var_p, double *var);


} // end of namespace

#endif
