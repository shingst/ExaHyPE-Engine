/* file    tov.c
 * author  Frank Loeffler, converted from fortran thorn by Ian Hawke
 * date    2002/10/21
 * desc    TOV initial data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <algorithm>
#include <string>

#include <cstdio>
#include <iostream>

#include "tov.h"

using namespace std;
using namespace TOV;

#define NUMVARS 6

bool CCTK_EQUALS(const char* a, const char* b) {
	std::string A(a), B(b);
	std::transform(A.begin(), A.end(), A.begin(), ::tolower);
	std::transform(B.begin(), B.end(), B.begin(), ::tolower);
	return A==B;
}

void CCTK_WARN(int shall_exit, const char* message) {
	puts(message);
	if(shall_exit) std::abort();
}

#define CCTK_THORNSTRING "whatever"
#define CCTK_INFO(m) printf(m "\n")
#define CCTK_VInfo(something, msg...) printf(msg)

#define CONSTANT_c_SI  2.99792458e8
#define CONSTANT_c_cgi 2.99792458e10
#define CONSTANT_c_C   1.0e0

#define CONSTANT_G_SI  6.6732e-11
#define CONSTANT_G_cgi 6.6732e-8
#define CONSTANT_G_C   1.0e0

#define CONSTANT_Msolar_SI  1.987e30
#define CONSTANT_Msolar_cgi 1.987e33
#define CONSTANT_Msolar_C   1.0d0



/*@@
   @routine    TOV_Source_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler - converted fortran routine by Ian Hawke
   @desc
      The source terms for the ODEs. These are equations (2), (3), (4)
      and (18) from the Baumgarte notes.
      That is the vector in order is (P, m, phi, rbar).
@@*/
void TOV_C_Source_RHS(double r, double K, double Gamma,
                      double old_data[NUMVARS], double source_data[NUMVARS])
{
  double LOCAL_TINY, PI;
  double press, rho, eps, mu, m;
  double r_minus_two_m;

  LOCAL_TINY = 1.0e-35;
  PI=4.0*atan(1.0);

  press           = old_data[0];
  if (press < LOCAL_TINY)
    press = LOCAL_TINY;
  m               = old_data[1];

  rho = pow(press / K, 1.0 / Gamma);
  eps = press / (Gamma - 1.0) / rho;
  mu  = rho * (1.0 + eps);

  r_minus_two_m = r - 2.0 * m;

  if ((r<=0.0) && (m<=0.0))
  {
    source_data[1] = 0.0;
    source_data[2] = 0.0;
    source_data[3] = 0.0;
    source_data[4] = 0.0;
    source_data[5] = 0.0;
  }
  else
  {
    source_data[2] = (m + 4*PI * r*r*r * press) / r_minus_two_m / r;
    /* source_data[0] = -(press + mu) * source_data[2]; */
    source_data[0] = -(press + mu) *
                     (m + 4*PI * r*r*r * press) / r_minus_two_m / r;
    source_data[1] = 4*PI * r*r * mu;
    source_data[3] = (sqrt(r) - sqrt(r_minus_two_m)) / r / sqrt(r_minus_two_m);
    source_data[5] = 1.0/sqrt(1.0-2.0*m/r);
    source_data[4] = source_data[5] * 4*PI * rho * r*r;
  }
}

/*@@
 * routine sum_helper
 * desc
 *      small helper function to make TOV_Integrate_RHS more readable
@@*/
double sum_helper(double* vec, int max_ind, double factor){
    if(max_ind <= 0)
        return 0.0;
    double sum;
    for(int i=0; i < max_ind; i++)
        sum += vec[i]*factor;
    return sum;
}

/*@@
   @routine    TOV_Integrate_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
      Integrates the ODEs using RK4.
      We rescale at the end to match to a Schwarzschild exterior.
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV::TOVSolver::TOV_C_Integrate_RHS()
{
  double LOCAL_TINY;

  int star, star_i, i, TOV_Surface_Index;
  double old_data[NUMVARS], source_data[NUMVARS],
            in_data[NUMVARS], new_data[NUMVARS],
            k1[NUMVARS], k2[NUMVARS], k3[NUMVARS], k4[NUMVARS];
  double Surface_Mass, factor, local_rho;

  LOCAL_TINY = 1.0e-20;

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);
  assert(TOV_RProp_Surface!=0);

  double TOV_m_1d_surf[TOV_Num_TOVs];
  double TOV_mbary_1d_surf[TOV_Num_TOVs];
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);

  /* do it for all stars */
  for (star=0; star < TOV_Num_TOVs; star++)
  {
    /* remember array index */
    star_i = star * TOV_Num_Radial;
    const double rho_central=TOV_Rho_Central[star];

    /* clear arrays first */
    TOV_C_fill(&(TOV_press_1d[star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_phi_1d  [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_rbar_1d [star_i]), TOV_Num_Radial, 0.0);

    /* set start values */
    TOV_press_1d[star_i] = TOV_K *
                            pow(rho_central, TOV_Gamma);
    double TOV_m_1d = 0.0;
    double TOV_mbary_1d = 0.0;
    double TOV_rprop_1d = 0.0;
    TOV_Surface[star] = -1.0;
    TOV_Surface_Index = -1.0;

#define RKLOOP for (int rk=0; rk<NUMVARS; rk++)
    /* loop over all radii */
    for (i=star_i; (i < star_i+TOV_Num_Radial-1) &&
                            (TOV_Surface[star] < 0.0); i++)
    {
      double TOV_r_1d = i*TOV_dr[star] + sum_helper(TOV_dr, star-1, TOV_Num_Radial);
      /* set up RK arrays */
      old_data[0] = TOV_press_1d[i];
      old_data[1] = TOV_m_1d;
      old_data[2] = TOV_phi_1d[i];
      if (fabs(TOV_rbar_1d[i] - TOV_r_1d) < LOCAL_TINY)
        old_data[3] = 0.0;
      else
        old_data[3] = log(TOV_rbar_1d[i] / TOV_r_1d);
      old_data[4] = TOV_mbary_1d;
      old_data[5] = TOV_rprop_1d;

      /* usual RK4 */
      RKLOOP in_data[rk] = old_data[rk];
      TOV_C_fill(source_data, 6, 0.0);

      TOV_C_Source_RHS(TOV_r_1d,
                     TOV_K, TOV_Gamma,
                     in_data, source_data);

      RKLOOP k1[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + 0.5 * k1[rk];
      TOV_C_Source_RHS(TOV_r_1d+ 0.5 * TOV_dr[star],
                       TOV_K, TOV_Gamma,
                       in_data, source_data);

      RKLOOP k2[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + 0.5 * k2[rk];
      TOV_C_Source_RHS(TOV_r_1d+ 0.5 * TOV_dr[star],
                       TOV_K, TOV_Gamma,
                       in_data, source_data);

      RKLOOP k3[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + k3[rk];
      TOV_C_Source_RHS(TOV_r_1d+ TOV_dr[star],
                       TOV_K, TOV_Gamma,
                       in_data, source_data);
      RKLOOP k4[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP new_data[rk] = old_data[rk] + (k1[rk] + k4[rk] + 2.0 * (k2[rk] + k3[rk])) /6.0;

      TOV_press_1d[i+1] = new_data[0];
      TOV_phi_1d  [i+1] = new_data[2];
      TOV_rbar_1d [i+1] = (TOV_r_1d + TOV_dr[star]) * exp(new_data[3]);

      /* otherwise the code crashes later */
      if (TOV_press_1d[i+1] < 0.0)
          TOV_press_1d[i+1] = 0.0;

      local_rho = pow(TOV_press_1d[i+1] / TOV_K, 1.0 / TOV_Gamma);

      /* scan for the surface */
      if ( (local_rho <= 0.0) ||
           (TOV_press_1d[i+1] <= 0.0) )
      {
        TOV_Surface[star]   = TOV_r_1d;
        TOV_R_Surface[star] = TOV_rbar_1d[i];
        TOV_RProp_Surface[star] = TOV_rprop_1d;
        TOV_Surface_Index = i;
        TOV_m_1d_surf[star] =TOV_m_1d;
        TOV_mbary_1d_surf[star] =TOV_m_1d;
      }
      //update
      TOV_m_1d          = new_data[1];
      TOV_mbary_1d      = new_data[4];
      TOV_rprop_1d      = new_data[5];
    }
    if (TOV_Surface[star] < 0.0)
      CCTK_WARN(0, "Did not integrate out to surface of the star! "
                   "Increase TOV_dr or TOV_Num_Radial and rerun");

    Surface_Mass = TOV_m_1d_surf[star];
    factor = 0.5 * (sqrt(TOV_Surface[star] *
                         (TOV_Surface[star] - 2.00 * Surface_Mass)) +
                    TOV_Surface[star] - Surface_Mass) /
                    TOV_rbar_1d[TOV_Surface_Index];

    TOV_R_Surface[star] *= factor;
    for (i=star_i; i < star_i+TOV_Num_Radial; i++)
    {
      double TOV_r_1d = i*TOV_dr[star] + sum_helper(TOV_dr, star-1, TOV_Num_Radial);
      TOV_rbar_1d[i] *= factor;
      TOV_phi_1d[i]  -= TOV_phi_1d[TOV_Surface_Index] -
                        0.5 * log(1.0 - 2.0 * Surface_Mass / TOV_Surface[star]);
      /* match to Schwarzschield */
      if (i > TOV_Surface_Index)
      {
        TOV_press_1d[i] = 0.0;
        TOV_rbar_1d [i] = 0.5 *
                          (sqrt(TOV_r_1d*(TOV_r_1d - 2.0*Surface_Mass)) +
                          TOV_r_1d - Surface_Mass);
        //TOV_m_1d[i]     = Surface_Mass;
        TOV_phi_1d[i]   = 0.5 * log( 1.0 - 2.0 * Surface_Mass / TOV_r_1d);
        //TOV_mbary_1d[i] = TOV_mbary_1d[TOV_Surface_Index];
      }
    }
  }
  /*CCTK_INFO("Integrated TOV equation");
   * TODO save also these values for output
  // do some info */
  /*CCTK_VInfo(CCTK_THORNSTRING, "Information about the TOVs used:");
  CCTK_VInfo("", "TOV    radius    mass  bary_mass mass(g) cent.rho rho(cgi)        K   K(cgi)    Gamma");
  for (i=0; i<TOV_Num_TOVs; i++)
    if (fabs(TOV_Gamma - 2.0) < LOCAL_TINY)
      CCTK_VInfo("","  %d  %8g %8g %8g %8.3g %8g %8.3g %8g %8.3g %8g",
                 (int)i+1, TOV_R_Surface[i],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_mbary_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1]*CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i]/pow(CONSTANT_G_cgi,3.0)/
                                    pow(CONSTANT_Msolar_cgi,2.0)*
                                    pow(CONSTANT_c_cgi,6.0),
                 TOV_K,
                 TOV_K*pow(CONSTANT_G_cgi,3.0)*
                          pow(CONSTANT_Msolar_cgi,2.0)/
                          pow(CONSTANT_c_cgi,4.0),
                 TOV_Gamma);
    else
      CCTK_VInfo("","  %d  %8g %8g %8.3g %8g %8.3g %8g %8g",
                 (int)i+1, TOV_R_Surface[i],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1]*CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i]/pow(CONSTANT_G_cgi,3.0)/
                                    pow(CONSTANT_Msolar_cgi,2.0)*
                                    pow(CONSTANT_c_cgi,6.0),
                 TOV_K, TOV_Gamma);*/

}

/*----------------------------------------------------------------------------*/

/* utility routine
 * recursive search-routine for arrays
 * here used to look for the last index in an ordered array with its
 * value < goal
 */
int TOV_C_find_index(int   array_size,
                          double *array,
                          double  goal,
                          int   lower_index,
                          int   upper_index)
{
  int middle_index;

  if (lower_index >= (upper_index-1))
      return lower_index;

  middle_index = (lower_index + upper_index) /2;

  if (array[middle_index] < goal)
    return TOV_C_find_index(array_size, array, goal, middle_index, upper_index);
  else
    return TOV_C_find_index(array_size, array, goal, lower_index, middle_index);
}


/* utility rountine
 * interpolates from (thorn-internal) 1D-data to Cactus 3D-grid */
/* input is all but *press_point *phi_point and *r_point */
void TOV::TOVSolver::TOV_C_interp_tov_isotropic(
                                int  star,
                                double *TOV_press_1d_local,
                                double *TOV_phi_1d_local,
                                double *TOV_rbar_1d_local,
                                double *r,
                                double surface,
                                double *press_point,
                                double *phi_point,
                                double *r_point)
{
  int  left_index;
  double h, M;

  if (*r < 0.0)
    CCTK_WARN(0, "Negative radius found");
  if (*r < TOV_rbar_1d_local[1])
    *r=TOV_rbar_1d_local[1];
  if (*r > TOV_rbar_1d_local[TOV_Num_Radial-2])
  {
      double TOV_r_1d = (TOV_Num_Radial-1)*TOV_dr[star] + sum_helper(TOV_dr, star-1, TOV_Num_Radial);
      *press_point= 0.0;
      M = 0.5 * TOV_r_1d *
        (1.0 - exp(2.0*TOV_phi_1d_local[TOV_Num_Radial-1]));
      *r_point=(2* *r+M)*(2* *r+M)*0.25/ *r;
      *phi_point=0.5*log(1-2*M/ *r_point);
      return;
  }

  if (TOV_Fast_Interpolation)
    left_index = TOV_C_find_index(TOV_Num_Radial-1, TOV_rbar_1d_local, *r, 0,
        TOV_Num_Radial-1);
  else
  {
    left_index=0;
    while( (left_index < TOV_Num_Radial-2) &&
        (TOV_rbar_1d_local[left_index+1] < *r) )
      left_index++;
  }

  h = (*r - TOV_rbar_1d_local[left_index]) / (TOV_rbar_1d_local[left_index+1] - TOV_rbar_1d_local[left_index]);

  *phi_point   = (1.0 - h) * TOV_phi_1d_local[left_index] + h  * TOV_phi_1d_local[left_index+1];

  double TOV_r_1d = (left_index)*TOV_dr[star] + sum_helper(TOV_dr, star-1, TOV_Num_Radial);
  *r_point     = (1.0 - h) * TOV_r_1d + h  * (TOV_r_1d + TOV_dr[star]);
  if (*r_point < surface)
    *press_point = (1.0 - h) * TOV_press_1d_local[left_index] + h  * TOV_press_1d_local[left_index+1];
  else
    *press_point = 0.0;
}

/*@@
   @routine    TOV_Exact
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
       Schedule routine for interpolation of 1D to 3D grid
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV::TOVSolver::Interpolate(const double pos[3], idvars& id)
{
  double *press_point, *rho_point, *eps_point,
         *mu_point, *phi_point, *r_point;

  int star;
  double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  double *r_to_star;
  double g_diag, max_g_diag, max_rho;
  double my_velx, my_vely, my_velz, my_psi4;
  double PI, local_tiny;

  PI=4.0*atan(1.0);
  local_tiny=1.0e-14;

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);

  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);

  /* allocate local arrays */ 
  r_to_star   = new double[TOV_Num_TOVs];
  press_point = new double[TOV_Num_TOVs];
  rho_point   = new double[TOV_Num_TOVs];
  eps_point   = new double[TOV_Num_TOVs];
  mu_point    = new double[TOV_Num_TOVs];
  phi_point   = new double[TOV_Num_TOVs];
  r_point     = new double[TOV_Num_TOVs];

  /* use the fast interpolation? only useful for testing this */
  if (TOV_Fast_Interpolation == 0) {
    CCTK_INFO("Interpolating the slow way.");
  }


  /* remember the old conformal factor to the power of 4 */
  my_psi4=pow(id.psi, 4.0);

  for (star=0; star<TOV_Num_TOVs; star++)
  {
    r_to_star[star] =
      sqrt( (pos[0]-TOV_Position_x[star]) *
          (pos[0]-TOV_Position_x[star]) +
          (pos[1]-TOV_Position_y[star]) *
          (pos[1]-TOV_Position_y[star]) +
          (pos[2]-TOV_Position_z[star]) *
          (pos[2]-TOV_Position_z[star]) );
    int star_i = star * TOV_Num_Radial;

    /* do the actual interpolation */
    TOV_C_interp_tov_isotropic(star,
        &(TOV_press_1d[star_i]), &(TOV_phi_1d[star_i]),
        &(TOV_rbar_1d[star_i]),
        &(r_to_star[star]), TOV_Surface[star],
        &(press_point[star]),
        &(phi_point[star]), &(r_point[star]));

    if (Perturb_Pressure[star])
      press_point[star] *= 1 + Pert_Press_Amplitude[star];

    /* is some perturbation wanted? */
    if (Perturb[star] == 0)
      rho_point[star] = pow(press_point[star]/TOV_K,
          1.0/TOV_Gamma);
    else
      rho_point[star] = pow(press_point[star]/TOV_K,
          1.0/TOV_Gamma) *
        (1.0 +
         Pert_Amplitude[star] *
         cos(PI/2.0 * r / TOV_R_Surface[star]));

    if (rho_point[star] > local_tiny)
      eps_point[star] = press_point[star] / (TOV_Gamma - 1.0)
        /  rho_point[star];
    else
      eps_point[star] = 0.0;
    mu_point[star]  = rho_point[star] * (1.0 + eps_point[star]);
  }
    /* to do this, we use here simply the max of the gxx-value */
    star=0;
    max_g_diag = 0.0;
    max_rho = rho_point[0];
    for (int star_i=0; star_i<TOV_Num_TOVs; star_i++)
    {
      g_diag = (r_point[star_i] / (r_to_star[star_i] + 1.0e-30)) *
        (r_point[star_i] / (r_to_star[star_i] + 1.0e-30));
      if ((g_diag - max_g_diag) > local_tiny)
      {
        max_g_diag=g_diag;
        star=star_i;
      }
      if ((rho_point[star_i] - max_rho) > local_tiny)
      {
        max_rho=rho_point[star_i];
        star=star_i;
      }
    }
    /* no psi, since it is 1.0 here */
    /* but maybe we want to have it != 1.0 */
    if (TOV_Conformal_Flat_Three_Metric)
    {
      id.psi = pow(max_g_diag, 0.25);
      my_psi4 = max_g_diag;
      id.gam[0][0] = id.gam[1][1] = id.gam[2][2] = 1.0;
      id.gam[0][1] = id.gam[0][2] = id.gam[1][2] = 0.0;
    }
    else
    {
      id.gam[0][0] = max_g_diag;
      id.gam[1][1] = max_g_diag;
      id.gam[2][2] = max_g_diag;
      id.gam[0][1] = id.gam[0][2] = id.gam[1][2] = 0.0;
    }
    id.alp = exp(phi_point[star]);
    id.beta[0] = 0.0;
    id.beta[1] = 0.0;
    id.beta[2] = 0.0;
    my_velx=TOV_Velocity_x[star];
    my_vely=TOV_Velocity_y[star];
    my_velz=TOV_Velocity_z[star];

    /* set to defined velocity. default is 0.0 because other velocities
     * violate Einsteins equations */
    id.vel[0] = my_velx;
    id.vel[1] = my_vely;
    id.vel[2] = my_velz;

    id.w_lorentz = 1/sqrt(1.0-(
          id.gam[0][0] * id.vel[0] * id.vel[0]+
          id.gam[1][1] * id.vel[1] * id.vel[1]+
          id.gam[2][2] * id.vel[2] * id.vel[2]+
          2*id.gam[0][1] * id.vel[0] * id.vel[1]+
          2*id.gam[0][2] * id.vel[0] * id.vel[2]+
          2*id.gam[1][2] * id.vel[1] * id.vel[2])*
        my_psi4);

    id.rho = max_rho;
    id.eps = eps_point[star];
    id.press = press_point[star];


  /* free local arrays */
  delete r_to_star;
  delete press_point;
  delete rho_point;
  delete eps_point;
  delete mu_point;
  delete phi_point;
  delete r_point;
} // end of function 

