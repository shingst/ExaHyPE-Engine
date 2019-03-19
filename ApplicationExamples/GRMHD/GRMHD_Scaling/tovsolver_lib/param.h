#ifndef TOVOSOLVER_STANDALONE_PARAMETERS
#define TOVOSOLVER_STANDALONE_PARAMETERS

#include <string>

namespace TOV {

  struct Parameters {
    int TOV_Num_TOVs = 1; // "The number of TOVs"
    int TOV_Solve_for_TOVs = 3; // "Solve for TOVs even if no TOV initial data was requested?"
    // 0:3 :: "depreciated in favour of TOVSolver::TOV_Enforce_Interpolation"

    bool TOV_Enforce_Interpolation = false; // "Enforce the interpolation of the data onto the Hydro GFs even without tov as specified initial data" STEERABLE=always

    int TOV_Num_Radial = 100000; // "The number of radial points for the ODE integration"
    double TOV_Rho_Central[10]; // "The central density"
    // default: 1e-3

    // the default atmosphere value
    double atmo_rho =  1.e-13; 

    // the default atmosphere value of pressure
    double atmo_press = 1.e-7;

    double TOV_Gamma = 2.0; // "The polytropic constant in P = K rho^Gamma"
    // 1.0: :: "The physical range at high Lorentz factors is [1,2], but otherwise higher values of gamma can also be used"

    double TOV_K = 100.0; // "The polytropic constant in P = K rho^Gamma"

    double TOV_dr[10]; // "The spacing in the radial direction on the 1d grid"
    // default: 5.e-4

    bool Perturb[10]; // "Add density perturbation (you should solve the IVP if true)"
    // default: no

    double Pert_Amplitude[10]; // "Amplitude of perturbation"
    // default: 0


    bool Perturb_Pressure[10]; // "Add pressure perturbation (you should solve the IVP if true)"
    // default: no

    double Pert_Press_Amplitude[10]; // "Amplitude of pressure perturbation"
    // default: 0


    double TOV_Position_x[10]; // "Position of neutron star, x coordinate" STEERABLE=always
    // default: 0

    double TOV_Position_y[10]; // "Position of neutron star, y coordinate" STEERABLE=always
    // 0.0

    double TOV_Position_z[10]; // "Position of neutron star, z coordinate" STEERABLE=always
    // 0.0

    // contravariant fluid three velocity as measured by the Eulerian observers: v^i = u^i / (alpha u^t) + beta^i / alpha
    // as used by HydroBase. Follows the Valencia formulation eg. eqs. 26 and 27 of Font et al's paper (gr-qc/9811015) or 
    // below Equ. 31 in http://relativity.livingreviews.org/Articles/lrr-2008-7/articlesu1.html#x6-30002.1

    double TOV_Velocity_x[10]; // "(fixed) Velocity of neutron star, x coordinate (caution!)" STEERABLE=always
    // 0.0

    double TOV_Velocity_y[10]; // "(fixed) Velocity of neutron star, y coordinate (caution!)" STEERABLE=always
    // 0.0

    double TOV_Velocity_z[10]; // "(fixed) Velocity of neutron star, z coordinate (caution!)" STEERABLE=always
    // 0.0


    bool TOV_ProperPosition = false; // "For use only with two NSs, atm only handles equal mass"

    bool TOV_Fast_Interpolation = true; // "Use faster interpolation algorithm? Default is yes."

    bool TOV_Clear_Initial_Data = true; // "Clear initial data (spacetime)? Default is yes."

    bool TOV_Use_Old_Initial_Data = false; // "Take old initial data into account (spacetime)? Default is no."

    bool TOV_Use_Old_Matter_Initial_Data = false; // "Use also old matter initial data? Default is no."

    bool TOV_Conformal_Flat_Three_Metric = false; // "Use conformal factor to get the 3-metric flat. default is no"

    std::string TOV_Combine_Method = "average"; // "Which combine method should be used."
    //"maximum" :: "Take the maximum of rho and gxx as clue for the rest as clue."
    // "average" :: "Take the average of all available parts."

    int TOV_Populate_Timelevels = 1; // "Populate that amount of timelevels" STEERABLE=always
    // :: "1 (default) to 3"

    int TOV_Momentum_Psi_Power = 0; // "Power of Psi to be multiplied with J^i for Mom"

    int TOV_fake_evolution = 0; // "Fake evolution by setting ID at every step" STEERABLE=always
    // "anything, 0 as off (default), everything else as on"

    std::string TOV_save_to_datafile = ""; // "Only save data to file and exit"

    inline Parameters() {
      for(int i=0; i<10; i++) {
        TOV_Rho_Central[i] = 1e-3;
        TOV_dr[i] = 5.e-4;
        Perturb[i] = false;
        Pert_Amplitude[i] = 0;
        Perturb_Pressure[i] = false;
        Pert_Press_Amplitude[i] = 0;
        TOV_Position_x[i] = 0;
        TOV_Position_y[i] = 0;
        TOV_Position_z[i] = 0;
        TOV_Velocity_x[i] = 0;
        TOV_Velocity_y[i] = 0;
        TOV_Velocity_z[i] = 0;
      }
    }
  }; 

} // end of namespace

#endif
