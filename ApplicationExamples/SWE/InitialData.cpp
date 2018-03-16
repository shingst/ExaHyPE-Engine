#include "InitialData.h"
#include "MySWESolver_Variables.h"
#include "MySWESolver.h"

#include <cmath>

using namespace std;

const double grav=9.81;
///// 2D /////

#ifdef Dim2

/*
 * Constant water height with "exponential" hump in bathymetry.  
 */
void SWE::ExpBreakProblem(const double* const x,double* Q) {
  MySWESolver::Variables vars(Q);

  if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
    vars.h() = 4.0;
    vars.hu()= 0.0;
    vars.hv()= 0.0;
    vars.b() = exp(-pow((x[0]-5),2)-pow((x[1]-5),2));
  } else {
    vars.h() = 4.0;
    vars.hu()= 0.0;
    vars.hv()= 0.0;
    vars.b() = 0.0;
  }
}

/*
 * Constant water height with discontinuous jump in bathymetry.  
 */
void SWE::DamBreakProblem(const double* const x,double* Q) {
  MySWESolver::Variables vars(Q);

  if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
    vars.h() = 4.0;
    vars.hu()= 0.0;
    vars.hv()= 0.0;
    vars.b() = 1;
  } else {
    vars.h() = 4.0;
    vars.hu()= 0.0;
    vars.hv()= 0.0;
    vars.b() = 0.0;
  }
}

/*
 * Sea at rest (steady state).  
 */
void SWE::SeaAtRestProblem(const double* const x,double* Q) {
  MySWESolver::Variables vars(Q);

  if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
    vars.h() = 3.0;
    vars.hu()= 0.0;
    vars.hv()= 0.0;
    vars.b() = 1;
  } else {
    vars.h() = 4.0;
    vars.hu()= 0.0;
    vars.hv()= 0.0;
    vars.b() = 0.0;
  }
}

/*
* Simulates a channel with a linearly rising ramp at the left end of the channel.
* To be used with WALL boundaries.
*/
void SWE::SteadyRunUpLinear(const double* const x, double* Q) {
   MySWESolver::Variables vars(Q);

   const double d =1;
   const double xr=19.85;

   vars.b()  = (x[0]>xr) ? 0 : (-d/xr)*x[0] +d;
   vars.h()  = (x[0]>xr) ? d : (d/xr)*x[0]; 
   vars.hu() = 0.0;
   vars.hv() = 0.0;

}

/*
* Simulates a wave which approaches a linearly rising ramp at the left end of the channel.
* To be used with WALL boundaries.
*/
void SWE::RunUpLinear(const double* const x, double* Q) {
   MySWESolver::Variables vars(Q);

   const double as=0.3;
   const double d =1;
   const double xr=19.85;
   const double alpha=atan(d/xr);
   const double xs= d/tan(alpha) + sqrt(4.0/(3.0*as))*acosh(sqrt(20.0));

   vars.b()  = (x[0]>xr) ? 0 : (-d/xr)*x[0] +d;
   vars.h()  = (x[0]>xr) ? d : (d/xr)*x[0]; 
   vars.h()  += as*(1.0/pow(cosh(sqrt(3.0*as/4.0)*(x[0]-xs)),2));
   vars.hu() = -vars.h()* as*(1.0/pow(cosh(sqrt(3*as/4)*(x[0]-xs)),2))*sqrt(grav/d);
   vars.hv() = 0.0;

}

/*
* Steady-state test case for artificial continental shelf.
*/
void SWE::SteadyRunUpShelf(const double* const x, double* Q) {
   MySWESolver::Variables vars(Q);

   const double d1 = 200.0;
   const double d2 = 2000.0;
   const double xb = 10.0;
   const double xe = 50.0;

   if(x[0]>xe)
   {
      vars.h()  = d2;
      vars.b()  = 0.0;
   }
   else if ( x[0]>xb)
   {
      vars.h()  = ((d1-d2)/(xb-xe))*(x[0]-xb) +d1;
      vars.b()  = d2 - vars.h();
   }
   else
   {
      vars.h()  = d1;
      vars.b()  = d2-d1;
   }

   vars.hu() = 0.0;
   vars.hv() = 0.0;
}

/*
* Run-up on artificial continental shelf.
*/
void SWE::RunUpShelf(const double* const x, double* Q) {
   MySWESolver::Variables vars(Q);

   const double d1 = 200.0;
   const double d2 = 2000.0;
   const double xb = 10.0;
   const double xe = 50.0;
   const double xw = 500.0;
   const double aw = 1.0;

   vars.hu()=0.0;
   vars.hv()=0.0;

   if(x[0]>xe)
   {
      vars.h()  = d2;
      vars.b()  = 0.0;
   }
   else if (x[0]>xb)
   {
      vars.h()  = ((d1-d2)/(xb-xe))*(x[0]-xb) +d1;
      vars.b()  = d2 - vars.h();
   }
   else
   {
      vars.h()  = d1;
      vars.b()  = d2-d1;
   }

   vars.h()  += aw*(1.0/pow(cosh(x[0]-xw),2));
   vars.hu() +=-vars.h()*aw*(1.0/pow(cosh(x[0]-xw),2))*sqrt(grav/d2);
   //vars.hu() = -vars.h()* aw*(1.0/pow(cosh(sqrt(x[0]-xw),2))*sqrt(grav/d);
   //vars.hu() = -vars.h()* sqrt(grav*vars.h());
}

#endif

void SWE::initialData(const double* const x,double* Q) {
 //DamBreakProblem(x,Q);
 //ExpBreakProblem(x,Q);
 //SteadyRunUpLinear(x,Q);
 //RunUpLinear(x,Q);
 //SteadyRunUpShelf(x,Q);
  //DamBreakProblem(x,Q);
  SeaAtRestProblem(x,Q);
}
