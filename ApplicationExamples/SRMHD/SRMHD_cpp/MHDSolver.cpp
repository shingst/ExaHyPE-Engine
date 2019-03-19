#include "MHDSolver.h"
//#include "fortran.h" _ltob

#include <memory>
#include <cstring>

#include <cmath>
#include <stdio.h>

#include "C2P-MHD.h"
#include "C2PRoutines.h"

#include "tarch/la/ScalarOperations.h"

// Fortran->C++:
// Regex Patterns:
// Semicolons:     (.$)                       ,
// Square:         ([a-z]+)\*\*2              ,      \1\*\1
// Vector indices: \(([0-9]+)\)               ,      \[\1\]     ,  don't forget to adjust indices; -1
// Variables       ^  ([a-z]+)                ,      double \1
// Col->Row        \(([0-9]+),\s*([0-9]+)\)   ,      [\2][\1]

// uncomment to enable
//#define USE_FORTRAN_HELPER_FUNCS

#define MHD_VARIABLES 9
#define MHD_PARAMETERS 0

#define DivCleaning_a 0.5 // 0 M. Dumbser, 0.5 Vasco
#define gamma 5.0

void MHDSolver::MHDSolver::init() {
}


void FUNC_C2P_RMHD1(const double x,double* const f,double* const df,const double gam,const double d,const double e,const double s2,const double b2,
    const double sb2,double* const w_out) {
  //
  // This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  // and it corresponds to their choice 3 in Section 3.2
  //
  double third = 1./3.;

  double v2=x;
  double rho=d*std::sqrt(1.-v2);

  double c3=1.-gam*(1.-v2);
  double c2=gam*rho+.5*b2*(1.+v2)-e;
  double c0=-0.5*sb2;

  //  For every x=v2, we solve a cubic in W of the form:
  //  c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  //  W=y of the paper. If sb=0 ( meaning c0 = 0),
  //  w = -c2/c3 > 0 and dw = 0 in the do loop below.
  //  If -c2/c3 < 0 when sb=0, which makes w=0,
  //  this is a signature that something was wrong before.
  double w;
  if ( abs ( c0) < 1.0e-20) {
    w = -c2 / c3;
  } else {
    w = std::max ( - c2 / c3, std::pow(( -c0 / c3), third));   // ( -c0 / c3)**third
    for(int iter = 0; iter < 200; iter++) {
      double dw = -((c3*w + c2)*w*w + c0)/((3*c3*w + 2*c2)*w);
      if (std::abs(dw/w)<1.e-10) {
        break; // happy breakdown
      }
      w = w + dw;
    }
  };
  *w_out = w;

  double dc3   = gam;
  double dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2));
  double dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2);
  double wb    = w + b2;
  double vb2 = sb2 / (w*w);
  *f   = wb*wb * v2 - ( 2.0 * w + b2) * vb2 - s2;
  *df  = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2));

//  std::cout << "w_out="<<*w_out;
//  std::cout << ",f="<<*f;
//  std::cout << ",df="<<*df << std::endl;
}

double RTSAFE_C2P_RMHD1(double X1,double X2,double XACC,double gam,double d,
    double e,double s2,double b2,double sb2,double* const w,bool* failed) {
  int MAXIT=200;
  *failed = false;
  double FL;
  double DF;

#ifdef USE_FORTRAN_HELPER_FUNCS
  func_c2p_rmhd1_(&X1,&FL,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
  FUNC_C2P_RMHD1(X1,&FL,&DF,gam,d,e,s2,b2,sb2,w);
#endif
  double FL_test,DF_test,w_test;
  FUNC_C2P_RMHD1(X1,&FL_test,&DF_test,gam,d,e,s2,b2,sb2,&w_test);
//  assertionNumericalEquals(FL,FL_test);
//  assertionNumericalEquals(DF,DF_test);
//  assertionNumericalEquals(*w,w_test);

  double RTSAFE_C2P_RMHD1_result = 0; // don't know if this is a good init value but flag failed should safe us
  if(FL==0) {
    RTSAFE_C2P_RMHD1_result=X1;
    return RTSAFE_C2P_RMHD1_result;
  }

  double FH;
#ifdef USE_FORTRAN_HELPER_FUNCS
  func_c2p_rmhd1_(&X2,&FH,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
  FUNC_C2P_RMHD1(X2,&FH,&DF,gam,d,e,s2,b2,sb2,w);
#endif
  double FH_test;
  FUNC_C2P_RMHD1(X2,&FH_test,&DF_test,gam,d,e,s2,b2,sb2,&w_test);
//  assertionNumericalEquals(FH,FH_test);
//  assertionNumericalEquals(DF,DF_test);
//  assertionNumericalEquals(*w,w_test);

  if(FH==0) {
     RTSAFE_C2P_RMHD1_result=X2;
     return RTSAFE_C2P_RMHD1_result;
  }
  if(FL*FH>0) {
     *failed = true;
     return RTSAFE_C2P_RMHD1_result;
  }

  double XL,XH;
  double SWAP;
  if(FL<0) {
    XL=X1;
    XH=X2;
  } else {
    XH=X1;
    XL=X2;

    SWAP=FL;
    FL=FH;
    FH=SWAP;
  }
  RTSAFE_C2P_RMHD1_result=.5*(X1+X2);
  double DXOLD=std::abs(X2-X1);
  double DX=DXOLD;

  double F;
#ifdef USE_FORTRAN_HELPER_FUNCS
  func_c2p_rmhd1_(&RTSAFE_C2P_RMHD1_result,&F,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
  FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1_result,&DF,&DF,gam,d,e,s2,b2,sb2,w);
#endif
  //double F_test;
  //FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1_result,&F_test,&DF_test,gam,d,e,s2,b2,sb2,&w_test);
//  assertionNumericalEquals(F,F_test);
//  assertionNumericalEquals(DF,DF_test);
//  assertionNumericalEquals(*w,w_test);
  for (int J=1; J<MAXIT; ++J) {
     if(((RTSAFE_C2P_RMHD1_result-XH)*DF-F)*((RTSAFE_C2P_RMHD1_result-XL)*DF-F)>=0
          || std::abs(2.*F)>std::abs(DXOLD*DF) ) {
        DXOLD=DX;
        DX=0.5*(XH-XL);
        RTSAFE_C2P_RMHD1_result=XL+DX;
        if(XL==RTSAFE_C2P_RMHD1_result) {
          return RTSAFE_C2P_RMHD1_result;
        }
     } else {
        DXOLD=DX;
        DX=F/DF;
        double TEMP=RTSAFE_C2P_RMHD1_result;
        RTSAFE_C2P_RMHD1_result=RTSAFE_C2P_RMHD1_result-DX;
        if (TEMP==RTSAFE_C2P_RMHD1_result) {
          return RTSAFE_C2P_RMHD1_result;
        }
     }
     if (std::abs(DX)<XACC) {
       return RTSAFE_C2P_RMHD1_result;
     }
#ifdef USE_FORTRAN_HELPER_FUNCS
     func_c2p_rmhd1_(&RTSAFE_C2P_RMHD1_result,&F,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
     FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1_result,&F,&DF,gam,d,e,s2,b2,sb2,w);
     //FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1_result,&F_test,&DF_test,gam,d,e,s2,b2,sb2,&w_test);
#endif
//     assertionNumericalEquals(F,F_test);
//     assertionNumericalEquals(DF,DF_test);
//     assertionNumericalEquals(*w,w_test);
     if(F<0) {
        XL=RTSAFE_C2P_RMHD1_result;
        FL=F;
     } else {
        XH=RTSAFE_C2P_RMHD1_result;
        FH=F;
     }
  }
  *failed = true;
  return RTSAFE_C2P_RMHD1_result;
}

void MHDSolver::MHDSolver::Cons2Prim(const double* const Q,double* const V) {
  double tol       = 1e-8;
//  double third     = 1.0/3.0; // unused for this Cons2Prim
  double p_floor   = 1.0e-5;
  double rho_floor = 1.0e-4;
  //;
  // First option [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)];
  bool failed   = false;
  double gamma1 = gamma/(gamma - 1.0);
  double gam    = 1.0/gamma1;
  double d      = Q[0];
  double sx     = Q[1];
  double sy     = Q[2];
  double sz     = Q[3];
  double e      = Q[4];
  double bx     = Q[5];
  double by     = Q[6];
  double bz     = Q[7];
  //
  double s2     = sx*sx + sy*sy + sz*sz;
  double b2     = bx*bx + by*by + bz*bz;
  double sb     = sx*bx + sy*by + sz*bz;
  double sb2    = sb*sb;
  double eps    = 1.e-10;
  //
  double x1     = 0.;
  double x2     = 1.-eps;

  double w;     // input var
  double v2     = RTSAFE_C2P_RMHD1(x1,x2,tol,gam,2,e,s2,b2,sb2,&w,&failed); // w is input parameter
  //
  if (failed) {
    double p    = p_floor;
    double rho  = rho_floor;
    double vx   = 0.0;
    double vy   = 0.0;
    double vz   = 0.0;
    bx   = 0.0;
    by   = 0.0;
    bz   = 0.0;

    V[0] = rho;
    V[1] = vx;
    V[2] = vy;
    V[3] = vz;
    V[4] = p;
    V[5] = bx;
    V[6] = by;
    V[7] = bz;
    V[8] = Q[8];
  } else {
    double den  = 1.0/(w+b2);
    double vb   = sb/w;
    //
    double rho  = d*sqrt(1.-v2);
    double vx   = (sx + vb*bx)*den;
    double vy   = (sy + vb*by)*den;
    double vz   = (sz + vb*bz)*den;
    double p    = std::max(1.e-15, gam*(w*(1.-v2)-rho));

    V[0] = rho;
    V[1] = vx;
    V[2] = vy;
    V[3] = vz;
    V[4] = p;
    V[5] = bx;
    V[6] = by;
    V[7] = bz;
    V[8] = Q[8];
  }
}


// V: prim
// Q: cons
void MHDSolver::MHDSolver::Prim2Cons(const double* const V,double* const Q) {
  double rho    = V[0];
  double vx     = V[1];
  double vy     = V[2];
  double vz     = V[3];
  double p      = V[4];
  double bx     = V[5];
  double by     = V[6];
  double bz     = V[7];

  double ex     = - (vy*bz - vz*by);
  double ey     = - (vz*bx - vx*bz);
  double ez     = - (vx*by - vy*bx);

  double v2     = vx*vx + vy*vy + vz*vz;
  double b2     = bx*bx + by*by + bz*bz;
  double e2     = ex*ex + ey*ey + ez*ez;

  if (v2 > 1.0) {
     std::cerr << "Superluminal velocity in PDEPrim2Cons!!" << std::endl;
     exit(1);
  }
  double lf     = 1.0 / std::sqrt(1.0 - v2);
  double gamma1 = gamma/(gamma-1.0);
  double w      = rho + gamma1*p;
  double ww     = w*lf*lf;
  double uem    = 0.5*(b2+e2);

  Q[0]   = rho*lf;
  Q[1]   = ww*vx + (ey*bz - ez*by);
  Q[2]   = ww*vy + (ez*bx - ex*bz);
  Q[3]   = ww*vz + (ex*by - ey*bx);
  Q[4]   = ww - p + uem;
  //
  Q[5]   = bx;
  Q[6]   = by;
  Q[7]   = bz;
  //
  Q[8]   = V[8];
}

void MHDSolver::MHDSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* const lambda) {
  // These are not the exact eigenvalues, instead of Lambda(1..9)
  // we compute only two eigenvalues: Approximate magnetosonics // cons2prim
  double V[MHD_VARIABLES];
//  int i;
//  pdecons2prim_(V,const_cast<double*>(Q),&i);
  Cons2Prim(Q,V);
  double rho = V[0];
  double vx  = V[1];
  double vy  = V[2];
  double vz  = V[3];
  double p   = V[4];
  double bx  = V[5];
  double by  = V[6];
  double bz  = V[7];

  double ex     = - (vy*bz - vz*by);
  double ey     = - (vz*bx - vx*bz);
  double ez     = - (vx*by - vy*bx);
  double gamma1 = gamma/(gamma-1.0);
  double b2     = bx*bx + by*by + bz*bz;
  // v2     = vx*vx + vy*vy + vz*vz;
  double e2     = ex*ex + ey*ey + ez*ez;
  double w      = rho + gamma1*p;
  // vn     = vx*nv[1] + vy*nv[2] + vz*nv[3];

  // the velocity vector is dimension agnostic;
  double vel[3] = { vx,vy,vz };
  double v2     = vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2];
  double vn     = vel[normalNonZeroIndex];

  double cs2    = gamma * p / w;
  double ca2    = (b2 - e2) / ( w + b2 - e2);
  double a2     = cs2 + ca2 - cs2 * ca2;
  double den    = 1.0/( 1.0 - v2 * a2);
  double vf1    = den * ( 1.0 - a2) * vn;
  double vf2    = den * std::sqrt(a2 * ( 1.0 - v2) * ((1.0 - v2 * a2) - (1.0 - a2) * vn*vn));

  lambda[0] = vf1 + vf2;
  lambda[1] = vf1 - vf2;
  lambda[2] = 0.0;
  lambda[3] = 0.0;
  lambda[4] = 0.0;
  lambda[5] = 0.0;
  lambda[6] = 0.0;
  lambda[7] = 0.0;
  lambda[8] = 0.0;
}

void MHDSolver::MHDSolver::flux(const double* const Q, double** const F) {
  double V[MHD_VARIABLES];
  Cons2Prim(Q,V);
  double gamma1 = gamma/(gamma-1.0);
  double rho    = V[0];
  double vx     = V[1];
  double vy     = V[2];
  double vz     = V[3];
  double p      = V[4];
  double bx     = V[5];
  double by     = V[6];
  double bz     = V[7];
  double ex     = - (vy*bz - vz*by);
  double ey     = - (vz*bx - vx*bz);
  double ez     = - (vx*by - vy*bx);
  double v2     = vx*vx + vy*vy + vz*vz;
  double b2     = bx*bx + by*by + bz*bz;
  double e2     = ex*ex + ey*ey + ez*ez;
  double lf     = 1.0/std::sqrt(1.0 - v2);
  double w      = rho + gamma1*p;
  double ww     = w*lf*lf;
  double uem    = 0.5*(b2 + e2);
  double wwx    = ww*vx;
  double wwy    = ww*vy;

  //
  F[0][0] = vx*rho*lf;
  F[0][1] = wwx*vx - bx*bx - ex*ex + p + uem;
  F[0][2] = wwx*vy - bx*by - ex*ey;
  F[0][3] = wwx*vz - bx*bz - ex*ez;
  F[0][4] = wwx + (ey*bz - ez*by);
  F[0][5] = V[8];
  F[0][6] = -ez;
  F[0][7] = ey;
  F[0][8] = DivCleaning_a*DivCleaning_a*bx;
  //
  F[1][0] = vy*rho*lf;
  F[1][1] = wwy*vx - by*bx - ey*ex;
  F[1][2] = wwy*vy - by*by - ey*ey + p + uem;
  F[1][3] = wwy*vz - by*bz - ey*ez;
  F[1][4] = wwy + (ez*bx - ex*bz);
  F[1][5] = ez;
  F[1][6] = V[8];
  F[1][7] = -ex;
  F[1][8] = DivCleaning_a*DivCleaning_a*by;
  //
#if DIMENSIONS==3
  double wwz    = ww*vz;

  F[2][0] = vz*rho*lf;
  F[2][1] = wwz*vx - bz*bx - ez*ex;
  F[2][2] = wwz*vy - bz*by - ez*ey;
  F[2][3] = wwz*vz - bz*bz - ez*ez + p + uem;
  F[2][4] = wwz + (ex*by - ey*bx);
  F[2][5] = -ey;
  F[2][6] = ex;
  F[2][7] = V[8];
  F[2][8] = DivCleaning_a*DivCleaning_a*bz;
#endif
}

bool MHDSolver::MHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {
  return tarch::la::equals(t,0.0);

//  // Fixed the following.
//  bool refine;
//  hastoadjustsolution_(&t, &refine);
//  //return _ltob(refine); // something like this is needed
}

void MHDSolver::MHDSolver::AlfvenWave(const double* const x, double* const Q, const double t) {
  // Computes the AlfvenWave at a given time t.
  // Use it ie. with t=0 for initial data
  // Use it for any other time ie. for comparison
  double rho0 = 1.0;
  double p0   = 1.0;
  double eta  = 1.0;
  double B0   = 1.0;
  //;
  double hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0;
  double tempaa = rho0 * hh + B0*B0 * ( 1.0 + eta*eta);
  double tempab = 2.0 * eta * B0*B0 / tempaa;
  double tempac = 0.5 * ( 1.0 + std::sqrt ( 1.0 - tempab*tempab));
  double va2 = B0*B0 / ( tempaa * tempac);
  double vax = std::sqrt ( va2);
  //;
  double BV[3];
  BV[0] = B0;
  BV[1] = eta * B0 * std::cos(2*M_PI*( x[0] - vax*t));
  BV[2] = eta * B0 * std::sin(2*M_PI *( x[0] - vax*t));
  //;
  double VV[3];
  VV[0]   = 0.0;
  VV[1] = - vax * BV[1] / B0;
  VV[2] = - vax * BV[2] / B0;
  //;
  //  Now convert to conservative variables;
  //;
  double V[MHD_VARIABLES];
  V[0] = rho0;
  V[1] = VV[0];
  V[2] = VV[1];
  V[3] = VV[2];
  V[4] = p0;
  V[5] = BV[0];
  V[6] = BV[1];
  V[7] = BV[2];
  V[8] = 0.0;
  Prim2Cons(V,Q);
}

void MHDSolver::MHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* const Q) {
  assertion1(tarch::la::equals(t,0.0),t);

  AlfvenWave(x, Q, t);
}

void MHDSolver::MHDSolver::algebraicSource(const double* const Q, double* const S) {
  // TODO: pass this to Fortran.
  for(int i=0; i < MHD_VARIABLES; i++) {
    S[i] = 0.0;
  }
}



exahype::solvers::Solver::RefinementControl MHDSolver::MHDSolver::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void MHDSolver::MHDSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // TODO: Pass this to Fortran

  // Impose exact boundary conditions
//  AlfvenWave(x, stateOut, t);
//
//  double* F[3];
//  double f[MHD_VARIABLES];
//  double g[MHD_VARIABLES];
//  F[0] = f;
//  F[1] = g;
//
//  #if DIMENSIONS==3
//  double h[MHD_VARIABLES];
//  F[2] = h;
//  #endif
//
//  F[normalNonZero] = fluxOut;
//  flux(stateOut, F);

  // These are the no-boundary conditions:

  for(int i=0; i < MHD_VARIABLES; i++) {
      fluxOut[i] = fluxIn[i];
      stateOut[i] = stateIn[i];
  }

}

