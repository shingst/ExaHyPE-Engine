/**
 * This file is a C++ drop-in replacement for C2PRoutines.f90.
 * You can switch between these two files by putting to one of the files a suffix ".disabled",
 * wipe the *.mk lists and recompile.
 * 
 * This code was written by Dominic in Oct2016 for the SRMHD_cpp application and adopted for the GRMHD
 * by Sven at Jun2017.
 **/

#include <cmath> // abs
#include <algorithm> // max

#if 0

// For the actual Fortran interface, look at the end of this file.

namespace C2PinC {

// This is the helper function for RTSAFE_C2P_RMHD1, we never expect it to be
// called from Fortran so we stick to an ordinary C signature.
void FUNC_C2P_RMHD1(const double x,double* f,double* df,const double gam,const double d,const double e,const double s2,const double b2, const double sb2,double* w_out) {
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
  if ( std::abs ( c0) < 1.0e-20) {
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

// This is the C implemenation of RTSAFE_C2P_RMHD1 which is not directly callable by Fortran
// due to the return value and the parameter handling. Instead, use the RTSAFE_C2P_RMHD1_
// function given below
double RTSAFE_C2P_RMHD1(double X1,double X2,double XACC,double gam,double d,
    double e,double s2,double b2,double sb2,double* w,bool* failed) {
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
  FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1_result,&F,&DF,gam,d,e,s2,b2,sb2,w);
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

} // ns C2PinC


typedef const double* const in;

// Use this to switch between Fortran and C++ version.
#if 1

extern "C" {
// Dropin replacement for
// RECURSIVE FUNCTION RTSAFE_C2P_RMHD1(X1,X2,XACC,gam,d,e,s2,b2,sb2,w,FAILED)
// Thanks to the conversion between double* <-> double, we don't have to change aynthing in Fortran
// when calling this function
double rtsafe_c2p_rmhd1_(in X1,in X2,in XACC,in gam,in d,in e,in s2,in b2,in sb2,double* w,bool* failed) {
	return C2PinC::RTSAFE_C2P_RMHD1(*X1, *X2, *XACC, *gam, *d, *e, *s2, *b2, *sb2, w, failed);
}

} /* Extern C */

#endif

#endif
