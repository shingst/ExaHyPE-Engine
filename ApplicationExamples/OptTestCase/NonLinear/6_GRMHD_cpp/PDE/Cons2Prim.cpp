#include "PDE.h"
#include <cmath> // max, sqrt
#include <algorithm> // max
#include <utility> // swap

#include <stdio.h> // printf

using namespace SVEC;
using namespace std; // sqrt, max

// The C2P-Routines ported from Fortran for the SRMHD by Dominic, fully correct also for GRMHD.
double RTSAFE_C2P_RMHD1(const double X1, const double X2, const double XACC, const double gam, const double d, const double e, const double s2, const double b2, const double sb2, double& w,bool& failed);
void FUNC_C2P_RMHD1(const double x,double& f,double& df,const double gam,const double d,const double e,const double s2,const double b2, const double sb2,double& w);

// At this stage, we can prepare conservative variables but do not yet have access to primitive ones.
void GRMHD::Cons2Prim::prepare() {
	// B_i: Needed for Sij and preparation:
	// v^i: is needed for the PDE system (Sij) and is either computed in the followup() or by the Cons2Prim operation.
	// S^i: is needed in both flux and ncp
	Bmag.lo=0; DFOR(i) CONTRACT(j) Bmag.lo(i) += gam.lo(i,j) * Bmag.up(j);
	Si.up  =0; DFOR(i) CONTRACT(j) Si.up(i)   += gam.up(i,j) * Si.lo(j);
	BmagBmag = 0; CONTRACT(k) BmagBmag += Bmag.lo(k)*Bmag.up(k); // B^j * B_j // needed for ptot
	SconScon = 0; CONTRACT(k) SconScon +=   Si.lo(k)*Si.up(k);   // S^j * S_j // needed for c2p
	BmagScon = 0; CONTRACT(k) BmagScon += Bmag.lo(k)*Si.up(k);   // B^j * S_j // needed for c2p
}

/**
 * Debug the C2P computation.
 **/
constexpr bool debug_c2p = false;

/**
 * This is a *major* switch whether to quietly ignore C2P failures (as needed for
 * a limiting scheme where failures shall be caught) or fail verbosely. We could
 * do anything beyond (like: Passing errors etc).
 *
 * -> locally, it can be chosen on the constructor.
 **/
constexpr bool global_force_crash_on_failure = false;

#define S_STR_HELPER(x) #x
#define S_STR(x) S_STR_HELPER(x)
#define SL(x,l) printf(l " = %e\n", x);
#define S(x) SL(x, #x)
// https://stackoverflow.com/a/5459929 for expanding before insertion
#define SI(x) { DFOR(i) SL(x(i), #x "(" S_STR(i) ")"); }
#define SIsym(x) { SYMFOR(i,j) SL(x(i,j), #x "(" S_STR(i) "," S_STR(j) ")"); }

// At this stage, we have all primitives recovered and can postcompute some quantities. This is as a
// service or can be used for 
void GRMHD::Cons2Prim::followup() {
	BmagVel = 0;  CONTRACT(j) BmagVel  += Bmag.up(j)*vel.lo(j);  // B^j * v_j // needed for ptot
	WLorentz = Dens/rho;
	WW = SQ(WLorentz); // W^2: Lorentz factor squared
	ptot = press + 0.5*(BmagBmag/WW + SQ(BmagVel)); // total pressure incl magn. field, needed in 3-energy-mom-tensor


	double VelVel2=0; CONTRACT(j) VelVel2  += vel.up(j)*vel.lo(j);
	double WLorentz2 = 1.0 / std::sqrt(1.0 - VelVel2); // another way of computing WLorentz

	if(debug_c2p) { S(WLorentz); S(WW); S(Dens/rho); S(WLorentz2); S(VelVel2); }
}

void GRMHD::Cons2Prim::copyFullStateVector() {
	// 1) Assume that hydro variables have been set
	// 2) Copy only:
	copy_magneto(V);
	copy_adm(V);
}

// small helper
bool isEnvSet(const char* str) { return std::getenv(str); }

void GRMHD::Cons2Prim::report_c2p(const char* const context) {
	printf("%s: C2P Failure at x=[%f,%f,%f]\n", context, Q[19],Q[20],Q[21]);
}

void GRMHD::Cons2Prim::perform() {
	// Experiment: Checking whether atmo values are suitable for flooring.
	constexpr double p_floor   = atmo_press; // 1.0e-5;
	constexpr double rho_floor = atmo_rho; // 1.0e-4;
	//constexpr double gamma = gamma;

	// RTSAFE gives us x = v^2, y = rho * h * Gamma^2.
	failed = rtsafe(VelVel, RhoEnthWW);
	if(debug_c2p) { printf("C2P Outcome: "); S(VelVel); S(RhoEnthWW);}
	if (failed) {
		if(crash_on_failure || global_force_crash_on_failure) {
			// We should raise an error instead, the c2p failed.
			printf("C++ C2P FAILED.\n");
			printf("Input conserved State Vector Q: \n");
			// makes no sense with densitied state
			//TDO(d,MHD::size) printf("Q[%d] = %e\n", d, Q[d]);
			S(Dens); S(tau); SI(Si.lo); SI(Si.up);
			SI(Bmag.up); S(phi);
			S(BmagBmag); S(SconScon);
			
			// metric stuff which should not matter anyway because it is
			// constant. Anyway let us print the coordinates where this happens:
			S(gam.det); S(gam.sqdet); S(alpha);
			SI(beta.up); SIsym(gam.up);
			// debugging: Coordinates
			SL(Q[19], "Coord:x");
			SL(Q[20], "Coord:y");
			SL(Q[21], "Coord:z");
			SL(Q[22], "proof");
			
			std::abort();
		} else {
			// Silently do nothing.
			// use report() to print something.
			
			// if(!isEnvSet("EXAHYPE_SILENT_C2P"))
			//	printf("C2P Failure at x=[%f,%f,%f]\n", Q[19],Q[20],Q[21]);
		}
		
		// this is (probably!) the atmosphere treatment:
		rho = rho_floor;
		press = p_floor;
		DFOR(i) vel.up(i) = 0;
	} else {
		rho  = Dens*std::sqrt(1.-VelVel);
		DFOR(i) vel.up(i) = (Si.up(i) + BmagScon/RhoEnthWW * Bmag.up(i)) / (RhoEnthWW+BmagBmag);
		//DFOR(i) vel.lo(i) = (Si.lo(i) + BmagScon/RhoEnthWW * Bmag.lo(i)) / (RhoEnthWW+BmagBmag);
		//vel.up=0; DFOR(i) CONTRACT(j)  vel.up(i) +=  vel.lo(j) * gam.up(i,j);
		vel.lo=0; DFOR(i) CONTRACT(j)  vel.lo(i) +=  vel.up(j) * gam.lo(i,j);
		press     = (gamma-1)/gamma*((1.-VelVel)*RhoEnthWW - rho); // Ideal EOS
		press     = max(p_floor, press); // bracketing // TODO: was atmo_press before.
		if(debug_c2p) { S(rho); SI(vel.up); SI(vel.lo); S(press); }
	}
}

bool GRMHD::Cons2Prim::rtsafe(double& x, double& y) {
	constexpr double tol    = 1e-8;
	constexpr double eps    = 1.e-10;
	constexpr double x1     = 0.;
	constexpr double x2     = 1.-eps;
	
	constexpr double gamma = GRMHD::Parameters::gamma;
	//double e      = tau;// + Dens; // sic! My naming is probably sick.
	double e      = tau + Dens;

	// This is ECHO paper, first option
	// [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)]
	// x = v^2, y = rho * h * Gamma^2.
	bool failed=false;
	y = 0;
	x = RTSAFE_C2P_RMHD1(x1,x2,tol,(gamma-1.0)/gamma,Dens,e,SconScon,BmagBmag,BmagScon*BmagScon,y,failed);
	return failed;
}

double RTSAFE_C2P_RMHD1(const double X1, const double X2, const double XACC, const double gam, const double d, const double e, const double s2, const double b2, const double sb2, double& w,bool& failed) {
  constexpr int MAXIT=200;
  double FL, DF;

  failed = false;
  FUNC_C2P_RMHD1(X1,FL,DF,gam,d,e,s2,b2,sb2,w);

  double v2 = NAN;
  if(FL==0) {
    v2=X1;
    return v2;
  }

  double FH;
  FUNC_C2P_RMHD1(X2,FH,DF,gam,d,e,s2,b2,sb2,w);

  if(FH==0) {
     v2=X2;
     return v2;
  }
  if(FL*FH>0) {
     failed = true;
     return v2;
  }

  double XL,XH;
  if(FL<0) {
    XL=X1;
    XH=X2;
  } else {
    XH=X1;
    XL=X2;
    std::swap(FL,FH);
  }
  v2=.5*(X1+X2);
  double DXOLD=std::abs(X2-X1);
  double DX=DXOLD;

  double F;
  FUNC_C2P_RMHD1(v2,F,DF,gam,d,e,s2,b2,sb2,w); // Dominic made an error here: &DF,&DF
  for (int J=1; J<MAXIT; ++J) {
     if(((v2-XH)*DF-F)*((v2-XL)*DF-F)>=0
          || std::abs(2.*F)>std::abs(DXOLD*DF) ) {
        DXOLD=DX;
        DX=0.5*(XH-XL);
        v2=XL+DX;
        if(XL==v2) {
          return v2;
        }
     } else {
        DXOLD=DX;
        DX=F/DF;
        double TEMP=v2;
        v2=v2-DX;
        if (TEMP==v2) {
          return v2;
        }
     }
     if (std::abs(DX)<XACC) {
       return v2;
     }
     FUNC_C2P_RMHD1(v2,F,DF,gam,d,e,s2,b2,sb2,w);
     if(F<0) {
        XL=v2;
        FL=F;
     } else {
        XH=v2;
        FH=F;
     }
  }
  failed = true;
  return v2;
}


void FUNC_C2P_RMHD1(const double x,double& f,double& df,const double gam,const double d,const double e,const double s2,const double b2, const double sb2,double& wo) {
  //
  // This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  // and it corresponds to their choice 3 in Section 3.2
  //
  constexpr int maxiter = 200; // in F it is 100

  const double v2=x;
  const double rho=d*std::sqrt(1.-v2);

  const double c3=1.-gam*(1.-v2);
  const double c2=gam*rho+.5*b2*(1.+v2)-e;
  const double c0=-0.5*sb2;

  //  For every x=v2, we solve a cubic in W of the form:
  //  c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  //  W=y of the paper. If sb=0 ( meaning c0 = 0),
  //  w = -c2/c3 > 0 and dw = 0 in the do loop below.
  //  If -c2/c3 < 0 when sb=0, which makes w=0,
  //  this is a signature that something was wrong before.
  double w;
  if (std::abs(c0) < 1.0e-20) {
    w = -c2 / c3;
  } else {
    w = std::max(-c2/c3, std::pow(-c0/c3, 1./3.));
    for(int iter = 0; iter < maxiter; iter++) {
      double dw = -((c3*w + c2)*w*w + c0)/((3*c3*w + 2*c2)*w);
      if (std::abs(dw/w)<1.e-10) {
        break; // happy breakdown
      }
      w = w + dw;
    }
  }

  const double dc3   = gam;
  const double dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2));
  const double dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2);
  const double wb    = w + b2;
  const double vb2 = sb2 / (w*w);
  
  // set output:
  f   = wb*wb * v2 - ( 2.0 * w + b2) * vb2 - s2;
  df  = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2));
  wo  = w;

//  std::cout << "wo="<<*wo;
//  std::cout << ",f="<<*f;
//  std::cout << ",df="<<*df << std::endl;
}
