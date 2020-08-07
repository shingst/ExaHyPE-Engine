#include "PDE.h"

 double GAMMA = 1.4;
 double EPSalpha = 1.e-5;
 const int nVar = 9;
 
void PDEPrim2Cons(double* Q,double* V){
   Q[0] = V[0]*V[5];
   Q[1] = V[1]*V[0]*V[5];
   Q[2] = V[2]*V[0]*V[5];
   Q[3] = V[3]*V[0]*V[5];
   Q[4] = V[5]*(V[4]/(GAMMA-1.0)+V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])/2.0);
   Q[5] = V[5];
   Q[6] = V[6];
   Q[7] = V[7];
   Q[8] = V[8];
}

void PDECons2Prim(double* V, const double* const Q){
  const double ialpha = Q[5]/(Q[5]*Q[5] + EPSalpha);
  V[0]   = Q[0]*ialpha;
  const double irho = V[0]/(V[0]*V[0] + EPSalpha);
  V[1] = Q[1]*irho;
  V[2] = Q[2]*irho;
  V[3] = Q[3]*irho;
  V[4] = (GAMMA-1.)*( Q[4]*ialpha - 0.5*V[0]*( V[1]*V[1] + V[2]*V[2] +V[3]*V[3] ) ) ;
  V[5] = Q[5];
  V[6] = Q[6];
  V[7] = Q[7];
}

#include "peano/utils/Loop.h"
void initialdata_(const double* const x,const double t,double* const Q){
    typedef tarch::la::Vector<DIMENSIONS,double> vecNd;
    vecNd xvec(x[0],x[1]);
    vecNd x0(0.75, 1.0/6.0);

      Q[0] = 1.;
      Q[1] = 0.;
      Q[2] = 0.;
      Q[3] = 0.;
      Q[4] = (GAMMA - 1) + exp(-tarch::la::norm2(xvec - x0)/0.25/0.25 ) * 2;
      Q[6] = 0.0;// 3.0*x[1]; //psi
      Q[7] = 0.0;//-3.0*x[0]; //psi
      Q[8] = 0.0; //psi
      double r = x[0] - 1.0/2.0; //signed distance from discontinuity
      auto smoothInterface = [](double r, double dist){
       if(r > dist)
          return 1.0 - 0.01;
       if(r<-dist)
            return 0.0;
       return 0.5*(std::sin(M_PI/2.0/dist*r)+1.0) - 0.01;
      };
      Q[5] = smoothInterface(r,0.01);
      Q[0] = Q[0]*Q[5];
      Q[1] = Q[1]*Q[5];
      Q[4] = Q[4]*Q[5];

}

//standard symmetric 4-digit NACA airfoil
double t = 2*0.594689181;
double symmetric_NACA_airfoil(double x){
    return t/2.0*(0.298222773*std::sqrt(x) - 0.127125232*x - 0.357907906*x*x + 0.291984971*x*x*x - 0.105174606*x*x*x*x); //1.2*(0.2969*std::sqrt(x)-0.126*x-0.3516*x*x+0.2843*x*x*x-0.1036*x*x*x*x);
}

//Use 4612 airfoil values
double m = 0.04; //maximum camber
double p =0.6; //maximum camber location
double c = 1.0;
//introduce cambering
double cambered_NACA_airfoil(double x){
    if(x >=0 && x < p*c)
        return m/p/p*(2*p*(x/c)-(x/c)*(x/c));
    if(x <= c)
        return m/(1.0-p)/(1.0-p)*((1-2*p)+2*p*(x/c)-(x/c)*(x/c));
}

//angle for cambering
double theta_NACA(double x){
    if(x>=0 && x < p*c)
        return 2*m/p/p*(p-(x/c));
    else if(x<= c)
        return 2*m/(1-p)/(1.-p)*(p-x/c);
}

void initialdata_0012(const double* const x,const double t,double* const Q){
    typedef tarch::la::Vector<DIMENSIONS,double> vecNd;
    vecNd xvec(x[0],x[1]);
    vecNd x0(0.75, 1.0/6.0);

    double steepness = 20;
    Q[0] = 0.1   + 0.5 * ( 1.0 + std::tanh(steepness*(x[0]-0.25)) ) * (1.0-0.1);
    Q[4] = 0.125 + 0.5 * ( 1.0 + std::tanh(steepness*(x[0]-0.25)) ) * (1.0-0.125);

    Q[1] = 0.;
    Q[2] = 0.;
    Q[3] = 0.;
    Q[6] = 0.0;
    Q[7] = 0.0;
    Q[8] = 0.0; //psi
    double yu = symmetric_NACA_airfoil(1-x[0]);
    double ru = (x[0]<0 || x[0]>1) ? -1.0 : (yu-std::abs(x[1]));
    auto smoothInterface = [](double r, double dist){
        if(r > dist)
            return 1.0;
        if(r<-dist)
            return 0.0;
        return 0.5*(std::sin(M_PI/2.0/dist*r)+1.0);
    };
    Q[5] =  1.0 - smoothInterface(ru,0.05);

    Q[0] = Q[0]*Q[5];
    Q[1] = Q[1]*Q[5];
    Q[2] = Q[2]*Q[5];
    Q[4] = Q[4]*Q[5];
}


void initialdata(const double* const x,const double t,double* const Q){
    typedef tarch::la::Vector<DIMENSIONS,double> vecNd;
    vecNd xvec(x[0],x[1]);
    vecNd x0(0.75, 1.0/6.0);

    double V[nVar];
    V[0] = 1.4;
    V[1] = 1.0;
    V[2] = 0.;
    V[3] = 0.;
    V[4] = 1.0;
    V[6] = 0.0;
    V[7] = 0.0;
    V[8] = 0.0; //psi
    /* cambered airfoil NACA 2412 */
    double theta = theta_NACA(x[0]);
    double xu = x[0] - symmetric_NACA_airfoil(x[0])*theta/std::sqrt(theta*theta+1);
    double yu = cambered_NACA_airfoil(x[0])+symmetric_NACA_airfoil(x[0])/std::sqrt(theta*theta+1);
    double ru = std::sqrt((xu-x[0])*(xu-x[0])+(yu-x[1])*(yu-x[1]));
    int signu  = -(yu-x[1])/std::abs(yu-x[1]);

    double xl = x[0] + symmetric_NACA_airfoil(x[0])*theta/std::sqrt(theta*theta+1);
    double yl = cambered_NACA_airfoil(x[0])-symmetric_NACA_airfoil(x[0])/(theta*theta+1);
    double rl = std::sqrt((xl-x[0])*(xl-x[0])+(yl-x[1])*(yl-x[1]));
    int signl  = (yl-x[1])/std::abs(yl-x[1]);
    auto smoothInterface = [](double r, double dist){
        if(r > dist)
            return 1.0;
        if(r<-dist)
            return 0.0;
        return 0.5*(std::sin(M_PI/2.0/dist*r)+1.0);
    };

    /*std::cout << "yt at 1 " << symmetric_NACA_airfoil(1) << std::endl;
    std::cout << "yc at 1 " << cambered_NACA_airfoil(1) << std::endl;
    std::cout << "yt at 0 " << symmetric_NACA_airfoil(0) << std::endl;
    std::cout << "yc at 0 " << cambered_NACA_airfoil(0) << std::endl;*/

    V[5] = 1.0;
    if(x[0] >=0.0 && x[0] <= 1.0)
        V[5] = std::min(V[5], ((x[1]>=0) ? smoothInterface( signu*ru, 0.005) : smoothInterface(signl*rl, 0.005) ));
    PDEPrim2Cons(Q,V);
}

void PDEflux(const double* const Q,double** const F){
    double V[nVar];
    PDECons2Prim(V,Q);

    double p = V[4];
    F[0][0] = V[5]*V[0]*V[1];
    F[0][1] = V[5]*( V[0]*V[1]*V[1] + p );
    F[0][2] = V[5]*V[0]*V[1]*V[2];
    F[0][3] = V[5]*V[0]*V[1]*V[3];
    F[0][4] = V[1]*(Q[4] + V[5]*p);
    F[0][5] = 0.; 
    F[0][6] = 0.; 
    F[0][7] = 0.; 
    F[0][8] = 0.; 

    F[1][0] = V[5]*V[0]*V[2]; 
    F[1][1] = V[5]*V[0]*V[2]*V[1];
    F[1][2] = V[5]*( V[0]*V[2]*V[2] + p ); 
    F[1][3] = V[5]*V[0]*V[2]*V[3];
    F[1][4] = V[2]*(Q[4] + V[5]*p);
    F[1][5] = 0.; 
    F[1][6] = 0.; 
    F[1][7] = 0.; 
    F[1][8] = 0.; 
}


void PDEncp(const double* const Q,const double* const gradQ,double* const BgradQ){
    double V[nVar];
    PDECons2Prim(V,Q);

    BgradQ[0] = 0.0;
    BgradQ[1] = -V[4] * gradQ[0*nVar+5];
    BgradQ[2] = -V[4] * gradQ[1*nVar+5];
    BgradQ[3] = 0.0;
    BgradQ[4] = -V[4] * (V[6]*gradQ[0*nVar+5]+V[7]*gradQ[1*nVar+5]);
    BgradQ[5] = V[6]*gradQ[0*nVar+5]+V[7]*gradQ[1*nVar+5] ;
    BgradQ[6] = 0.0;
    BgradQ[7] = 0.0;
    BgradQ[8] = 0.0;
}

void PDEEigenvalues(const double* const Q,const int direction,double* const L){
    if(Q[5]<1.e-4){
        L[0] = 1.0;
        return;
    }
    double V[nVar];
    PDECons2Prim(V,Q);

    double u  =  V[1+direction];
    double vs =  V[6+direction];
    double c = std::sqrt( GAMMA*V[4]/V[0] );
    L[0]   = u-c;
    L[1]   = u;
    L[2]   = u;
    L[3]   = u;
    L[4]   = u+c;
    L[5]   = vs;
    L[6] = 0.;
    L[7] = 0.;
    L[8] = 0.;
    double sum = 0; // continue in spite of issues
    for(int i=0; i < nVar; i++) sum += std::abs(L[i]);
    if(sum < 1e-10) L[0] = 1;
}
