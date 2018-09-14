#include "MyNavierStokesSolver_FV.h"
#include "MyNavierStokesSolver_FV_Variables.h"

const double GAMMA = 1.4;
const double Pr = 0.75; //Prandtl
const double cv = 1.0;
const double c0 = 10.0; //adiabatic sound speed, choose high for incompressible limit
const double kappa = 0; //for now
const double mu = 1.0/1600;

tarch::logging::Log NavierStokesADERDG::MyNavierStokesSolver_FV::_log( "NavierStokesFV::MyNavierStokesSolver" );

void NavierStokesADERDG::MyNavierStokesSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void NavierStokesADERDG::MyNavierStokesSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
    if ( tarch::la::equals( t,0.0 ) ) {
        Variables vars(Q);
        vars.rho() = 1.0;
        if(std::abs(x[1]-1.0)<1e-6){
            vars.j(0, 1.0 ,0);
        }
        double p = c0*c0/GAMMA;
        vars.E() = p/(GAMMA-1) + 0.5 * (vars.j()*vars.j())/vars.rho();
    }
}

void NavierStokesADERDG::MyNavierStokesSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 + #parameters
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double irho = vars.rho()/(vars.rho()*vars.rho() + 1e-14);
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() ); //pressure

  double u_n = Q[dIndex + 1] * irho;
  double c  = std::sqrt(GAMMA * p * irho);

  eigs.rho() = u_n - c;
  eigs.j(u_n,u_n,u_n);
  eigs.E()  = u_n + c;
}

void NavierStokesADERDG::MyNavierStokesSolver_FV::viscousEigenvalues(const double* const Q, const int dIndex, double* lambda) {
    // Dimensions             = 2
    // Number of variables    = 5 + #parameters
    ReadOnlyVariables vars(Q);
    Variables eigs(lambda);
    const double irho = vars.rho()/(vars.rho()*vars.rho() + 1e-14);
    // @todo Please implement/augment if required
    lambda[0] = 0.0;
    lambda[1] = 0.0;
    lambda[2] = (4.0/3.0)*mu*irho;
    lambda[3] = mu*irho;
    lambda[4] = mu*GAMMA/Pr*irho;
}

void NavierStokesADERDG::MyNavierStokesSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside){
    // Dimensions             = 2
    // Number of variables    = 5 + #parameters
    ReadOnlyVariables varsInside(stateInside);
    Variables         varsOutside(stateOutside);

    varsOutside = varsInside;
    if(faceIndex == 3){
        stateOutside[1] = 0.0 - stateInside[1];
        stateOutside[2] = 2.0 - stateInside[2];
        stateOutside[3] = 0.0 - stateInside[3];
    }
    else{
        stateOutside[1] = 0.0 - stateInside[1];
        stateOutside[2] = 0.0 - stateInside[2];
        stateOutside[3] = 0.0 - stateInside[3];
    }
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

void NavierStokesADERDG::MyNavierStokesSolver_FV::viscousFlux(const double* const Q, const double* const gradQ, double** F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 5 + 0
  ReadOnlyVariables vars(Q);

  double irho = Q[0]/(Q[0]*Q[0] + 1e-14);
  const double p = (GAMMA-1) * (vars.E() - 0.5 * irho * vars.j()*vars.j() );

  //Non-viscous part
  F[0][0] = Q[1];
  F[0][1] = irho*Q[1]*Q[1] + p;
  F[0][2] = irho*Q[1]*Q[2];
  F[0][3] = irho*Q[1]*Q[3];
  F[0][4] = irho*Q[1]*(Q[4] + p);

  F[1][0] = Q[2];
  F[1][1] = irho*Q[2]*Q[1];
  F[1][2] = irho*Q[2]*Q[2] + p;
  F[1][3] = irho*Q[2]*Q[3];
  F[1][4] = irho*Q[2]*(Q[4] + p);

  irho  = 1./Q[0];
  double uu    = Q[1]*irho;
  double vv    = Q[2]*irho;
  double ww    = Q[3]*irho;

  double uux  = irho*(gradQ[1] - uu*gradQ[0]);
  double vvx  = irho*(gradQ[2] - vv*gradQ[0]);
  double wwx  = irho*(gradQ[3] - ww*gradQ[0]);

  double uuy  = irho*(gradQ[1+5] - uu*gradQ[0+5]);
  double vvy  = irho*(gradQ[2+5] - vv*gradQ[0+5]);
  double wwy  = irho*(gradQ[3+5] - ww*gradQ[0+5]);

  double iRho2 = irho*irho;
  double iRho3 = iRho2*irho;

  double dTdW1 = - Q[4]*iRho2 + iRho3*( Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] );
  double dTdW2 = - Q[1]*iRho2;
  double dTdW3 = - Q[2]*iRho2;
  double dTdW4 = - Q[3]*iRho2;

  const double icv = 1.0/cv;
  double Tx = icv*( dTdW1*gradQ[0]   + dTdW2*gradQ[1]   + dTdW3*gradQ[2]   + dTdW4*gradQ[3]   + irho*gradQ[4] );
  double Ty = icv*( dTdW1*gradQ[0+5] + dTdW2*gradQ[1+5] + dTdW3*gradQ[2+5] + dTdW4*gradQ[3+5] + irho*gradQ[4+5] );

  double divV23  = 2./3.*(uux + vvy);

  //viscous part
  //F[0][0] += 0.;
  F[0][1] += mu*( 2*uux - divV23 );
  F[0][2] += mu*(   uuy + vvx    );
  F[0][3] += mu*(  0 + wwx    ); //uuz <-0
  F[0][4] += F[0][1]*uu + F[0][2]*vv + F[0][3]*ww + kappa*Tx;

  //F[1][0] += 0.
  F[1][1] += F[0][2];
  F[1][2] += mu*( 2*vvy - divV23 );
  F[1][3] += mu*(   0 + wwy    ); //vvz<-0
  F[1][4] += F[1][1]*uu + F[1][2]*vv + F[1][3]*ww + kappa*Ty;
}



