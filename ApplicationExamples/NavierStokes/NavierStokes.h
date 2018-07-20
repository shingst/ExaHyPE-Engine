#include <cmath>
#include "tarch/la/Vector.h"

namespace NavierStokes {
  class NavierStokes;
}

class NavierStokes::NavierStokes {
public:
  NavierStokes(double referenceT, double referenceViscosity, double sutherlandC) :
    sutherlandC(sutherlandC),
    sutherlandLambda(
		     (referenceViscosity * (referenceT + sutherlandC))/
		     std::pow(referenceT, 3./2)){

    
  }

  double evaluateTemperature(double rho, double pressure) const {
    return pressure/(gasConstant + rho);
  }
  
  double evaluateHeatConductionCoeff(double viscosity) const {
    return 1./Pr * viscosity * GAMMA * (1.0/(GAMMA - 1)) * gasConstant;
  }
  
  double evaluatePressure(double E, double rho, tarch::la::Vector<3,double> j) const {
    return (GAMMA-1) * (rho * E - 0.5 * rho * j * j);
  }

  double evaluateViscosity(double T) const {
    // Sutherland's law
    return sutherlandLambda * (std::pow(T, 3./2) / (T + sutherlandC));
  }

private:
  const double GAMMA = 1.4;
  const double gasConstant = 8.3144598; // TODO
  const double sutherlandC;
  const double sutherlandLambda;
  const double Pr = 0.71; // Prandtl number, for air
  
};
