#ifndef EXAHYPE_GRMHD_PIZZA_TOV_ID
#define EXAHYPE_GRMHD_PIZZA_TOV_ID

#include "InitialData/InitialData.h"

// Forward declaration for the actual PizzaTOV implementation.
namespace Pizza {
namespace TOV {
	class pointwise_tov;
}// ns Pizza
}// ns TOV

/**
 * The PizzaTOV initial data are only accessible if support is compiled into the
 * ExaHyPE binary.
 * 
 * This class uses the pimpl mechanism and has a stub in the PizzaTOV.cpp implementation
 * in order to allow a seamless compilation in any way.
 **/
struct pizzatov : public idobj {
	Pizza::TOV::pointwise_tov* tov;
	
	pizzatov();
	void Interpolate(const double* x, double t, double* Q);
};

#endif /* EXAHYPE_GRMHD_PIZZA_TOV_ID */
