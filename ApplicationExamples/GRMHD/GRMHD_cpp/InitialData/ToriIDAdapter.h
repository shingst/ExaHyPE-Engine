#ifndef __EXAHYPE_GRMHD_CPP_TORIID_ADAPTER__
#define __EXAHYPE_GRMHD_CPP_TORIID_ADAPTER__

#include "InitialData.h"
#include "tarch/logging/Log.h"

class ToriIDAdapter : public InitialDataCode {
	tarch::logging::Log _log;

public:
	bool hasBeenPrepared; ///< a guard to ensure the preparation took place
	
	ToriIDAdapter() : _log("TorIIDAdapter"), hasBeenPrepared(false) {}
	
	void prepare() override;
	void Interpolate(const double* x, double t, double* Q) override;
};


#endif /* __EXAHYPE_GRMHD_CPP_TORIID_ADAPTER__ */
