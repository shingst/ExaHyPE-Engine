#ifndef __INITIAL_DATA_ADAPTER_GRMHD_cpp__
#define __INITIAL_DATA_ADAPTER_GRMHD_cpp__

#include "PDE/PDE-GRMHD-ExaHyPE.h"
#include "Parameters/mexa.h"
#include <cmath> 
#include <string>

/// An Interface to initial data creation.
struct InitialDataCode {
	/// ID codes can use this to read and validate parameters
	virtual void readParameters(const mexa::mexafile& parameters) {}
	/// ID codes can use this to run or prepare ID. This is called once.
	virtual void prepare() {}
	/// Give the initial data at a spacetime point. This is called massively parallel
	/// and the only function which must be implemented by this interface.
	virtual void Interpolate(const double* x, double t, double* Q) = 0;
	
	/// This is the actual initial data registration
	static InitialDataCode* getInstanceByName(std::string name);
};

// a shorthand, also for debugging
void InitialData(const double* x, double t, double* Q);

// global initial data
class GlobalInitialData {
	tarch::logging::Log _log;
public:
	/// The pointer is by default null, meaning no ID can be used. The setIdByName method
	/// sets this id pointer, finally.
	InitialDataCode* id;
	std::string name; ///< set by setIdByName
	bool alreadyPrepared; ///< set by prepare()
	
	GlobalInitialData(InitialDataCode* _id, std::string _name) : _log("GlobalInitialData"), id(_id), name(_name), alreadyPrepared(false) {}
	
	/// Returns true in case of success
	bool setIdByName(const std::string& name);
	
	/// Tries to read ID from parameters. In case of problems, raises exceptions.
	/// Chainable.
	GlobalInitialData& setByParameters(const mexa::mexafile& parameters);
	
	/// Run ID codes. Chainable.
	/// This tracks whether the preparation already took place, so you can
	/// call it multiple times.
	GlobalInitialData& prepare();
	
	/// Gives a singleton instance
	static GlobalInitialData& getInstance();
	
	// as a shorthand
	static InitialDataCode& getInitialDataCode() { return *(getInstance().id); }
	
	/// Ensure that globalID are available. Returns if available or not.
	/// aborts the code if doCrash is set.
	static bool ensureGlobalIDareSet(bool doCrash=true);
};

#include "PizzaTOV.h"
#include "RNSID.h"
#include "AnalyticalID.h"



#endif /* __INITIAL_DATA_ADAPTER_GRMHD_cpp__ */
