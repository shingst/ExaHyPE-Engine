#ifndef __INITIAL_DATA_ADAPTER_GRMHD__
#define __INITIAL_DATA_ADAPTER_GRMHD__

#include <string>

// Interface to initial data
struct idobj {
	virtual void Interpolate(const double* x, double t, double* Q) = 0;
};

/* note this is never deleted, so it is a memory leak, but it's only 
created once, so that's okay right? */
extern idobj *id; // a global accessible pointer


/**
 *a function in InitialData.cpp which prepares the ID, both accessible
 * from a pure ADERDG, pure FV or limiting application
 * 
 * @returns True when the preperation had success/was not neccessary, false
 *          in case of failure.
 **/
bool prepare_id(std::string idname);

#endif /* __INITIAL_DATA_ADAPTER_GRMHD__ */
