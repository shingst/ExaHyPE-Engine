/**
 * This file is a pure bugfix to circumvent C to FORTRAN parameter passing (especially
 * string passing). As I have not found a better way to push the outcome of
 * constants->getValueAsString() to FORTRAN, the FORTRAN initial data function calls
 * this C function which again in turn calls the FORTRAN functions. This is an ugly
 * workaround for a better future solution.
 **/

#include "InitialDataAdapter.h"
#include "MHDSolver.h"

// storage for the constants pointer
exahype::parser::ParserView* constants;

static const char* specfile_initialdata_key = "initialdata";

typedef void (*idfunc)(double* x, double* Q);

static tarch::logging::Log _log("MHDSolver");

// a workaround way to access parameters via environment variables while
// the exahype parameter stuff is not working.
static const std::string envPrefix = "exahype_";
//char* mygetenv( const char* env_var ) { printf("Trying to access key '%s'\n", env_var); return std::getenv(env_var); }
bool envKeyExists(const std::string& key) { return std::getenv(key.c_str()) != nullptr; }
std::string getEnvValue(const std::string& key) { assert(envKeyExists(key)); return std::getenv(key.c_str()); }
bool shallUseEnvWorkaround() { return envKeyExists(envPrefix+"parameter_workaround"); }
bool wrap_isValueValidString(const std::string& key) {
	if(shallUseEnvWorkaround()) return envKeyExists(envPrefix+key);
	assert(constants != nullptr);
	return constants->isValueValidString(key);
}
	
std::string wrap_getValueAsString(const std::string& key) {
	if(shallUseEnvWorkaround()) return getEnvValue(envPrefix+key);
	assert(constants != nullptr);
	return constants->getValueAsString(key);
}
// end of workaround.

// log some stuff only once. Another hacky thing
int atLeastOnceWarnedAboutEnv = 0;
int atLeastOnceWarnedAboutFallback = 0;
int atLeastOnceInformAboutSuccess = 0;
#define HASNOT(counter) (!counter++)

extern "C" {
	
// C function called by FORTRAN
void initialdatabyexahypespecfile(double* x, double* Q) {
	if(shallUseEnvWorkaround() && HASNOT(atLeastOnceWarnedAboutEnv))
		logWarning("InitialDatabyExahyPESpecFile()", "Using ENV workaround to determine parameters");
	if(!shallUseEnvWorkaround() && !constants && HASNOT(atLeastOnceWarnedAboutEnv)) {
		logError("InitialDatabyExahyPESpecFile()", "Parser instance is Null!");
	} else if(!wrap_isValueValidString(specfile_initialdata_key)) {
		logError("InitialDatabyExahyPESpecFile()", "Initial Data Value is not valid");
	} else {
		std::string id_name = wrap_getValueAsString(specfile_initialdata_key);
		// c++ string to lower
		std::transform(id_name.begin(), id_name.end(), id_name.begin(), ::tolower);
		
		     if(id_name == "alfenwave") initialalfenwave_(x,Q);
		else if(id_name == "blast")     initialblast_(x,Q);
		else if(id_name == "orsagtang") initialorsagtang_(x,Q);
		else if(id_name == "rotor")     initialrotor_(x,Q);
		else {
			logError("InitialDatabyExahyPESpecFile()", "Unknown Initial Data value '"<< id_name <<"'");
		}
		
		if(HASNOT(atLeastOnceInformAboutSuccess)) {
			logWarning("InitialDatabyExahyPESpecFile()", "Successfully loaded Initial Data '"<< id_name << "'");
		}
		
		return; // success. or so.
	}
	// failure: Do fallback initial data
	if(HASNOT(atLeastOnceWarnedAboutFallback))
		logWarning( "InitialDatabyExahyPESpecFile()", "Falling back to AlfenWave");
	initialalfenwave_(x,Q);
}

}/* extern "C" */
