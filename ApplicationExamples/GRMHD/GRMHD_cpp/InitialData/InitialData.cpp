
#include "InitialData/InitialData.h"
#include "tarch/logging/Log.h"

#include <algorithm> // transform
#include <cctype> // tolower

/// String to lowercase, inplace.
namespace IDtools {
	
void toLower(std::string& data) {
	std::transform(data.begin(), data.end(), data.begin(), ::tolower);
}

} // IDtools

using namespace IDtools;

InitialDataCode* InitialDataCode::getInstanceByName(std::string idname) {
	toLower(idname);
	if(idname == "alfenwave") 	return new AnalyticalID<AlfenWaveCons>;
	if(idname == "pizzatov")	return new pizzatov();
	if(idname == "rnsid") 		return new rnsid();
	if(idname == "shocktube") 	return new GRMHD_Shocktube();
	if(idname == "tovsolver")	return new TovSolverAdapter();
	if(idname == "toriid")		return new ToriIDAdapter();
	return nullptr;
}

GlobalInitialData& GlobalInitialData::getInstance() {
	static GlobalInitialData* me = nullptr;
	if(!me) me = new GlobalInitialData(nullptr, "null");
	return *me;
}

bool GlobalInitialData::setIdByName(const std::string& _name) {
	//static tarch::logging::Log _log("GlobalInitialData");
	name = _name;
	if(id) logWarning("setIdByName()", "Have already set global initial data, now overwriting with "<< name << ".");
	
	id = InitialDataCode::getInstanceByName(name);
	if(id) {
		logInfo("setIdByName()", "Successfully loaded "<< name << " Initial Data.");
		return true;
	} else {
		logError("setIdByName()", "Requested Initial Data '"<< name << "' not known.");
		return false;
	}
}

GlobalInitialData& GlobalInitialData::setByParameters(const mexa::mexafile& parameters) {
	if(alreadySetParameters) {
		// happens anways only in the LimitingContext. Can be removed once we are sure everything runs correctly.
		logInfo("setByParameters(mexafile&)", "Global unique initial data have already read the parameters. In order to avoid inconsistencies, we don't allow this two times.");
		return *this;
	} else {
		alreadySetParameters = true;
	}
	
	std::string idseckey = "initialdata";
	std::string idnamekey = "name";
	mexa::mexafile idparam = parameters(idseckey);
	if(!idparam.contains(idnamekey)) {
		logError("setByParameters()", "For setting up the initila data, I need a section " << idseckey  << " in the parameters, as well as the key " << idseckey << "/" << idnamekey << " to be set to a valid initial data function. Instead, I got these parameters: " << parameters.toString());
		std::abort();
	}
	bool couldSetId = setIdByName(idparam(idnamekey).get_string());
	if(!couldSetId) {
		logError("setByParameters()", "Could not create Initial Data. Cannot solve an initial value problem without initial data.");
		std::abort();
	} else {
		id->readParameters(idparam);
	}
	return *this;
}

GlobalInitialData& GlobalInitialData::prepare() {
	if(!alreadyPrepared) {
		id->prepare();
	} else {
		logInfo("prepare()", "Global unique inital data have already been prepared (for instance by the FV or DG counterpart in the LimitingADERDGSolver)");
	}
	alreadyPrepared = true;
	return *this;
}

bool GlobalInitialData::ensureGlobalIDareSet(bool doCrash) {
	bool isSet = (getInstance().id != nullptr);
	if(doCrash && !isSet) {
		static tarch::logging::Log _log("InitialData");
		logError("InitialData()", "Cannot access InitialData because no initial Data has been defined yet.");
		std::abort();
	}
	return isSet;
}

#include "PDE/PDE.h"
#include "GRMHDSolver_ADERDG_Variables.h"
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;

/// a shorthand
void InitialData(const double* x, double t, double* Q) {
	for(int i=0;i<nVar;i++) Q[i] = 0.0; // Zeroize
	
        GlobalInitialData::ensureGlobalIDareSet(/*doCrash=*/true);	
	GlobalInitialData::getInitialDataCode().Interpolate(x,t,Q);
	
	  // also store the positions for debugging
	GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts Qi;
	constexpr int posI = 19; // position 19
	for(int i=0; i<DIMENSIONS; i++) Q[posI+i] = x[i];
	if(DIMENSIONS<3) Q[posI+2] = -1;
	var.check() = 47110815;
	//zero2Din3D(Q);

	for(int i=0;i<nVar;i++) { if(!std::isfinite(Q[i])) { printf("Qid[%d] = %e\n", i, Q[i]); std::abort(); } }
}
