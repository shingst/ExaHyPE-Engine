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

/* BEGIN Fortran definitions */

/// An enum-like structure for named magic numbers in Fortran.
/// The storages have values there.
const extern struct F_tInitialDataEnum {
	int CCZ4GaugeWave;
	int CCZ4Puncture;
	int CCZ4TwoPunctures;
	int CCZ4GRMHDAccretion;
	int CCZ4Kerr2D;
} typesDef_ICType;

/// The actual storage to tell Fortran what ID to run.
extern int typesDef_ic;

/* END Fortran definitions */

InitialDataCode* InitialDataCode::getInstanceByName(std::string idname) {
	toLower(idname);
	if(idname == "gaugewave") 	return new FortranPointwiseID(typesDef_ICType.CCZ4GaugeWave);
	if(idname == "puncture")	return new FortranPointwiseID(typesDef_ICType.CCZ4Puncture);
	if(idname == "twopunctures") 	return new FortranPointwiseID(typesDef_ICType.CCZ4TwoPunctures);
	if(idname == "grmhdaccretion") 	return new FortranPointwiseID(typesDef_ICType.CCZ4GRMHDAccretion);
	if(idname == "kerr2d") 		return new FortranPointwiseID(typesDef_ICType.CCZ4Kerr2D);
	return nullptr;
}

FortranPointwiseID::FortranPointwiseID(int ic_type) {
	typesDef_ic = ic_type;
}

GlobalInitialData& GlobalInitialData::getInstance() {
	static GlobalInitialData* me = nullptr;
	if(!me) me = new GlobalInitialData(nullptr, "null");
	return *me;
}

bool GlobalInitialData::setIdByName(const std::string& _name) {
	static tarch::logging::Log _log("GlobalInitialData");
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

void GlobalInitialData::setByParameters(const mexa::mexafile& parameters) {
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
}

/// a shorthand
void InitialData(const double* x, double t, double* Q) {
	if(GlobalInitialData::getInstance().id == nullptr) {
		static tarch::logging::Log _log("InitialData");
		logError("InitialData()", "Cannot access InitialData because no initial Data has been defined yet.");
		std::abort();
	}
	
	GlobalInitialData::getInitialDataCode().adjustPointSolution(x,t,Q);
}

