/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
#include "exahype/parser/Parser.h"

#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>
#include <string>
#include <regex>
#include <cstdlib> // getenv, exit
#include <sstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <exception>
#include <cstdarg>

#include "tarch/la/ScalarOperations.h"

#include "exahype/parser/ParserView.h"


#include "json.hpp" // this is an ~1MB header file which is included here exclusively
using json = nlohmann::json;

tarch::logging::Log exahype::parser::Parser::_log("exahype::parser::Parser");

namespace exahype {
namespace parser {
namespace tools {

// Using a local namespace to avoid conflicts with these common function names.

// a buffer-overflow-safe version of sprintf
// source: http://stackoverflow.com/a/69911
std::string vformat (const char *fmt, va_list ap) {
  using namespace std;
  // Allocate a buffer on the stack that's big enough for us almost
  // all the time.  Be prepared to allocate dynamically if it doesn't fit.
  size_t size = 1024;
  char stackbuf[1024];
  std::vector<char> dynamicbuf;
  char *buf = &stackbuf[0];
  va_list ap_copy;

  while (1) {
      // Try to vsnprintf into our buffer.
      va_copy(ap_copy, ap);
      int needed = vsnprintf (buf, size, fmt, ap);
      va_end(ap_copy);

      // NB. C99 (which modern Linux and OS X follow) says vsnprintf
      // failure returns the length it would have needed.  But older
      // glibc and current Windows return -1 for failure, i.e., not
      // telling us how much was needed.

      if (needed <= (int)size && needed >= 0) {
          // It fit fine so we're done.
          return std::string (buf, (size_t) needed);
      }

      // vsnprintf reported that it wanted to write more characters
      // than we allotted.  So try again using a dynamic buffer.  This
      // doesn't happen very often if we chose our initial size well.
      size = (needed > 0) ? (needed+1) : (size*2);
      dynamicbuf.resize (size);
      buf = &dynamicbuf[0];
  }
}

std::string sformat(const char *fmt, ...) {
  using namespace std;
  va_list ap;
  va_start (ap, fmt);
  std::string buf = vformat (fmt, ap);
  va_end (ap);
  return buf;
}

} // ns tools
} // ns parser
} // ns exahype

using namespace exahype::parser::tools;

struct exahype::parser::ParserImpl {
  json data;

  void read(std::istream& is) {
    data = json::parse(is);
  }

  void print() {
      std::cout << std::setw(4) << data << std::endl;
  }
  
  bool isValid() {
    return true; // TODO: Wether json holds something or not
  }

  /**
   * Generic method to extract various types (like int,double,bool,string,vector),
   * throws std::runtime_error in case of error.
   **/
  template<typename T>
  T getFromPath(std::string path, T defaultValue, bool isOptional=false) const {
    std::stringstream ss;
    try {
      return data.at(json::json_pointer(path));
    } catch (json::type_error& e) {
      ss << path << " is not a " << typeid(T).name()   << " (" << e.what() << ")";
    } catch(json::out_of_range& e) {
      if(isOptional)
        return defaultValue;
      else
        ss << "Missing entry " << path << " (" << e.what() << ")";
    }
    throw std::runtime_error(ss.str());
  }

};

bool exahype::parser::Parser::_interpretationErrorOccured(false);

const std::string exahype::parser::Parser::_noTokenFound("notoken");

double exahype::parser::Parser::getValueFromPropertyString(
    const std::string& parameterString, const std::string& key) {
  std::size_t startIndex = parameterString.find(key);
  startIndex = parameterString.find(":", startIndex);
  std::size_t endIndex = parameterString.find_first_of("}, \n\r", startIndex + 1);

  std::string substring =
      parameterString.substr(startIndex + 1, endIndex - startIndex - 1);

  double result;
  std::istringstream ss(substring);
  ss >> result;

  if (ss) {
    return result;
  } else {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

exahype::parser::Parser::Parser() {
  _impl = new exahype::parser::ParserImpl();

  _identifier2Type.insert(
      std::pair<std::string, exahype::solvers::Solver::Type>(
          "ADER-DG", exahype::solvers::Solver::Type::ADERDG));
  _identifier2Type.insert(
      std::pair<std::string, exahype::solvers::Solver::Type>(
          "Finite-Volumes", exahype::solvers::Solver::Type::FiniteVolumes));
  _identifier2Type.insert(
        std::pair<std::string, exahype::solvers::Solver::Type>(
            "Limiting-ADER-DG", exahype::solvers::Solver::Type::LimitingADERDG));

  _identifier2TimeStepping.insert(
      std::pair<std::string, exahype::solvers::Solver::TimeStepping>(
          "global", exahype::solvers::Solver::TimeStepping::Global));
  _identifier2TimeStepping.insert(
      std::pair<std::string, exahype::solvers::Solver::TimeStepping>(
          "globalfixed", exahype::solvers::Solver::TimeStepping::GlobalFixed));
}



void exahype::parser::Parser::readFile(const std::string& filename) {
    std::ifstream inputFile;
    inputFile.open(filename.c_str());
    if (!inputFile.good()) {
      logError("readFile(String)", "cannot open file " << filename);
      invalidate();
      return;
    }
    readFile(inputFile, filename);
}

void exahype::parser::Parser::readFile(std::istream& inputFile, std::string filename) {
  _filename = filename;
  _impl->read(inputFile);

  //  For debugging purposes
  if(std::getenv("EXAHYPE_VERBOSE_PARSER")) { // runtime debugging
      _impl->print();
  }

  checkValidity();
  //  std::string configuration = getMPIConfiguration();
  //  int ranksPerNode = static_cast<int>(exahype::parser::Parser::getValueFromPropertyString(configuration,"ranks_per_node"));
  //  std::cout << "ranks_per_node="<<ranksPerNode << std::endl;
}

void exahype::parser::Parser::checkValidity() {
  // functions have side-effects: might set _interpretationErrorOccured
  getDomainSize();
  getOffset();
  if (foundSimulationEndTime()) {
    getSimulationEndTime();
  } else {
    getSimulationTimeSteps();
  }
}

bool exahype::parser::Parser::isValid() const {
  return !_impl->isValid() && !_interpretationErrorOccured;
}

void exahype::parser::Parser::invalidate() const {
  invalidate();
}

bool exahype::parser::Parser::hasPath(std::string path) const {
  // Q: I'm not sure whethe the iterator end() is also true for nested structures (path pointer)
  bool found = _impl->data.find(json::json_pointer(path)) != _impl->data.end();
  return found;
}

std::string exahype::parser::Parser::dumpPath(std::string path) const {
  try {
    return _impl->data.find(json::json_pointer(path))->dump();
  } catch (json::type_error& e) {
    logError("dumpPath()", "Path " << path << " contains a non-UTF-8-encodable string (" << e.what() << ")");
  } catch(json::out_of_range& e) {
    logError("dumpPath()", "Path " << path << " not found (" << e.what() << ")");
  }
  invalidate();
  return "{'dumpPath':'error'}";
}

bool exahype::parser::Parser::isValueValidBool(const std::string& path) const {
  return hasPath(path) && _impl->data.find(json::json_pointer(path))->is_boolean();
}

bool exahype::parser::Parser::isValueValidInt(const std::string& path) const {
  return hasPath(path) && _impl->data.find(json::json_pointer(path))->is_number_integer();
}

bool exahype::parser::Parser::isValueValidDouble(const std::string& path) const {
  return hasPath(path) && _impl->data.find(json::json_pointer(path))->is_number();
}

bool exahype::parser::Parser::isValueValidString(const std::string& path) const {
  return hasPath(path) && _impl->data.find(json::json_pointer(path))->is_string();
}

std::string exahype::parser::Parser::getStringFromPath(std::string path, std::string defaultValue, bool isOptional) const {
  assertion(isValid());
  try {
    return _impl->getFromPath(path, defaultValue, isOptional);
  } catch(std::runtime_error& e) {
    logError("getStringFromPath()", e.what());
    invalidate();
    return defaultValue; /* I don't like returning something here */
  }
  /*
  assertion(isValid());
  try {
    return _impl->data.at(json::json_pointer(path));
  } catch (json::type_error& e) {
    logError("getStringFromPath()", path << " is not a string (" << e.what() << ")");
  } catch(json::out_of_range& e) {
    logError("getStringFromPath()", "Missing entry " << path << " (" << e.what() << ")");
  }
  invalidate();
  return "";
  */
}

int exahype::parser::Parser::getIntFromPath(std::string path, int defaultValue, bool isOptional) const {
  assertion(isValid());
  try {
    return _impl->getFromPath(path, defaultValue, isOptional);
  } catch(std::runtime_error& e) {
    logError("getIntFromPath()", e.what());
    invalidate();
    return defaultValue; /* I don't like returning something here */
  }
}

double exahype::parser::Parser::getDoubleFromPath(std::string path, double defaultValue, bool isOptional) const {
  assertion(isValid());
  try {
    return _impl->getFromPath(path, defaultValue, isOptional);
  } catch(std::runtime_error& e) {
    logError("getIntFromPath()", e.what());
    invalidate();
    return defaultValue; /* I don't like returning something here */
  }
}

bool exahype::parser::Parser::getBoolFromPath(std::string path, bool defaultValue, bool isOptional) const {
  assertion(isValid());
  try {
    return _impl->getFromPath(path, defaultValue, isOptional);
  } catch(std::runtime_error& e) {
    logError("getBoolFromPath()", e.what());
    invalidate();
    return defaultValue; /* I don't like returning something here */
  }
}

tarch::la::Vector<DIMENSIONS,double> exahype::parser::Parser::getDimVectorFromPath(std::string path) const {
  assertion(isValid());
  tarch::la::Vector<DIMENSIONS,double> result;
  // Here we do not rely on _impl->getFromPath because we don't want to
  // copy from a STL vector. Somehow.
  try {
    json::json_pointer p(path);
    result(0) = _impl->data.at(p).at(0);
    if(_impl->data.at(p).size() != DIMENSIONS) {
      logError("getDimVectorFromPath()", path << " holds a vector of size " << _impl->data.at(p).size() << ", however we have " << DIMENSIONS << " spatial dimensions");
    }
    result(1) = _impl->data.at(p).at(1);
    if(DIMENSIONS == 3)
      result(2) = _impl->data.at(p).at(2);
    return result;
  } catch(json::type_error& e) {
    logError("getDimVectorFromPath()", path << " holds not a double-vector of size " << DIMENSIONS << " (" << e.what() << ")");
  } catch(json::out_of_range& e) {
    logError("getDimVectorFromPath()", "Missing entry " << path << " (" << e.what() << ")");
  }

  invalidate();
  return result;
}

std::vector<int> exahype::parser::Parser::getIntVectorFromPath(std::string path) const {
  assertion(isValid());
  std::vector<int> empty;
  try {
    return _impl->getFromPath(path, empty, isMandatory);
  } catch(std::runtime_error& e) {
    logError("getIntVectorFromPath()", e.what());
    invalidate();
    return empty; /* I don't like returning something here */
  }
}

bool exahype::parser::Parser::flagListContains(std::string path, std::string keyword) const {
  assertion(isValid());
  try {
    auto p = json::json_pointer(path);
    if(_impl->data.count(p)) {
      auto j = _impl->data.at(p);
      if(!j.is_array()) {
        logError("flagListContains()", "Expected " << path << " to hold an array, but this is not the case.");
        invalidate();
        return false;
      }
      for (json::iterator it = j.begin(); it != j.end(); ++it) {
        if(! it->is_string() ) {
          logError("flagListContains()", path << " holds a non-string element in its array-content.");
          invalidate();
          return false;
        }
        if( it->get<std::string>().compare(keyword) == 0) return true;
      }
      return false; // not found
    } else {
      logDebug("flagListContains()", "Flag array " << path << " is not existing, defaulting to empty.");
    }
  } catch(json::type_error& e) {
    logError("flagListContains()", path << " holds weird data (" << e.what() << ")");
  } catch(json::out_of_range& e) {
    logError("flagListContains()", "Missing something below or at  " << path << " (" << e.what() << ")");
  }
  
  invalidate();
  return false; // Should not return data
}

int exahype::parser::Parser::getNumberOfThreads() const {
  return getIntFromPath("/shared_memory/cores");
}

tarch::la::Vector<DIMENSIONS, double> exahype::parser::Parser::getDomainSize() const {
  assertion(isValid());
  tarch::la::Vector<DIMENSIONS, double> result;
  
  int dim = getIntFromPath("/computational_domain/dimension");
  
  if(dim != DIMENSIONS) {
    logError("getDomainSize()",
              "dimension: value "<< dim << " in specification file" <<
              " does not match -DDim"<<DIMENSIONS<<" switch in Makefile. Rerun toolkit!");
    invalidate();
    return result;
  }
  
  result = getDimVectorFromPath("/computational_domain/width");
  logDebug("getDomainSize()", "found size " << result);
  return result;
}

tarch::la::Vector<DIMENSIONS, double> exahype::parser::Parser::getOffset() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;
  result = getDimVectorFromPath("/computational_domain/offset");
  logDebug("getOffset()", "found offset " << result);
  return result;
}

std::string exahype::parser::Parser::getMulticorePropertiesFile() const {
  std::string result = getStringFromPath("/shared_memory/properties_file");
  logDebug("getMulticorePropertiesFile()", "found " << result);
  return result;
}

exahype::parser::Parser::MPILoadBalancingType exahype::parser::Parser::getMPILoadBalancingType() const {
  std::string token = getStringFromPath("/distributed_memory/identifier");
  exahype::parser::Parser::MPILoadBalancingType result = MPILoadBalancingType::Static;
  if (token.compare("static_load_balancing") == 0) {
    result = MPILoadBalancingType::Static;
  } else {
    logError("getMPILoadBalancingType()",
             "Invalid distributed memory identifier " << token);
    invalidate();
  }
  return result;
}


bool exahype::parser::Parser::MPIConfigurationContains(std::string flag) const {
  return flagListContains("/distributed_memory/configure", flag);
}

bool exahype::parser::Parser::SharedMemoryConfigurationContains(std::string flag) const {
  return flagListContains("/shared_memory/configure", flag);
}


double exahype::parser::Parser::getNodePoolAnsweringTimeout() const {
  const std::string path = "/distributed_memory/max_node_pool_answering_time";
  const double defaultResult = 1e-2;
  double result = getDoubleFromPath(path, defaultResult, true);
  if(tarch::la::equals(defaultResult, result)) {
    logWarning( "getNodePoolAnsweringTimeout()", path << " not specified for MPI configuration so use default timeout of " << defaultResult );
  }
  return result;
}


int exahype::parser::Parser::getMPIBufferSize() const {
  int result = getIntFromPath("distributed_memory/buffer_size");

  // Apparently, in former days an invalid value just yielded in a non-fatal error.
  // all non-castable ints resulted in negative numbers.

  if(result <= 0) {
    logError("getMPIBufferSize()", "Invalid MPI buffer size " << result);
    result = 64;
    invalidate();
  }

  return result;
}

int exahype::parser::Parser::getMPITimeOut() const {
  double result = getIntFromPath("distributed_memory/timeout");

  // Apparently, in former days an invalid value just yielded in a non-fatal error.
  // all non-castable doubles resulted in negative numbers.

  if (result <= 0) {
    logError("getMPIBufferSize()", "Invalid MPI timeout value " << result);
    result = 0;
    invalidate();
  }
  
  return result;
}

exahype::parser::Parser::MulticoreOracleType exahype::parser::Parser::getMulticoreOracleType() const {
  std::string token = getStringFromPath("shared-memory/identifier");
  exahype::parser::Parser::MulticoreOracleType result = MulticoreOracleType::Dummy;
  if (token.compare("dummy") == 0) {
    result = MulticoreOracleType::Dummy;
  } else if (token.compare("autotuning") == 0) {
    result = MulticoreOracleType::AutotuningWithRestartAndLearning;
  } else if (token.compare("autotuning_without_learning") == 0) {
    result = MulticoreOracleType::AutotuningWithoutLearning;
  } else if (token.compare("autotuning_without_restart") == 0) {
    result = MulticoreOracleType::AutotuningWithLearningButWithoutRestart;
  } else if (token.compare("sampling") == 0) {
    result = MulticoreOracleType::GrainSizeSampling;
  } else {
    logError("getMulticoreOracleType()", "Invalid shared memory identifier "
                                             << token);
    result = MulticoreOracleType::Dummy;
    invalidate();
  }
  return result;
}

double exahype::parser::Parser::getSimulationEndTime() const {
  double result = getDoubleFromPath("computational_domain/end_time");
  logDebug("getSimulationEndTime()", "found end time " << result);
  
  // Apparently, in former days an invalid value just yielded in a non-fatal error.
  // all non-castable doubles resulted in negative numbers.

  if (result <= 0) {
    logError("getSimulationEndTime()",
             "Invalid simulation end-time: " << result);
    result = 1.0;
    invalidate();
  }
  return result;
}

bool exahype::parser::Parser::foundSimulationEndTime() const {
  const double not_there = -43;
  double result = getDoubleFromPath("computational_domain/end_time", not_there, isOptional);
  bool found = !tarch::la::equals(result, not_there);
  return found;
}

int exahype::parser::Parser::getSimulationTimeSteps() const {
  int result = getIntFromPath("computational_domain/time_steps");
  logDebug("getSimulationEndTime()", "found result " << result);
  
  // Apparently, in former days an invalid value just yielded in a non-fatal error.
  // all non-castable ints resulted in negative numbers.

  if (result < 0) {
    logError("getSimulationEndTime()",
             "Invalid simulation timestep: " << result);
    invalidate();
  }
  return result;
}

bool exahype::parser::Parser::getFuseAlgorithmicSteps() const {
  const bool default_value = false;
  bool result = getBoolFromPath("optimisation/fuse_algorithmic_steps", default_value, isOptional);
  return result;
}



double exahype::parser::Parser::getFuseAlgorithmicStepsFactor() const {
  const double default_value = 0.0;
  double result = getDoubleFromPath("optimisation/fuse_algorithmic_steps_factor", default_value, isOptional);
  logDebug("getFuseAlgorithmicStepsFactor()", "found fuse-algorithmic-steps-factor " << result);
  if(result < 0.0 || result > 1.0) {
    logError("getFuseAlgorithmicStepsFactor()",
              "'fuse-algorithmic-steps-factor': Value must be greater than zero "
              "and smaller than one: "
                  << result);
    result = 0.0;
    invalidate();
  }
  return result;
}

bool exahype::parser::Parser::getSpawnPredictionAsBackgroundThread() const {
  return getBoolFromPath("/optimisation/spawn_predictor_as_background_thread", false, isOptional);
}

bool exahype::parser::Parser::getSpawnAMRBackgroundThreads() const {
  return getBoolFromPath("/optimisation/spawn_amr_background_threads", false, isOptional);
}

bool exahype::parser::Parser::getDisableMetadataExchangeInBatchedTimeSteps() const {
  return getBoolFromPath("/optimisation/disable_metadata_exchange_in_batched_time_steps", false, isOptional);
}

bool exahype::parser::Parser::getDisablePeanoNeighbourExchangeInTimeSteps() const {
  return getBoolFromPath("/optimisation/disable_vertex_exchange_in_time_steps", false, isOptional);
}

double exahype::parser::Parser::getTimestepBatchFactor() const {
  double result = getDoubleFromPath("/optimisation/time_step_batch_factor", 0.0, isOptional);
  logDebug("getFuseAlgorithmicStepsFactor()", "found time-step-batch-factor " << result);

  if (result < 0.0 || result > 1.0) {
    logError("getFuseAlgorithmicStepsFactor()",
              "'time-step-batch-factor': Value is required in global-optimisation "
              "section and must be greater than zero and smaller than one: "
                  << result);
    result = 0.0;
    invalidate();
  }
  return result;
}


bool exahype::parser::Parser::hasOptimisationSegment() const {
  return hasPath("/optimisation");
}


bool exahype::parser::Parser::getSkipReductionInBatchedTimeSteps() const {
  if (hasOptimisationSegment()) {
    return tarch::la::greater(getTimestepBatchFactor(),0.0);
  }
  else return false;
}


double exahype::parser::Parser::getDoubleCompressionFactor() const {
  double result = getDoubleFromPath("optimisation/double_compression", 0.0, isOptional);

  if (result < 0.0) {
    logError("getDoubleCompressionFactor()",
            "'double-compression': Value is required in global-optimisation "
            "section and must be greater than or equal to zero: " << result);
    result = 0.0;
    invalidate();
  }

  return result;
}


bool exahype::parser::Parser::getSpawnDoubleCompressionAsBackgroundTask() const {
  return getBoolFromPath("optimisation/spawn_double_compression_as_background_thread");
}


exahype::solvers::Solver::Type exahype::parser::Parser::getType(
    int solverNumber) const {
  exahype::solvers::Solver::Type result =
      exahype::solvers::Solver::Type::ADERDG;
  std::string token = getStringFromPath(sformat("solver/%d/type", solverNumber+1));
  if (_identifier2Type.find(token) != _identifier2Type.end()) {
    result = _identifier2Type.at(token);
    logDebug("getType()", "found type " << exahype::solvers::Solver::toString(result));
  } else {
    logError(
        "getType()",
        "'" << getIdentifier(solverNumber) << "': 'type': Value '" << token
            << "' is invalid. See the ExaHyPE documentation for valid values.");
    invalidate();
  }
  return result;
}

std::string exahype::parser::Parser::getIdentifier(int solverNumber) const {
  std::string token = getStringFromPath(sformat("solver/%d/name", solverNumber+1));
  logDebug("getIdentifier()", "found identifier " << token);
  return token;
}

int exahype::parser::Parser::getVariables(int solverNumber) const {
  std::string path = sformat("solver/%d/variables", solverNumber+1);
  auto j = _impl->data.at(path);
  if(j.is_primitive()) {
    // variables=N
    return getIntFromPath(path);
  } else if(j.is_array()) {
    // variables=["foo","bar","baz"], meaning variables=3
    return j.size();
  } else {
    // count multiplicities
    int result = 0;
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
      // the lazy way, constructing the global path again. Side note: We require multiplicity, we don't support an implicit 1.
      result += getIntFromPath(sformat("solver/%d/variables/%d/multiplicity", solverNumber+1, it - j.begin()));
    }
    if(result == 0) {
      logError("getVariables()",
               "'" << getIdentifier(solverNumber)
               << "': 'variables': Value must be greater than zero.");
          invalidate();
    }
    return result;
  }
}

int exahype::parser::Parser::getParameters(int solverNumber) const {
  // Copied the getVariables code here.
  
  std::string path = sformat("solver/%d/parameters", solverNumber+1);
  auto j = _impl->data.at(path);
  if(j.is_primitive()) {
    // variables=N
    return getIntFromPath(path);
  } else if(j.is_array()) {
    // variables=["foo","bar","baz"], meaning variables=3
    return j.size();
  } else {
    // count multiplicities
    int result = 0;
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
      // the lazy way, constructing the global path again. Side note: We require multiplicity, we don't support an implicit 1.
      result += getIntFromPath(sformat("solver/%d/parameters/%d/multiplicity", solverNumber+1, it - j.begin()));
    }
    if(result == 0) {
      logError("getParameters()",
               "'" << getIdentifier(solverNumber)
               << "': 'parameters': Value must be greater than zero.");
          invalidate();
    }
    return result;
  }
}

int exahype::parser::Parser::getOrder(int solverNumber) const {
  return getIntFromPath(sformat("solver/%d/order", solverNumber+1));
}


double exahype::parser::Parser::getMaximumMeshSize(int solverNumber) const {
  double result = getDoubleFromPath(sformat("solver/%d/maximum_mesh_size", solverNumber+1));

  if (tarch::la::smallerEquals(result, 0.0)) {
    logError("getMaximumMeshSize(int)",
             "'" << getIdentifier(solverNumber)
                 << "': 'maximum-mesh-size': Value must be greater than zero.");
    invalidate();
  }

  logDebug("getMaximumMeshSize()", "found maximum mesh size " << result);
  return result;
}

int exahype::parser::Parser::getMaximumMeshDepth(int solverNumber) const {
  int result = getIntFromPath(sformat("solver/%d/maximum_mesh_depth", solverNumber+1));

  if (tarch::la::smaller(result, 0)) {
    logError("getMaximumMeshDepth(int)",
             "'" << getIdentifier(solverNumber)
                 << "': 'maximum-mesh-depth': Value must be greater than or equal to zero.");
    invalidate();
  }

  logDebug("getMaximumMeshDepth()", "found maximum mesh size " << result);
  return result;
}

int exahype::parser::Parser::getHaloCells(int solverNumber) const {
  int result = getIntFromPath(sformat("solver/%d/halo_cells", solverNumber+1));

  if (tarch::la::smaller(result, 0)) {
    logError("getHaloCells(int)",
             "'" << getIdentifier(solverNumber)
                 << "': 'halo-cells': Value must be greater than or equal to zero.");
    invalidate();
  }

  logDebug("getHaloCells()", "found halo-cells " << result);
  return result;
}

int exahype::parser::Parser::getRegularisedFineGridLevels(int solverNumber) const {
  int result = getIntFromPath(sformat("solver/%d/regularised_fine_grid_levels", solverNumber+1));
  
  if (tarch::la::smaller(result, 0)) {
    logError("getRegularisedFineGridLevels(int)",
             "'" << getIdentifier(solverNumber)
                 << "': 'regularised-fine-grid-levels': Value must be greater than or equal to zero.");
    invalidate();
  }

  logDebug("getRegularisedFineGridLevels()", "found regularised-fine-grid-levels " << result);
  return result;
}

exahype::solvers::Solver::TimeStepping exahype::parser::Parser::getTimeStepping(
    int solverNumber) const {
  exahype::solvers::Solver::TimeStepping result;
  const std::string default_value = "global";
  std::string token = getStringFromPath(sformat("solver/%d/time_stepping", solverNumber+1), default_value, isOptional);
  if (_identifier2TimeStepping.find(token) != _identifier2TimeStepping.end()) {
    result = _identifier2TimeStepping.at(token);
    // logDebug("getTimeStepping()", "found TimeStepping " << result);
    logDebug("getTimeStepping()", "found TimeStepping "<< token);
    return result;
  } else {
    logError(
        "getTimeStepping()",
        "'" << getIdentifier(solverNumber) << "': 'time-stepping': Value '"
            << token
            << "' is invalid. See the ExaHyPE documentation for valid values.");
    invalidate();
  }
  return exahype::solvers::Solver::TimeStepping::Global; /* keep in sync with default_value above */
}

double exahype::parser::Parser::getDMPRelaxationParameter(int solverNumber) const {
  double result = getDoubleFromPath(sformat("solver/%d/dmp_relaxation_parameter", solverNumber+1));

  if (result < 0) {
    logError("getDMPRelaxationParameter()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-relaxation-parameter': Value must not be negative.");
    invalidate();
  }

  logDebug("getParameters()", "found dmp-relaxation-parameter " << result);
  return result;
}

double exahype::parser::Parser::getDMPDifferenceScaling(int solverNumber) const {
  double result = getDoubleFromPath(sformat("solver/%d/dmp_difference_scaling", solverNumber+1));

  if (result < 0) {
    logError("getDMPDifferenceScaling()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-difference-scaling': Value must not be negative.");
    invalidate();
  }

  logDebug("getDMPDifferenceScaling()", "found dmp-difference-scaling " << result);
  return result;
}

int exahype::parser::Parser::getDMPObservables(int solverNumber) const {
  int result = getIntFromPath(sformat("solver/%d/dmp_observables", solverNumber+1));

  if (result < 0) {
    logError("getDMPObservables()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-observables': Value must not be negative.");
    invalidate();
  }

  logDebug("getDMPObservables()", "found dmp-observables " << result);
  return result;
}

int exahype::parser::Parser::getStepsTillCured(int solverNumber) const {
  const int default_value = 0;
  int result = getIntFromPath(sformat("solver/%d/steps_till_cured", solverNumber+1), default_value, isOptional);

  if (result < 0) {
    logError("getStepsTillCured()",
              "'" << getIdentifier(solverNumber)
              << "': 'steps-till-cured': Value must be integral and not negative.");
    invalidate();
  }

  logDebug("getStepsTillCured()", "found steps-till-cured " << result);
  return result;
}

int exahype::parser::Parser::getLimiterHelperLayers(int solverNumber) const {
  return getIntFromPath(sformat("solver/%d/help_layers", solverNumber+1), 1, isOptional);
}

std::string exahype::parser::Parser::getIdentifierForPlotter(int solverNumber,
                                                     int plotterNumber) const {
  return getStringFromPath(sformat("solver/%d/plotters/%d/name", solverNumber+1, plotterNumber+1));
}

std::string exahype::parser::Parser::getNameForPlotter(int solverNumber,
                                               int plotterNumber) const {
  return getStringFromPath(sformat("solver/%d/plotters/%d/type", solverNumber+1, plotterNumber+1));
}

int exahype::parser::Parser::getUnknownsForPlotter(int solverNumber,
                                           int plotterNumber) const {
  return getIntFromPath(sformat("solver/%d/plotters/%d/variables", solverNumber+1, plotterNumber+1));
}

double exahype::parser::Parser::getFirstSnapshotTimeForPlotter(
    int solverNumber, int plotterNumber) const {
  return getDoubleFromPath(sformat("solver/%d/plotters/%d/time", solverNumber+1, plotterNumber+1));
}

double exahype::parser::Parser::getRepeatTimeForPlotter(int solverNumber,
                                                int plotterNumber) const {
  return getDoubleFromPath(sformat("solver/%d/plotters/%d/repeat", solverNumber+1, plotterNumber+1));
}

std::string exahype::parser::Parser::getFilenameForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  return getStringFromPath(sformat("solver/%d/plotters/%d/output", solverNumber+1, plotterNumber+1));
}

exahype::parser::ParserView exahype::parser::Parser::getSelectorForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  return exahype::parser::ParserView(this,sformat("/solver/%d/plotters/%d/select", solverNumber+1, plotterNumber+1));
}

std::string exahype::parser::Parser::getLogFileName() const {
  return getStringFromPath("/paths/log_file", "", isOptional);
}

std::string exahype::parser::Parser::getProfilerIdentifier() const {
  return getStringFromPath("/profiling/profiler", "NoOpProfiler", isOptional);
}

std::string exahype::parser::Parser::getMetricsIdentifierList() const {
  return getStringFromPath("/profiling/metrics", "{}", isOptional);
}

std::string exahype::parser::Parser::getProfilingOutputFilename() const {
  return getStringFromPath("/profiling/profiling_output", "", isOptional);
}

void exahype::parser::Parser::logSolverDetails(int solverNumber) const {
  logInfo("logSolverDetails()",
          "Solver " << exahype::solvers::Solver::toString(getType(solverNumber)) << " "
                    << getIdentifier(solverNumber) << ":");
  logInfo("logSolverDetails()", "variables:\t\t" << getVariables(solverNumber));
  logInfo("logSolverDetails()", "parameters:\t" << getParameters(solverNumber));
  logInfo("logSolverDetails()", "order:\t\t" << getOrder(solverNumber));
  logInfo("logSolverDetails()", "maximum-mesh-size:\t"
                                    << getMaximumMeshSize(solverNumber));
  logInfo("logSolverDetails()",
          "time-stepping:\t" << exahype::solvers::Solver::toString(getTimeStepping(solverNumber)));
}

void exahype::parser::Parser::checkSolverConsistency(int solverNumber) const {
  assertion1(solverNumber <
                 static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
             solverNumber);
  exahype::solvers::Solver* solver =
      exahype::solvers::RegisteredSolvers[solverNumber];

  bool recompile = false;
  bool runToolkitAgain = false;
  if (solver->getType() != getType(solverNumber)) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Solver type in specification file "
		 << "('" << exahype::solvers::Solver::toString(getType(solverNumber)) << "') "
                 << "differs from solver type used in implementation "
		 << "('" << exahype::solvers::Solver::toString(solver->getType()) << "').");
    recompile = true;
    invalidate();
  }

  if (solver->getIdentifier().compare(getIdentifier(solverNumber))) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Identifier in specification file "
                 << "('" << getIdentifier(solverNumber)
                 << "') differs from identifier used in implementation ('"
                 << solver->getIdentifier() << "').");
    recompile = true;
    invalidate();
  }

  if (solver->getNumberOfVariables() != getVariables(solverNumber)) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Value for 'variables' in specification file"
                 << "('" << getVariables(solverNumber)
                 << "') differs from number of variables used in "
                    "implementation file ('"
                 << solver->getNumberOfVariables() << "').");
    recompile = true;
    invalidate();
  }

  if (solver->getNumberOfParameters() != getParameters(solverNumber)) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Value for field 'parameters' in specification file"
                 << "('" << getParameters(solverNumber)
                 << "') differs from  number of parameters used in "
                    "implementation file ('"
                 << solver->getNumberOfParameters() << "').");
    recompile = true;
    invalidate();
  }

  if (solver->getType() == exahype::solvers::Solver::Type::ADERDG &&
      solver->getNodesPerCoordinateAxis() != getOrder(solverNumber) + 1) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Value for field 'order' in specification file "
                 << "('" << getOrder(solverNumber)
                 << "') differs from value used in implementation file ('"
                 << solver->getNodesPerCoordinateAxis() - 1 << "'). ");
    runToolkitAgain = true;
    invalidate();
  }

  // @todo We should add checks for FV as well
  
  // (Sven:) somehow for me the following lines are never printed. I don't
  // know why.

  if (runToolkitAgain) {
    logError("checkSolverConsistency",
             "Please (1) run the Toolkit again for " << getSpecfileName() << ", and (2) recompile!");
    invalidate();
  }

  if (recompile) {
    logError(
        "checkSolverConsistency",
        "Please (1) adjust the specification file (" << getSpecfileName() <<  ") or the file '"
            << solver->getIdentifier()
            << ".cpp' (and similar) accordingly, and (2) recompile!");
    invalidate();
  }
}

std::string exahype::parser::Parser::getSpecfileName() const {
  return _filename;
}

int exahype::parser::Parser::getRanksPerNode() {
  return getIntFromPath("/distributed_memory/ranks_per_node");
}


int exahype::parser::Parser::getNumberOfBackgroundTasks() {
  int result = getIntFromPath("/shared_memory/background_tasks");

  if (result<-2) {
    logWarning("getNumberOfBackgroundTasks()", "invalid number of background tasks (background-tasks field in configuration) " <<
      "set or no number at all. Use default (1). See BackgroundTasks.h for documentation.");
    result = 1;
  }
  return result;
}


bool exahype::parser::Parser::useManualPinning() {
  return flagListContains("/shared_memory/configure", "manual_pinning");
}

exahype::parser::ParserView exahype::parser::Parser::createParserView(const int solverNumberInSpecificationFile) {
  return exahype::parser::ParserView(this,sformat("/solver/%d/constants", solverNumberInSpecificationFile));
}


exahype::parser::Parser::TBBInvadeStrategy exahype::parser::Parser::getTBBInvadeStrategy() const {
  if (flagListContains("/shared_memory/configure", "no_invade")) return TBBInvadeStrategy::NoInvade;
  if (flagListContains("/shared_memory/configure", "analyse_optimal_static_distribution_but_do_not_invade")) return TBBInvadeStrategy::NoInvadeButAnalyseDistribution;
  if (flagListContains("/shared_memory/configure", "occupy_all_cores")) return TBBInvadeStrategy::OccupyAllCores;
  if (flagListContains("/shared_memory/configure", "invade_between_time_steps")) return TBBInvadeStrategy::InvadeBetweenTimeSteps;
  if (flagListContains("/shared_memory/configure", "invade_throughout_computation")) return TBBInvadeStrategy::InvadeThroughoutComputation;
  if (flagListContains("/shared_memory/configure", "invade_at_time_step_startup_plus_throughout_computation")) return TBBInvadeStrategy::InvadeAtTimeStepStartupPlusThroughoutComputation;

  return TBBInvadeStrategy::Undef;
}
