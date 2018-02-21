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

#include "tarch/la/ScalarOperations.h"

#include "exahype/parser/ParserView.h"

tarch::logging::Log exahype::parser::Parser::_log("exahype::parser::Parser");

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
      _tokenStream.clear();
      _interpretationErrorOccured = true;
      return;
    }
    readFile(inputFile, filename);
}

void exahype::parser::Parser::readFile(std::istream& inputFile, std::string filename) {
   try {
    const int MAX_CHARS_PER_LINE = 65536;

    std::regex COMMENT_BEGIN(R"((\/\*))"); // Covers all cases /*,/**,/***,... .
    std::regex COMMENT_END(R"((\*\/))"); //
    std::regex GROUP_BEGIN_OR_END(R"(^(\s|\t)*([a-zA-Z][^\=]+)+)");
    std::regex CONST_PARAMETER(R"(^(\s|\t)*([A-Za-z](\w|[^a-zA-Z\d\s\t])*)(\s|\t)+const(\s|\t)*=(\s|\t)*(([^\s\t]|\,\s*)+)(\s|\t)*$)");
    std::regex PARAMETER(R"(^(\s|\t)*([A-Za-z](\w|[^a-zA-Z\d\s\t])*)(\s|\t)*=(\s|\t)*(([^\s\t]|\,\s*)+)(\s|\t)*$)");
    std::regex NO_SPLITTING(R"(\}|\{)");
    std::regex COMMA_SEPARATED(R"((\w|[^a-zA-Z\,\s\t])+)");
    std::regex WHITESPACE_SEPARATED(R"(([^\s\t]+))");
    std::smatch match;

    _tokenStream.clear();
    _filename = filename;

    int currentlyReadsComment = 0;
    int lineNumber            = 0;
    while (!inputFile.eof() && inputFile) {
      char lineBuffer[MAX_CHARS_PER_LINE];
      inputFile.getline(lineBuffer, MAX_CHARS_PER_LINE);
      std::string line(lineBuffer);

      // parse the line
      if (std::regex_search(line, match, COMMENT_BEGIN) && match.size() > 1) {
        currentlyReadsComment += 1;
      }
      if (std::regex_search(line, match, COMMENT_END) && match.size() > 1) {
        currentlyReadsComment -= 1;
      }

      // Runtime parameters
      if (currentlyReadsComment==0 && std::regex_search(line, match, PARAMETER) && match.size() > 1) {
        _tokenStream.push_back(match.str(2)); // Subgroup 2 is left-hand side (trimmed)
        std::string rightHandSide = match.str(6);

        if (!std::regex_search(rightHandSide, match, NO_SPLITTING)) {
          std::regex_iterator<std::string::iterator> regex_it ( rightHandSide.begin(), rightHandSide.end(), COMMA_SEPARATED );
          std::regex_iterator<std::string::iterator> rend;
          while(regex_it!=rend) {
            _tokenStream.push_back(regex_it->str());
            ++regex_it;
          }
        } else {
          _tokenStream.push_back(rightHandSide);
        }
      // Compile time parameters (Do not push the token const on the stream)
      } else if (currentlyReadsComment==0 && std::regex_search(line, match, CONST_PARAMETER) && match.size() > 1) {
        _tokenStream.push_back(match.str(2)); // Subgroup 2 is left-hand side (trimmed)
        std::string rightHandSide = match.str(7);

        if (!std::regex_search(rightHandSide, match, NO_SPLITTING)) {
          std::regex_iterator<std::string::iterator> regex_it ( rightHandSide.begin(), rightHandSide.end(), COMMA_SEPARATED );
          std::regex_iterator<std::string::iterator> rend;
          while(regex_it!=rend) {
            _tokenStream.push_back(regex_it->str());
            ++regex_it;
          }
        } else {
          _tokenStream.push_back(rightHandSide);
        }
      } else if (currentlyReadsComment==0 && std::regex_search(line, match, GROUP_BEGIN_OR_END) && match.size() > 1) {
        std::regex_iterator<std::string::iterator> regex_it ( line.begin(), line.end(), WHITESPACE_SEPARATED );
        std::regex_iterator<std::string::iterator> rend;
        if (regex_it->str().compare("end")!=0) { // first token should not be end
          while(regex_it!=rend) {
            _tokenStream.push_back(regex_it->str());
            ++regex_it;
          }
        } // else do nothing
      } else if (currentlyReadsComment<0) {
        logError("readFile(istream)",
             "Please remove additional multi-line comment end(s) in line '" << lineNumber << "'.");
         _interpretationErrorOccured = true;
      }
      lineNumber++;
    }

    if (currentlyReadsComment>0) {
      logError("readFile(istream)",
               "A multi-line comment was not closed after line " << lineNumber);
      _interpretationErrorOccured = true;
    }
   }
  catch (const std::regex_error& e) {
    logError("readFile(istream)", "catched exception " << e.what() );
    _interpretationErrorOccured = true;
  }

  //  For debugging purposes
  if(std::getenv("EXAHYPE_VERBOSE_PARSER")) { // runtime debugging
      std::cout << getTokenStreamAsString() << std::endl;
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
  return !_tokenStream.empty() && !_interpretationErrorOccured;
}

void exahype::parser::Parser::invalidate() {
  _interpretationErrorOccured = true;
}

std::string exahype::parser::Parser::getTokenAfter(std::string token,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         _tokenStream[currentToken] != token) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

std::string exahype::parser::Parser::getTokenAfter(std::string token0,
                                           std::string token1,
                                           int additionalTokensToSkip) const {
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         _tokenStream[currentToken] != token0) {
    currentToken++;
  }
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         _tokenStream[currentToken] != token1) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

std::string exahype::parser::Parser::getTokenAfter(std::string token0, int occurance0,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  assertion(occurance0 > 0);
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         (_tokenStream[currentToken] != token0 || occurance0 > 1)) {
    if (_tokenStream[currentToken] == token0) occurance0--;
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

std::string exahype::parser::Parser::getTokenAfter(std::string token0, int occurance0,
                                           std::string token1, int occurance1,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  assertion(occurance0 > 0);
  assertion(occurance1 > 0);
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         (_tokenStream[currentToken] != token0 || occurance0 > 1)) {
    if (_tokenStream[currentToken] == token0) occurance0--;
    currentToken++;
  }
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         (_tokenStream[currentToken] != token1 || occurance1 > 1)) {
    if (_tokenStream[currentToken] == token1) occurance1--;
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

int exahype::parser::Parser::getNumberOfThreads() const {
  assertion(isValid());
  std::string token = getTokenAfter("shared-memory", "cores");
  logDebug("getNumberOfThreads()", "found token " << token);

  int result = 0;
  try {
    result = std::stoi(token);
  }
  catch (const std::invalid_argument& ia) {}

  if (result <= 0) {
    logError("getNumberOfThreads()",
        "Invalid number of cores set or token shared-memory missing: "
        << token);
    result = 1;
    _interpretationErrorOccured = true;
  }
  return result;
}

tarch::la::Vector<DIMENSIONS, double> exahype::parser::Parser::getDomainSize() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;

  token = getTokenAfter("computational-domain", "dimension", 0);	

  int dim = -1;
  try {
    dim = std::stoi(token);
  }
  catch (const std::invalid_argument& ia) {}

  if (dim < DIMENSIONS) {
    logError("getDomainSize()",
             "dimension: value "<< token << " in specification file" <<
             " does not match -DDim"<<DIMENSIONS<<" switch in Makefile");
    _interpretationErrorOccured = true;
    return result;
  }

  try {
    token = getTokenAfter("computational-domain", "width", 0);
    result(0) = std::stod(token);
    token = getTokenAfter("computational-domain", "width", 1);
    result(1) = std::stod(token);
    #if DIMENSIONS == 3
    token = getTokenAfter("computational-domain", "width", 2);
    if (token.compare("offset")==0) {
      logError("getDomainSize()",
          "width: not enough values specified for " <<
          DIMENSIONS<< " dimensions");
      _interpretationErrorOccured = true;
      return result;
    }
    result(2) = std::stod(token);
    #endif
  } catch (const std::invalid_argument& ia) {}

  logDebug("getDomainSize()", "found size " << result);
  return result;
}

tarch::la::Vector<DIMENSIONS, double> exahype::parser::Parser::getOffset() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;

  try {
    token = getTokenAfter("computational-domain", "offset", 0);
    result(0) = std::stod(token);
    token = getTokenAfter("computational-domain", "offset", 1);
    result(1) = std::stod(token);
    #if DIMENSIONS == 3
    token = getTokenAfter("computational-domain", "offset", 2);
    if (token.compare("end-time")==0) {
      logError("getOffset()",
               "offset: not enough values specified for " <<
               DIMENSIONS<< " dimensions");
      _interpretationErrorOccured = true;
      return result;
    }
    token = getTokenAfter("computational-domain", "offset", 2);
    result(2) = std::stod(token);
    #endif
  } catch (const std::invalid_argument& ia) {}

  logDebug("getOffset()", "found offset " << result);
  return result;
}

std::string exahype::parser::Parser::getMulticorePropertiesFile() const {
  std::string result = getTokenAfter("shared-memory", "properties-file");
  logDebug("getMulticorePropertiesFile()", "found token " << result);
  return result;
}

exahype::parser::Parser::MPILoadBalancingType exahype::parser::Parser::getMPILoadBalancingType() const {
  std::string token = getTokenAfter("distributed-memory", "identifier");
  exahype::parser::Parser::MPILoadBalancingType result = MPILoadBalancingType::Static;
  if (token.compare("static_load_balancing") == 0) {
    result = MPILoadBalancingType::Static;
  } else {
    logError("getMPILoadBalancingType()",
             "Invalid distributed memory identifier " << token);
    _interpretationErrorOccured = true;
  }
  return result;
}


std::string exahype::parser::Parser::getMPIConfiguration() const {
  return getTokenAfter("distributed-memory", "configure");
}


std::string exahype::parser::Parser::getSharedMemoryConfiguration() const {
  return getTokenAfter("shared-memory", "configure");
}


double exahype::parser::Parser::getNodePoolAnsweringTimeout() const {
  const std::string token         = "max-node-pool-answering-time";
  const double      defaultResult = 1e-2;
  if (getMPIConfiguration().find( token )!=std::string::npos ) {
    const double result = static_cast<double>(exahype::parser::Parser::getValueFromPropertyString(getMPIConfiguration(),token));
    if (result<0.0) {
      logWarning( "getNodePoolAnsweringTimeout()", "token " << token << " not specified for MPI configuration so use default timeout of " << defaultResult );
      return defaultResult;
    }
    else return result;
  }
  else {
    logWarning( "getNodePoolAnsweringTimeout()", "token " << token << " not specified for MPI configuration so use default timeout of " << defaultResult );
    return defaultResult;
  }
}


int exahype::parser::Parser::getMPIBufferSize() const {
  std::string token = getTokenAfter("distributed-memory", "buffer-size");

  int result = -1;
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (result <= 0) {
    logError("getMPIBufferSize()", "Invalid MPI buffer size " << token);
    result = 64;
    _interpretationErrorOccured = true;
  }
  return result;
}

int exahype::parser::Parser::getMPITimeOut() const {
  std::string token = getTokenAfter("distributed-memory", "timeout");

  int result = -1;
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (result <= 0) {
    logError("getMPIBufferSize()", "Invalid MPI timeout value " << token);
    result = 0;
    _interpretationErrorOccured = true;
  }
  return result;
}

exahype::parser::Parser::MulticoreOracleType exahype::parser::Parser::getMulticoreOracleType()
    const {
  std::string token = getTokenAfter("shared-memory", "identifier");
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
    _interpretationErrorOccured = true;
  }
  return result;
}

double exahype::parser::Parser::getSimulationEndTime() const {
  std::string token = getTokenAfter("computational-domain", "end-time");
  logDebug("getSimulationEndTime()", "found token " << token);

  double result = -1.0;
  try {
    result = std::stod(token);
  } catch (const std::invalid_argument& ia) {}

  if (result <= 0) {
    logError("getSimulationEndTime()",
             "Invalid simulation end-time: " << token);
    result = 1.0;
    _interpretationErrorOccured = true;
  }
  return result;
}

bool exahype::parser::Parser::foundSimulationEndTime() const {
  bool found = false;
  for (auto& token : _tokenStream ) {
    if ( token.compare("end-time")==0 ) {
      found = true;
      break;
    }
  }
  return found;
}

int exahype::parser::Parser::getSimulationTimeSteps() const {
  std::string token = getTokenAfter("computational-domain", "time-steps");
  logDebug("getSimulationEndTime()", "found token " << token);

  int result = -1;
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (result < 0) {
    logError("getSimulationEndTime()",
             "Invalid simulation timestep: " << token);
    _interpretationErrorOccured = true;
  }
  return result;
}

bool exahype::parser::Parser::getFuseAlgorithmicSteps() const {
  std::string token = getTokenAfter("global-optimisation", "fuse-algorithmic-steps");
  logDebug("getFuseAlgorithmicSteps()", "found fuse-algorithmic-steps"
                                            << token);

  bool result = token.compare("on") == 0;

  if (token.compare(_noTokenFound) == 0) {
    result = false;  // default value
  } else if (token.compare("on") != 0 && token.compare("off") != 0) {
    logError("getFuseAlgorithmicSteps()",
             "fuse-algorithmic-steps is required in the "
             "global-optimisation segment and has to be either on or off: "
                 << token);
    _interpretationErrorOccured = true;
  }
  return result;
}



double exahype::parser::Parser::getFuseAlgorithmicStepsFactor() const {
  if (hasOptimisationSegment()) {
      std::string token =
          getTokenAfter("global-optimisation", "fuse-algorithmic-steps-factor");

      char* pEnd;
      double result = std::strtod(token.c_str(), &pEnd); // TODO(Dominic)
      logDebug("getFuseAlgorithmicStepsFactor()",
               "found fuse-algorithmic-steps-factor " << token);

      if (result < 0.0 || result > 1.0 || pEnd == token.c_str()) {
        logError("getFuseAlgorithmicStepsFactor()",
                 "'fuse-algorithmic-steps-factor': Value must be greater than zero "
                 "and smaller than one: "
                     << result);
        result = 0.0;
        _interpretationErrorOccured = true;
      }

      return result;
  }
  else return false;
}

bool exahype::parser::Parser::getSpawnPredictionAsBackgroundThread() const {
  std::string token = getTokenAfter("global-optimisation", "spawn-predictor-as-background-thread");
  logDebug("getSpawnPredictorAsBackgroundTask()", "found spawn-predictor-as-background-thread"
                                            << token);

  bool result = token.compare("on") == 0;

  if (token.compare(_noTokenFound) == 0) {
    result = false;  // default value
  } else if (token.compare("on") != 0 && token.compare("off") != 0) {
    logError("getSpawnPredictorAsBackgroundTask()",
             "spawn-predictor-as-background-thread is required in the "
             "global-optimisation segment and has to be either on or off: "
                 << token);
    _interpretationErrorOccured = true;
  }
  return result;
}

bool exahype::parser::Parser::getExchangeBoundaryDataInBatchedTimeSteps() const {
  std::string token = getTokenAfter(
      "global-optimisation",
      "disable-amr-if-grid-has-been-stationary-in-previous-iteration");
  if (token.compare("on") != 0 && token.compare("off") != 0) {
    logError("getExchangeBoundaryDataInBatchedTimeSteps()",
             "disable-amr-if-grid-has-been-stationary-in-previous-iteration is "
             "required in the "
             "global-optimisation segment and has to be either on or off: "
                 << token);
    _interpretationErrorOccured = true;
  }
  return token.compare("off") == 0;
}

double exahype::parser::Parser::getTimestepBatchFactor() const {
  std::string token = getTokenAfter("global-optimisation", "timestep-batch-factor");
  char* pEnd;
  double result = std::strtod(token.c_str(), &pEnd); // TODO(Dominic)
  logDebug("getFuseAlgorithmicStepsFactor()", "found timestep-batch-factor "
                                                  << token);

  if (result < 0.0 || result > 1.0 || pEnd == token.c_str()) {
    logError("getFuseAlgorithmicStepsFactor()",
             "'timestep-batch-factor': Value is required in global-optimisation "
             "section and must be greater than zero and smaller than one: "
                 << result);
    result = 0.0;
    _interpretationErrorOccured = true;
  }

  return result;
}


bool exahype::parser::Parser::hasOptimisationSegment() const {
  std::string token = getTokenAfter("global-optimisation");
  return token.compare(_noTokenFound)!=0;
}


bool exahype::parser::Parser::getSkipReductionInBatchedTimeSteps() const {
  if (hasOptimisationSegment()) {
    std::string token =
      getTokenAfter("global-optimisation", "skip-reduction-in-batched-time-steps");
    logDebug("getSkipReductionInBatchedTimeSteps()",
           "found skip-reduction-in-batched-time-steps " << token);
    if (token.compare("on") != 0 && token.compare("off") != 0) {
      logError("getSkipReductionInBatchedTimeSteps()",
             "skip-reduction-in-batched-time-steps is required in the "
             "global-optimisation segment and has to be either on or off: "
                 << token);
      _interpretationErrorOccured = true;
    }

    return token.compare("on") == 0;
  }
  else return false;
}


double exahype::parser::Parser::getDoubleCompressionFactor() const {
  std::string token = getTokenAfter("global-optimisation", "double-compression");

  if (token.compare(_noTokenFound) == 0) {
    return 0.0;  // default value
  }
  else {
    char* pEnd;
    double result = std::strtod(token.c_str(), &pEnd); // TODO(Dominic)
    logDebug("getDoubleCompressionFactor()", "found double-compression "
                                                  << token);

    if (result < 0.0 || pEnd == token.c_str()) {
      logError("getDoubleCompressionFactor()",
             "'double-compression': Value is required in global-optimisation "
             "section and must be greater than or equal to zero: " << result);
      result = 0.0;
      _interpretationErrorOccured = true;
    }

    return result;
  }
}


bool   exahype::parser::Parser::getSpawnDoubleCompressionAsBackgroundTask() const {
  std::string token =
      getTokenAfter("global-optimisation", "spawn-double-compression-as-background-thread");

  if (token.compare(_noTokenFound) == 0) {
    return false;  // default value
  }
  else {
    logDebug("getSpawnDoubleCompressionAsBackgroundTask()",
           "found spawn-double-compression-as-background-thread " << token);
    if (token.compare("on") != 0 && token.compare("off") != 0) {
      logError("getSpawnDoubleCompressionAsBackgroundTask()",
             "spawn-double-compression-as-background-thread is required in the "
             "global-optimisation segment and has to be either on or off: "
                 << token);
      _interpretationErrorOccured = true;
    }

    return token.compare("on") == 0;
  }
}


exahype::solvers::Solver::Type exahype::parser::Parser::getType(
    int solverNumber) const {
  std::string token;
  exahype::solvers::Solver::Type result =
      exahype::solvers::Solver::Type::ADERDG;
  token = getTokenAfter("solver", solverNumber + 1, 0);
  if (_identifier2Type.find(token) != _identifier2Type.end()) {
    result = _identifier2Type.at(token);
    // logDebug("getType()", "found type " << result);
    logDebug("getType()", "found type ");
  } else {
    logError(
        "getType()",
        "'" << getIdentifier(solverNumber) << "': 'type': Value '" << token
            << "' is invalid. See the ExaHyPE documentation for valid values.");
    _interpretationErrorOccured = true;
  }
  return result;
}

std::string exahype::parser::Parser::getIdentifier(int solverNumber) const {
  std::string token;
  token = getTokenAfter("solver", solverNumber + 1, 1);
  logDebug("getIdentifier()", "found identifier " << token);
  return token;
}

int exahype::parser::Parser::getVariables(int solverNumber) const {
  std::string token;
  std::regex COLON_SEPARATED(R"(([A-Za-z]\w*):([0-9]+))");
  std::smatch match;

  // first check if we read in a number
  int result = -1;
  token = getTokenAfter("solver", solverNumber + 1, "variables", 1);
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (result < 1) { // token is not a number
    result = 0;
    int i = 1;
    std::regex_search(token, match, COLON_SEPARATED);
    while (match.size() > 1) {
      int multiplicity = atoi(match.str(2).c_str());
      result +=multiplicity; // std::string name = match.str(1);

      // logInfo("getVariables(...)","token="<<token<<",name="<<match.str(1)<<",n="<<match.str(2));
      token = getTokenAfter("solver", solverNumber + 1, "variables", 1, i++);
      std::regex_search(token, match, COLON_SEPARATED);
    }

    if (result < 1) { // token is still 0
      logError("getVariables()",
               "'" << getIdentifier(solverNumber)
               << "': 'variables': Value must be greater than zero.");
          _interpretationErrorOccured = true;
    }
  }

  logDebug("getVariables()", "found variables " << result);
  return result;
}

int exahype::parser::Parser::getParameters(int solverNumber) const {
  std::string token;
  std::regex COLON_SEPARATED(R"(([A-Za-z]\w*):([0-9]+))");
  std::smatch match;

  // first check if we read in a number
  int result = -1;
  token = getTokenAfter("solver", solverNumber + 1, "parameters", 1);
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (result < 1) { // token is not a number
    result = 0;

    int i = 1;
    std::regex_search(token, match, COLON_SEPARATED);
    while (match.size() > 1) {
      int multiplicity = atoi(match.str(2).c_str());
      result +=multiplicity; // std::string name = match.str(1);

      // logInfo("getVariables(...)","token="<<token<<",name="<<match.str(1)<<",n="<<match.str(2));
      token = getTokenAfter("solver", solverNumber + 1, "parameters", 1, i++);
      std::regex_search(token, match, COLON_SEPARATED);
    }

    if (result < 0) { // token is still 0
      logError("getParameters()",
               "'" << getIdentifier(solverNumber)
               << "': 'parameters': Value must be non-negative.");
          _interpretationErrorOccured = true;
    }
  }

  logDebug("getVariables()", "found variables " << result);
  return result;
}

int exahype::parser::Parser::getOrder(int solverNumber) const {
  std::string token;

  int result;
  token = getTokenAfter("solver", solverNumber + 1, "order", 1);
  try {
    result = std::stoi(token);
  }
  catch (const std::invalid_argument& ia) {
    logError("getOrder()", "'" << getIdentifier(solverNumber)
        << "': 'order': Value must not be negative.");
    _interpretationErrorOccured = true;
    result = -1;
  }

  logDebug("getOrder()", "found order " << result);
  return result;
}


double exahype::parser::Parser::getMaximumMeshSize(int solverNumber) const {
  std::string token;

  double result = -1.0;
  token = getTokenAfter("solver", solverNumber + 1, "maximum-mesh-size", 1, 0);
  try {
    result = std::stod(token);
  } catch (const std::invalid_argument& ia) {}

  if (tarch::la::smallerEquals(result, 0.0)) {
    logError("getMaximumMeshSize(int)",
             "'" << getIdentifier(solverNumber)
                 << "': 'maximum-mesh-size': Value must be greater than zero.");
    _interpretationErrorOccured = true;
  }

  logDebug("getMaximumMeshSize()", "found maximum mesh size " << result);
  return result;
}

double exahype::parser::Parser::getMaximumMeshDepth(int solverNumber) const {
  std::string token;

  int result = 0;
  token = getTokenAfter("solver", solverNumber + 1, "maximum-mesh-depth", 1, 0);
  if (token==_noTokenFound) {
    return result;
  }

  result = -1;
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (tarch::la::smaller(result, 0)) {
    logError("getMaximumMeshDepth(int)",
             "'" << getIdentifier(solverNumber)
                 << "': 'maximum-mesh-depth': Value must be greater than or equal to zero.");
    _interpretationErrorOccured = true;
  }

  logDebug("getMaximumMeshDepth()", "found maximum mesh size " << result);
  return result;
}

exahype::solvers::Solver::TimeStepping exahype::parser::Parser::getTimeStepping(
    int solverNumber) const {
  std::string token;
  exahype::solvers::Solver::TimeStepping result;
  token = getTokenAfter("solver", solverNumber + 1, "time-stepping", 1);
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
    _interpretationErrorOccured = true;
  }
  return exahype::solvers::Solver::TimeStepping::Global;
}

double exahype::parser::Parser::getDMPRelaxationParameter(int solverNumber) const {
  std::string token;

  double result = -1.0;
  token = getTokenAfter("solver", solverNumber + 1, "dmp-relaxation-parameter", 1);
  try {
    result = std::stod(token);
  } catch (const std::invalid_argument& ia) {}

  if (result < 0) {
    logError("getDMPRelaxationParameter()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-relaxation-parameter': Value must not be negative.");
    _interpretationErrorOccured = true;
  }

  logInfo("getParameters()", "found dmp-relaxation-parameter " << result);
  return result;
}

double exahype::parser::Parser::getDMPDifferenceScaling(int solverNumber) const {
  std::string token;

  double result = -1.0;
  token = getTokenAfter("solver", solverNumber + 1, "dmp-difference-scaling", 1);
  try {
    result = std::stod(token);
  } catch (const std::invalid_argument& ia) {}

  if (result < 0) {
    logError("getDMPDifferenceScaling()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-difference-scaling': Value must not be negative.");
    _interpretationErrorOccured = true;
  }

  logInfo("getDMPDifferenceScaling()", "found dmp-difference-scaling " << result);
  return result;
}

int exahype::parser::Parser::getDMPObservables(int solverNumber) const {
  std::string token;

  int result = -1; // is required
  token = getTokenAfter("solver", solverNumber + 1, "dmp-observables", 1);
  try {
    result = std::stoi(token);
  } catch (const std::invalid_argument& ia) {}

  if (result < 0) {
    logError("getDMPObservables()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-observables': Value must not be negative.");
    _interpretationErrorOccured = true;
  }

  logInfo("getDMPObservables()", "found dmp-observables " << result);
  return result;
}

int exahype::parser::Parser::getStepsTillCured(int solverNumber) const {
  std::string token;

  int result = 0; // default value
  token = getTokenAfter("solver", solverNumber + 1, "steps-till-cured", 1);
  if (token.compare(_noTokenFound)!=0) {
    try {
      result = std::stoi(token);
    } catch (const std::invalid_argument& ia) {
      result = -1;
    }

    if (result < 0) {
      logError("getStepsTillCured()",
               "'" << getIdentifier(solverNumber)
               << "': 'steps-till-cured': Value must be integral and not negative.");
      _interpretationErrorOccured = true;
    }

    logInfo("getStepsTillCured()", "found steps-till-cured " << result);
  }
  return result;
}

int exahype::parser::Parser::getLimiterHelperLayers(int solverNumber) const {
  std::string token;
  int result = 1; // default value
  token = getTokenAfter("solver", solverNumber + 1, "helper-layers", 1);

  if (token.compare(_noTokenFound)!=0) {
    try {
      result = std::stoi(token);
    } catch (const std::invalid_argument& ia) {
      result = -1;
    }

    if (result < 1) {
      logError("getLimiterHelperLayers()",
               "'" << getIdentifier(solverNumber)
               << "': 'helper-layers': Value must be integral and greater or equal to 1.");
      _interpretationErrorOccured = true;
    }
    logInfo("getLimiterHelperLayers()", "found helper-layers " << result);
  }
  return result;
}

std::string exahype::parser::Parser::getIdentifierForPlotter(int solverNumber,
                                                     int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1);
  logDebug("getIdentifierForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return token;
}

std::string exahype::parser::Parser::getNameForPlotter(int solverNumber,
                                               int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 1);
  logDebug("getIdentifierForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return token;
}

int exahype::parser::Parser::getUnknownsForPlotter(int solverNumber,
                                           int plotterNumber) const {
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 3);
  logDebug("getUnknownsForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  try {
    return std::stoi(token);
  } catch (const std::invalid_argument& ia) {
    return 0;
  }
}

double exahype::parser::Parser::getFirstSnapshotTimeForPlotter(
    int solverNumber, int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 5);
  logDebug("getFirstSnapshotTimeForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);

  try {
    return std::stod(token);
  } catch (const std::invalid_argument& ia) {
    logError("getFirstSnapshotTimeForPlotter()",
        "'" << getIdentifier(solverNumber)
        << "' - plotter "<<plotterNumber<<": 'time' value must be a float.");
    _interpretationErrorOccured = true;
    return std::numeric_limits<double>::max();
  }
}

double exahype::parser::Parser::getRepeatTimeForPlotter(int solverNumber,
                                                int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 7);
  logDebug("getRepeatTimeForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);

  try {
    return std::stod(token);
  } catch (const std::invalid_argument& ia) {
    logError("getRepeatTimeForPlotter()",
        "'" << getIdentifier(solverNumber)
        << "' - plotter "<<plotterNumber<<": 'repeat' value must be a float.");
    _interpretationErrorOccured = true;
    return std::numeric_limits<double>::max();
  }
}

std::string exahype::parser::Parser::getFilenameForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 9);
  logDebug("getFilenameForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return token;
}

std::string exahype::parser::Parser::getSelectorForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 11);
  token += "," + getTokenAfter("solver", solverNumber + 1, "plot",
                                      plotterNumber + 1, 12);
  #if DIMENSIONS==3
  token += "," + getTokenAfter("solver",solverNumber + 1, "plot",
                                 plotterNumber + 1, 13);
  #endif

  logDebug("getSelectorForPlotter()", "found token " << token);
  return (token != _noTokenFound) ? token : "";
}

std::string exahype::parser::Parser::getLogFileName() const {
  std::string token = getTokenAfter("log-file");
  logDebug("getLogFileName()", "found token " << token);
  return (token != _noTokenFound) ? token : "";
}

std::string exahype::parser::Parser::getProfilerIdentifier() const {
  std::string token = getTokenAfter("profiling", "profiler");
  logDebug("getProfilerIdentifier()", "found token " << token);
  return (token != _noTokenFound) ? token : "NoOpProfiler";
}

std::string exahype::parser::Parser::getMetricsIdentifierList() const {
  std::string token = getTokenAfter("profiling", "metrics");
  logDebug("getMetricsIdentifierList()", "found token " << token);
  return (token != _noTokenFound) ? token : "{}";
}

std::string exahype::parser::Parser::getProfilingOutputFilename() const {
  std::string token = getTokenAfter("profiling", "profiling-output");
  logDebug("getProfilingOutputFilename()", "found token " << token);
  return (token != _noTokenFound) ? token : "";
}

void exahype::parser::Parser::logSolverDetails(int solverNumber) const {
  logInfo("logSolverDetails()",
          "Solver " << getTokenAfter("solver", solverNumber + 1, 0) << " "
                    << getIdentifier(solverNumber) << ":");
  logInfo("logSolverDetails()", "variables:\t\t" << getVariables(solverNumber));
  logInfo("logSolverDetails()", "parameters:\t" << getParameters(solverNumber));
  logInfo("logSolverDetails()", "order:\t\t" << getOrder(solverNumber));
  logInfo("logSolverDetails()", "maximum-mesh-size:\t"
                                    << getMaximumMeshSize(solverNumber));
  logInfo("logSolverDetails()",
          "time-stepping:\t" << getTokenAfter("solver", solverNumber + 1,
                                              "time-stepping", 1));
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
    _interpretationErrorOccured = true;
  }

  if (solver->getIdentifier().compare(getIdentifier(solverNumber))) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Identifier in specification file "
                 << "('" << getIdentifier(solverNumber)
                 << "') differs from identifier used in implementation ('"
                 << solver->getIdentifier() << "').");
    recompile = true;
    _interpretationErrorOccured = true;
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
    _interpretationErrorOccured = true;
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
    _interpretationErrorOccured = true;
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
    _interpretationErrorOccured = true;
  }

  // @todo We should add checks for FV as well
  
  // (Sven:) somehow for me the following lines are never printed. I don't
  // know why.

  if (runToolkitAgain) {
    logError("checkSolverConsistency",
             "Please (1) run the Toolkit again, and (2) recompile!");
    _interpretationErrorOccured = true;
  }

  if (recompile) {
    logError(
        "checkSolverConsistency",
        "Please (1) adjust the specification file (*.exahype) or the file '"
            << solver->getIdentifier()
            << ".cpp' accordingly, and (2) recompile!");
    _interpretationErrorOccured = true;
  }
}

std::string exahype::parser::Parser::getSpecfileName() const {
  return _filename;
}

std::string exahype::parser::Parser::getTokenStreamAsString() const {
  std::stringstream ret;
  ret << "Parser _tokenStream=" << std::endl;
  for (std::string str : _tokenStream) {
    ret << "["<<str<<"]" << std::endl;
  }
  return ret.str();
}

int exahype::parser::Parser::getRanksPerNode() {
  const std::string RanksPerNode = "ranks-per-node";

  return static_cast<int>(exahype::parser::Parser::getValueFromPropertyString(getMPIConfiguration(),RanksPerNode));
}


int exahype::parser::Parser::getNumberOfBackgroundTasks() {
  const std::string Search = "background-tasks";

  int result = static_cast<int>(exahype::parser::Parser::getValueFromPropertyString(getSharedMemoryConfiguration(),Search));
  if (result<-2) {
    logWarning("getNumberOfBackgroundTasks()", "invalid number of background tasks (background-tasks field in configuration) " <<
      "set or no number at all. Use default (1). See BackgroundTasks.h for documentation.");
    result = 1;
  }
  return result;
}


bool exahype::parser::Parser::useManualPinning() {
  const std::string Search = "manual-pinning";

  return getSharedMemoryConfiguration().find(Search) != std::string::npos;
}

exahype::parser::ParserView exahype::parser::Parser::createParserView(const int solverNumberInSpecificationFile) {
  return exahype::parser::ParserView(*this,solverNumberInSpecificationFile);
}