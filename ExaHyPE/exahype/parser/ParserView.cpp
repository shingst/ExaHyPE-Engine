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
#include "exahype/parser/ParserView.h"

#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>
#include <string>
#include <regex>
#include <cstdlib> // getenv, exit
#include <sstream>

#include "tarch/la/ScalarOperations.h"

#include "exahype/parser/Parser.h"

tarch::logging::Log exahype::parser::ParserView::_log( "exahype::parser::ParserView" );

exahype::parser::ParserView::ParserView(exahype::parser::Parser& parser,
                                        int solverNumberInSpecificationFile)
    : _parser(parser),
      _solverNumberInSpecificationFile(solverNumberInSpecificationFile) {}

std::string exahype::parser::ParserView::getValue(const std::string& key) const {
  assertion(_parser.isValid());

  std::string token;
  std::regex  COLON_SEPARATED(R"((.+):(.+))");
  std::smatch match;

  token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, 0);
  std::regex_search(token, match, COLON_SEPARATED);

  int i = 1;
  while (match.size() > 1) {
    if (match.str(1).compare(key)==0) {
      logDebug("getValue()", "solver " << _solverNumberInSpecificationFile + 1 << ": found constant '" << key << "' with value '" << match.str(2) << "'.");
      return match.str(2);
    }
    token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, i++);
    std::regex_search(token, match, COLON_SEPARATED);
  }
  logDebug("hasKey()", "solver " << _solverNumberInSpecificationFile + 1 << ": cannot find constant '" << key << "'.");
  return "";
}

bool exahype::parser::ParserView::hasKey(const std::string& key) const {
  assertion(_parser.isValid());

  std::string token;
  std::regex  COLON_SEPARATED(R"(.+):(.+))");
  std::smatch match;

  token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, 0);
  std::regex_search(token, match, COLON_SEPARATED);

  int i = 1;
  while (match.size() > 1) {
    if (match.str(1).compare(key)) {
      logDebug("hasKey()", "solver " << _solverNumberInSpecificationFile + 1 << ": found constant '" << key << "' with value '" << match.str(2) << "'.");
      return true;
    }
    token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, i++);
    std::regex_search(token, match, COLON_SEPARATED);
  }
  logDebug("hasKey()", "solver " << _solverNumberInSpecificationFile + 1 << ": cannot find constant '" << key << "'.");
  return false;
}

int exahype::parser::ParserView::getValueAsInt(const std::string& key) const {
  std::string value = getValue(key);

  int result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return result;
  } else {
    assertion(!isValueValidInt(key));
    assertionMsg(false, "shall not happen. Please call isValueValidXXX before");
    return -1;
  }
}

bool exahype::parser::ParserView::getValueAsBool(const std::string& key) const {
  std::string value = getValue(key);

  // We use 'on' and 'off' for multiple switches in the specification file
  // which would lead the java parser to object
  if (value.compare("true") == 0 ||
      value.compare("1") == 0) {
    return true;
  }
  else if (value.compare("false") == 0 ||
           value.compare("0") == 0) {
    return false;
  } else {
    assertion(!isValueValidBool(key));
    assertionMsg(false, "shall not happen. Please call isValueValidXXX before");
    return false;
  }
}

double exahype::parser::ParserView::getValueAsDouble(
    const std::string& key) const {
  std::string value = getValue(key);

  double result;
  std::istringstream ss(value); // TODO(Dominic)
  ss >> result;

  if (ss) {
    return result;
  } else {
    assertion(!isValueValidDouble(key));
    assertionMsg(false, "shall not happen. Please call isValueValidXXX before");
    return -1.0;
  }
}

std::string exahype::parser::ParserView::getValueAsString(
    const std::string& key) const {
  return getValue(key);
}

bool exahype::parser::ParserView::isValueValidInt(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(key);

  int result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    return false;
  }
}

bool exahype::parser::ParserView::isValueValidDouble(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(key);

  double result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    return false;
  }
}


bool exahype::parser::ParserView::isValueValidBool(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(key);

  // We use 'on' and 'off' for multiple switches in the specification file
  if (value.compare("true")  == 0 ||
      value.compare("false") == 0 ||
      value.compare("1") == 0 ||
      value.compare("0") == 0) {
    return true;
  } else {
    return false;
  }
}

std::vector< std::pair<std::string, std::string> > exahype::parser::ParserView::getAllAsOrderedMap() const {
  std::vector< std::pair<std::string, std::string> > retvec;

  std::string token;
  std::regex  COLON_SEPARATED(R"(^(.+?):(.+?),?$)");
  std::smatch match;

  token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, 0);
  std::regex_search(token, match, COLON_SEPARATED);

  int i = 1;
  while (match.size() > 1) {
    retvec.push_back( make_pair(match.str(1), match.str(2)) );
    token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, i++);
    std::regex_search(token, match, COLON_SEPARATED);
  }
  return retvec;
}

std::map<std::string, std::string> exahype::parser::ParserView::getAllAsMap() const {
  std::map<std::string, std::string> retmap;

  std::string token;
  std::regex  COLON_SEPARATED(R"(^(.+?):(.+?),?$)");
  std::smatch match;

  token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, 0);
  std::regex_search(token, match, COLON_SEPARATED);

  int i = 1;
  while (match.size() > 1) {
    retmap[match.str(1)] = match.str(2);
    token = _parser.getTokenAfter("solver", _solverNumberInSpecificationFile + 1, "constants", 1, i++);
    std::regex_search(token, match, COLON_SEPARATED);
  }
  return retmap;
}

bool exahype::parser::ParserView::isValueValidString(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  return getValue(key) != "";
}

const exahype::parser::Parser& exahype::parser::ParserView::getParser() const {
  return _parser;
}


std::string exahype::parser::ParserView::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::parser::ParserView::toString(std::ostream& out) const {
  out << "ParserView(";
  out << "specfile:" << _parser.getSpecfileName();
  out << ",";
  out << "solver:" << _parser.getIdentifier(_solverNumberInSpecificationFile);
  out <<  ")";
}

