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

#ifndef EXAHYPE_PARSER_ParserView_H
#define EXAHYPE_PARSER_ParserView_H

namespace exahype {
namespace parser {
class Parser;
class ParserView;
}
}

#include <iostream>

#include <map>
#include <vector>
#include <utility> // pair
#include <istream>
#include <ostream>

#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

/**
 * View on the parser
 *
 * An instance of this class is a parser view. While the parser sees the
 * whole specification file, a view only 'sees' the fragment that is
 * specific to one solver. As such, we do pass it to solvers (that hold
 * constants) and then allow these solvers to read their data from the
 * config file.
 *
 * From the user's point of view, this class provides access to key-value
 * pairs. If you have an instruction alike
 * <pre>
 *   constants         = {rho:0.4567,gamma:-4,alpha:4.04e-5,file:output}
 * </pre>
 * in your
 *
 * @author Tobias Weinzierl
 */
class exahype::parser::ParserView {
private:
  static tarch::logging::Log _log;
  exahype::parser::Parser& _parser;
  int _solverNumberInSpecificationFile;

  /**
   * @return Value for given key. Returns the empty string if there is no
   *         value. Returns empty string "" if the key does not exist.
   */
  std::string getValue(const std::string& key) const;

public:
  ParserView(exahype::parser::Parser& parser, int solverNumberInSpecificationFile);

  /**
   * You may use keys without a value. This operation allows you to check
   * whether there are such keys. Furthermore, you might use this guy as
   * a preamble to the other getters.
   */
  bool hasKey(const std::string& key) const;

  /**
   * Please ensure that isValueValidXXX holds before you invoke this
   * operation.
   */
  bool getValueAsBool(const std::string& key) const;

  /**
   * Please ensure that isValueValidXXX holds before you invoke this
   * operation.
   */
  int getValueAsInt(const std::string& key) const;

  /**
   * Please ensure that isValueValidXXX holds before you invoke this
   * operation.
   */
  double getValueAsDouble(const std::string& key) const;

  /**
   * Please ensure that isValueValidXXX holds before you invoke this
   * operation.
   */
  std::string getValueAsString(const std::string& key) const;

  /**
   * Returns the contents of this ParserView as a map, i.e. the keys
   * mapped to the values which remain as strings (uncasted to their
   * native description, i.e. bool,int,double).
   *
   * @see getAllAsOrderedMap()
   **/
  std::map<std::string, std::string> getAllAsMap() const;

  /**
   * Returns the contents of this ParserView as an "ordered map", i.e
   * a vector of key -> value pairs.
   *
   * @see getAllAsMap()
   **/
  std::vector< std::pair<std::string, std::string> > getAllAsOrderedMap() const;

  bool isValueValidBool(const std::string& key) const;
  bool isValueValidInt(const std::string& key) const;
  bool isValueValidDouble(const std::string& key) const;
  bool isValueValidString(const std::string& key) const;

  /**
   * Access the Parser. Can be helpful for obtaining information beyond the
   * solver constants which are represented by this object.
   **/
  const exahype::parser::Parser& getParser() const;

  /**
   * Returns a short string representation of this ParserView. It mentions the
   * Parser Filename as well as the ParserView's solver name.
   **/
  std::string toString() const;

  /**
   * \see toString()
   **/
  void toString(std::ostream& out) const;
};

#endif
