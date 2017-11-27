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
 
#ifndef __EXAHYPE_PLOTTERS_ASCII_CSV_LIB__
#define __EXAHYPE_PLOTTERS_ASCII_CSV_LIB__

#include <fstream> 

namespace exahype {
  namespace plotters {
    class ASCIIWriter;

  }
}
/*
class exahype::plotters::ASCIIWriter {
public:
  std::string   seperator;
  std::string   filename;
  std::ofstream ofs;

  struct DataLine;

  ASCIIWriter(const std::string& filename, const std::string& seperator);
  virtual ~ASCIIWriter();

  void writeCommentLine(const std::string& line);
  void startHeaderLine();
  void writeHeader
  void finishHeaderLine();

  void startDataLine();

  void finishDataLine();

};

class exahype::plotters::ASCIIWriter::DataLine {
	const exahype::plotters::ASCIIWriter& ref;
	template<typename T>
	void write(T stuff) { ref.ofs << stuff; }
};
*/

// To be done: A generalization for ASCII writers.

#endif /* __EXAHYPE_PLOTTERS_ASCII_CSV_LIB__ */
