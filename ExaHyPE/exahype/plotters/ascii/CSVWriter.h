#ifndef __EXAHYPE_PLOTTERS_ASCII__CSVWRITER__
#define __EXAHYPE_PLOTTERS_ASCII__CSVWRITER__

namespace exahype {
  namespace plotters {
    namespace ascii {
      class CSVWriter;
    }
  }
}

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>

/**
 * CSVWriter is a semi-reflecting writer for ASCII/CSV files (most likely used in time series
 * "plots" in ExaHyPE).
 * 
 * The central idea is that you have a structure in your code/plotter, such as
 *
 *   struct Line { int time; double x, y, z; bool active_flag; }
 *
 * and CSVWriter provides you the library code to write out this structure every now and then
 * without much registration overhead. This is archived by storing relative positions of the
 * data fields in the structures (C macro offsetof). I call this "semi-reflection".
 *
 * What CSVWriter does not:
 *  - Generation of machine readable self-describing data (this should be part of another layer ontop)
 *  - Global reductions or MPI communication (CSVWriter does not know much about MPI)
 **/
struct exahype::plotters::ascii::CSVWriter {
	using voidptr = char*; ///< A pointer to a 1-byte-width thing

	struct Column {
		std::string name;   ///< A short name for a CSV column such as "l1norm"
		std::string format; ///< A formatstring
		std::string description; ///< A one line human-readable description such as "Defines the l1Norm for a certain quantity"

		enum class Type { INT, DOUBLE, BOOL, STRING, UNDEF } type;
		uintptr_t offset; ///< Storage position in structure where it is read off
		std::string get(voidptr lineStruct); ///< get value formatted as string, given a structure
	};

	std::vector<Column> columns; ///< Column specifications, go set them!
	std::ofstream os; ///< The Stream you want to write to	
	
	std::string seperator = "\t"; ///< Column seperator
	std::string newline = "\n"; ///< Newline character to use
	std::string commentIntro = "# "; ///< Characters how to introduce a comment line
	bool rawcolumns = true; ///< Display the column names without a comment (=machine readable)
	
	void writeCommentLine(const std::string& line);
	void writeHeader(); ///< Write a header block
	
	/**
	 * Write an actual row. Here you should pass the pointer to your structure.
	 **/
	void writeRow(voidptr line);
	
	// Two helper routines if you want to delegate the file handling to this class:
	
	/// Compose a typical ExaHyPE output filename
	static std::string combineFilename(std::string base, std::string suffix=".asc");
	void openFile(std::string base, std::string suffix=".asc");

};


// Syntactic sugar for standard stuff

#include <cstddef> // offsetof macro (not member of std::)

#define CSVWRITER_INTEGER_COLUMN(userStructure, fieldName, description) \
	{ #fieldName, "%d", description, exahype::plotters::ascii::CSVWriter::Column::Type::INT, offsetof(userStructure, fieldName) }
		
#define CSVWRITER_DOUBLE_COLUMN(userStructure, fieldName, description) \
	{ #fieldName, "%e", description, exahype::plotters::ascii::CSVWriter::Column::Type::DOUBLE, offsetof(userStructure, fieldName) }
		
		
#define CSVWRITER_WRITE_ROW(csvWriterInstance, userStructureInstance) \
	(csvWriterInstance).writeRow( (exahype::plotters::ascii::CSVWriter::voidptr) &userStructureInstance )

// etc.

#endif /* __EXAHYPE_PLOTTERS_ASCII__CSVWRITER__ */
