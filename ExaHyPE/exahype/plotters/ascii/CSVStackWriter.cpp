#include "CSVStackWriter.h"
#include <exception>

void exahype::plotters::ascii::CSVStackWriter::writeHeader() {
	// Writes a header block which goes at the top of the file
	writeCommentLine("exahype::plotters::ascii::CSVStackWriter ASCII output");
	#ifdef PARALLEL
	os << commentIntro <<
	"MPI Rank "
	tarch::parallel::Node::getInstance()::getRank()
	<< " of " <<
	tarch::parallel::Node::getInstance()::getNumberOfNodes()
	<< "total ranks" << newline;
	#else
	writeCommentLine("No MPI Setup -- single rank here.");
	#endif
	os << commentIntro << "Created on " << str_time() << newline;
	
	// Human readable column list including description, in a format
	// suitable for gnuplotting (kind of similar to CarpetASCII)
	writeCommentLine(""); int c=1;
	for(auto& writer : writers) {
		for(auto&& col : writer.columns) {
			os << commentIntro
			<< c++ << ":" << col.name
			<< " -- "
			<< col.description
			<< newline;
		}
		os << newline;
	}
	writeCommentLine("");
	
	// Columns list
	for(auto& writer : writers) {
		if(!rawcolumns) os << commentIntro;
		auto &last = *(--writer.columns.end());
		for(auto&& col : writer.columns) {
			os << col.name;
			if(&col != &last) os << seperator;
		}
		if(!rawcolumns) os << newline;
	}
	os << newline;
}

void exahype::plotters::ascii::CSVStackWriter::writeColumns(voidptr line) {
	auto& writer = writers[numberOfWrittenCSVWritersInCurrentLine++];
	writer.newline = ""; // to be safe, do it every time
	writer.writeRow(line);
}


void exahype::plotters::ascii::CSVStackWriter::finishRow() {
	numberOfWrittenCSVWritersInCurrentLine = 0;
        os << newline;
	os.flush(); // important
}

void exahype::plotters::ascii::CSVStackWriter::writeRow(voidptr line) {
	throw std::runtime_error("You're using CSVStackWriter wrong. Use writeColumns() instead.");
}
