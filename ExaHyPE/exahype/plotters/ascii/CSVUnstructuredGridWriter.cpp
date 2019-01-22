#include "CSVUnstructuredGridWriter.h"

using VertexDataWriter = exahype::plotters::ascii::CSVUnstructuredGridWriter::VertexDataWriter;
using CellDataWriter   = exahype::plotters::ascii::CSVUnstructuredGridWriter::CellDataWriter;
using VertexWriter     = exahype::plotters::ascii::CSVUnstructuredGridWriter::VertexWriter;
using CellWriter       = exahype::plotters::ascii::CSVUnstructuredGridWriter::CellWriter;

VertexDataWriter* exahype::plotters::ascii::CSVUnstructuredGridWriter::createVertexDataWriter( const std::string& identifier, int recordsPerVertex ) {
	VertexDataWriter* v = new VertexDataWriter(*this, identifier, recordsPerVertex);
	return v;
}

/*

	CellDataWriter*   createCellDataWriter( const std::string& identifier, int recordsPerCell ) override
		{ return new CellDataWriter(*this, identifier, recordsPerCell); }

	VertexWriter*   createVertexWriter() override { return new VertexWriter(); }
	CellWriter*     createCellWriter() override { return new CellWriter(); }
	
*/
