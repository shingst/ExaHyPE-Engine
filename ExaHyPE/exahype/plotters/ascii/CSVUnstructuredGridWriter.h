#ifndef __EXAHYPE_PLOTTERS_ASCII_CSV_UNSTRUCTURED_GRID_WRITER_SVEN__
#define __EXAHYPE_PLOTTERS_ASCII_CSV_UNSTRUCTURED_GRID_WRITER_SVEN__

namespace exahype {
namespace plotters {
namespace ascii {
	class CSVUnstructuredGridWriter;
}
}
}

#include "tarch/plotter/griddata/unstructured/UnstructuredGridWriter.h"

/**
 * An CSV/ASCII interface to Peanos UnstructuredGridWriter interface,
 * suitable for use in any plotters used here.
 * 
 * This class should move to tarch/plotter/griddata/unstructured/csv/CSVTextFileWriter once,
 * i.e. being part of Peano. Feel free to do so.
 * 
 * @author SvenK
 **/
class exahype::plotters::ascii::CSVUnstructuredGridWriter:
	public tarch::plotter::griddata::unstructured::UnstructuredGridWriter {

	public:
	
	// @return Write has been successful
	bool writeToFile( const std::string& filename ) override;

	// @return Whether writer is ready to accept data.
	bool isOpen() override;

	/**
	 * Clear the writer, i.e. erase all the data. However, as the writer does
	 * not track how many vertex and cell writers you've created, it's up to
	 * you to ensure that none of these instances is left.
	 * 
	 * -> in fact, mine will do
	 */
	void clear() override;
	
	// A writer for scalar data on points (vertices).
	struct VertexDataWriter: 
		public tarch::plotter::griddata::Writer::VertexDataWriter {
			
		CSVUnstructuredGridWriter &parent;
		const std::string identifier;
		const int recordsPerVertex;
		
		VertexDataWriter(CSVUnstructuredGridWriter& parent, const std::string& identifier, int recordsPerVertex)
			: parent(parent), identifier(identifier), recordsPerVertex(recordsPerVertex) {}
		
		void plotVertex( int index, double* values, int numberOfValues ) override;
		
		void plotVertex( int index, double value ) override { plotVertex(index, &value, 1); }
		void plotVertex( int index, const tarch::la::Vector<2,double>& value ) override { plotVertex(index, (double*)value.data(), 2); }
		void plotVertex( int index, const tarch::la::Vector<3,double>& value ) override { plotVertex(index, (double*)value.data(), 3); }

		void close() override;
		void assignRemainingVerticesDefaultValues() override { close(); }
	};
	
	// A writer for scalar data on elements.
	struct CellDataWriter:
		public tarch::plotter::griddata::Writer::CellDataWriter {
			
		CSVUnstructuredGridWriter &parent;
		const std::string identifier;
		const int recordsPerCell;
		
		CellDataWriter(CSVUnstructuredGridWriter& parent, const std::string& identifier, int recordsPerCell)
			: parent(parent), identifier(identifier), recordsPerCell(recordsPerCell) {}

		void plotCell(int index, double* values, int numberOfValues ) override;
		
		void plotCell(int index, double value ) override { plotCell(index, &value, 1); }
		void plotCell(int index, const tarch::la::Vector<2,double>& value ) override { plotCell(index, (double*)value.data(), 2); }
		void plotCell(int index, const tarch::la::Vector<3,double>& value ) override { plotCell(index, (double*)value.data(), 3); }

		void close() override;
		void assignRemainingCellsDefaultValues() override { close(); }
	};
	
	struct VertexWriter:
		public tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter {
		int dimension = 0;

		int plotVertex(const tarch::la::Vector<2,double>& position) override;
		int plotVertex(const tarch::la::Vector<3,double>& position) override;

		void close() override;
	};

	struct CellWriter:
		public tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter {
		int vertexIndexLength = 0;

		int plot(int* vertexIndices, int size);
			
		int plotPoint(int vertexIndex) override { return plot(&vertexIndex, 1); }
		int plotLine(int vertexIndex[2]) override { return plot(vertexIndex, 2); }
		int plotTriangle(int vertexIndex[3]) override { return plot(vertexIndex, 3); }
		int plotQuadrangle(int vertexIndex[4]) override { return plot(vertexIndex, 4); }
		int plotHexahedron(int vertexIndex[8]) override { return plot(vertexIndex, 8); }

		void close() override;
	};
	
	// Vertex/Cell managament, from the UnstructuredGridWriter.
	int vertexIndex, cellIndex;

	// whatver... this is useless.
	//std::list<VertexDataWriter*> vertex_data_writers;
	//std::list<CellDataWriter*> cell_data_writers;
	//std::list<VertexDataWriter*> vertex_data_writers;

	
	VertexDataWriter* createVertexDataWriter( const std::string& identifier, int recordsPerVertex ) override;
		//{ return new VertexDataWriter(*this, identifier, recordsPerVertex); }
	CellDataWriter*   createCellDataWriter( const std::string& identifier, int recordsPerCell ) override;
		//{ return new CellDataWriter(*this, identifier, recordsPerCell); }

	VertexWriter*   createVertexWriter() override { return new VertexWriter(); }
	CellWriter*     createCellWriter() override { return new CellWriter(); }
};

#endif /* __EXAHYPE_PLOTTERS_ASCII_CSV_UNSTRUCTURED_GRID_WRITER_SVEN__ */
