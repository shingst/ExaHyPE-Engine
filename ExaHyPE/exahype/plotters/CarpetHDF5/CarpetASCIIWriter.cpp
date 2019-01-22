#include "exahype/plotters/CarpetHDF5/CarpetASCIIWriter.h"
#include "peano/utils/Loop.h" // dfor
#include "tarch/parallel/Node.h" // for basic MPI rank determination only
#include <sstream>
#include <stdexcept>

typedef tarch::la::Vector<DIMENSIONS, double> dvec;
typedef tarch::la::Vector<DIMENSIONS, bool> boolvec;
typedef tarch::la::Vector<DIMENSIONS, int> ivec;

struct exahype::plotters::CarpetASCIIDatasets::CompatFull { int it, tl, rl, c, ml; };
struct exahype::plotters::CarpetASCIIDatasets::Coord1D { int ix; double time, x; };
struct exahype::plotters::CarpetASCIIDatasets::Coord2D { int ix, iy; double time, x, y; };
struct exahype::plotters::CarpetASCIIDatasets::Coord3D { int ix, iy, iz; double time, x, y, z; };

using namespace exahype::plotters::CarpetASCIIDatasets;

static std::vector<exahype::plotters::ascii::CSVWriter::Column> carpet_writer_compat_full = {
	CSVWRITER_INTEGER_COLUMN(CompatFull, it, "iteration (in ExaHyPE: plotting step, not time integration step)"),
	CSVWRITER_INTEGER_COLUMN(CompatFull, tl, "time level (in ExaHyPE: no use)"),
	CSVWRITER_INTEGER_COLUMN(CompatFull, c,  "component (in ExaHYPE: cell number)"),
	CSVWRITER_INTEGER_COLUMN(CompatFull, ml, "multigrid level (in ExaHYPE: Limiter status)")
};
static std::vector<exahype::plotters::ascii::CSVWriter::Column> carpet_writer_coord1d = {
	CSVWRITER_INTEGER_COLUMN(D, ix, "X-Index of data point (in ExaHYPE: index within a patch)"),
	CSVWRITER_DOUBLE_COLUMN(D, time, "Time of vertex"),
	CSVWRITER_DOUBLE_COLUMN(D, x, "Coordinate of vertex")
};
static std::vector<exahype::plotters::ascii::CSVWriter::Column> carpet_writer_coord2d = {
	CSVWRITER_INTEGER_COLUMN(D, ix, "X-Index of data point (in ExaHYPE: index within a patch)"),
	CSVWRITER_INTEGER_COLUMN(D, iy, "Y-Index of data point (in ExaHYPE: index within a patch)"),
	CSVWRITER_DOUBLE_COLUMN(D, time, "Time of vertex"),
	CSVWRITER_DOUBLE_COLUMN(D, x, "Coordinate of vertex"),
	CSVWRITER_DOUBLE_COLUMN(D, y, "Coordinate of vertex")
};
static std::vector<exahype::plotters::ascii::CSVWriter::Column> carpet_writer_coord3d = {
	CSVWRITER_INTEGER_COLUMN(D, ix, "X-Index of data point (in ExaHYPE: index within a patch)"),
	CSVWRITER_INTEGER_COLUMN(D, iy, "Y-Index of data point (in ExaHYPE: index within a patch)"),
	CSVWRITER_INTEGER_COLUMN(D, iz, "Z-Index of data point (in ExaHYPE: index within a patch)"),
	CSVWRITER_DOUBLE_COLUMN(D, time, "Time of vertex"),
	CSVWRITER_DOUBLE_COLUMN(D, x, "Coordinate of vertex"),
	CSVWRITER_DOUBLE_COLUMN(D, y, "Coordinate of vertex"),
	CSVWRITER_DOUBLE_COLUMN(D, z, "Coordinate of vertex")
};


/**
 * Opens or switchs the currently active H5 file or the list of H5 files.
 * 
 * ATTENTION: The composal of the Carpet filename has to be chosen carefully. Some (!) readers except
 *    for 1D files a pattern *.{x,y,z}.asc
 *    for 2D files a pattern *.{xy,yz,xz}.asc
 *    for 3D files a pattern *.xyz.asc
 *    (all patterns here are bash linux command line globbing patterns)
 * This is especially for the readers from David Radice (https://bitbucket.org/dradice/scidata)
 * and Wolfgang Kastaun (https://bitbucket.org/DrWhat/pycactuset).
 **/
void exahype::plotters::CarpetASCIIWriter::openFile() {
	// TODO: Call CarpetWriter::openFile()		_writer->columns=
	for(auto& file : files) {
		// @TODO: MPI Rank number should go in here.
		local_filename = prefix + (allUnknownsInOneFile ? "" : (sep + writtenQuantitiesNames[writtenUnknown])) + suffix;

		logInfo("open", "Opening File '"<< local_filename << "'");
		file = new exahype::plotters::ascii::CSVWriter(local_filename);
		switch(dim) {
			case 3:
			file->columns = carpet_writer_coord3d;
			break;

			case 2:
			file->columns = carpet_writer_coord2d;
			break;
			
			case 1:
			file->columns = carpet_writer_coord1d;
			break;
		}
		file->writeCommentLine("Created by ExaHyPE/CarpetASCIIWriter");
		file->writeHeader();
		writtenUnknown++;
	}
}

void exahype::plotters::CarpetASCIIWriter::closeFile() {
	// flush in any case, even when closing the file. Cf. http://stackoverflow.com/a/31301117
	flushH5();
	for(auto& file : files) {
		if(file) {
			delete file;
			file = nullptr;
		}
	}
}

void exahype::plotters::CarpetASCIIWriter::flushFile() {
	for(auto& file : files) {
		if(file) {
			file->os.flush();
		}
	}
}

/**
 * This is 2D and 3D, allows several unknowns, named fields and all that.
 **/
void exahype::plotters::CarpetASCIIWriter::plotPatch(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp, int limiterStatus) {

	
	
	for(int writtenUnknown=0; writtenUnknown < writtenUnknowns; writtenUnknown++) {
		exahype::plotters::ascii::CSVWriter* target = files[allUnknownsInOneFile ? 0 : writtenUnknown];
		plotPatchForSingleUnknown(offsetOfPatch, sizeOfPatch, dx, mappedCell, timeStamp, limiterStatus, writtenUnknown, target);
	} // for writtenUnknown
	component++;
}

/**
 * ATTENTION at  composing the HDF5 data field names for the Carpet file format!
 * Many (actually all I know) codes parse these field names with regexps which are sometimes not too flexible.
 * For instance, Visit requries the presence of "::" in the field name.
 * Scidata parses with this python regexp: r"(\w+:?:?\w*\[?\d*\]?) it=(\d+) tl=(\d+) rl=(\d+) c=(\d+)$"
 * That means the whitespace is crucial and has to be exactly like this.
 **/
void exahype::plotters::CarpetASCIIWriter::plotPatchForSingleUnknown(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp, int limiterStatus_data,
      int writtenUnknown, H5::H5File* target) {
	assertion(target != nullptr);
	
	const std::string name = qualifiedWrittenQuantitiesNames[writtenUnknown];
	// Attention: I removed "tl=0 m=0 rl=0" => "tl=0 rl=0" for scidata.
	std::stringstream component_name;
	component_name << name << " it=" << iteration << " tl=0 rl=0 c=" << component;
	
	// 1.) Adopt the continous storage computed before.
	double *componentPatch = new double[singleFieldIdx->size];
	// assume it is finished
	
	// 2. write the data
	switch(dim) {
		case 1:
			for(int i=0; i<basisSize; i++) {
				coord_data = new Coord1D;
				Coord1D.ix = i;
				Coord1D.x = offsetOfPatch + i*sizeOfPatch)
			}
			break;
	}
	
	dfor(i,basisSize) {
		dvec pos = offsetOfPatch + i.convertScalar<double>()* (sizeOfPatch(0)/(basisSize+1));

		
		switch(dim) {
			case 3:

			
			componentPatch[singleFieldIdx->get(i(2),i(1),i(0))] = mappedCell[writtenCellIdx->get(i(2),i(1),i(0),writtenUnknown)];
			break;

			case 2:
			componentPatch[singleFieldIdx->get(i(1),i(0))] = mappedCell[writtenCellIdx->get(i(1),i(0),writtenUnknown)];
			break;
			
			case 1:

			componentPatch[singleFieldIdx->get(i(0))] = mappedCell[writtenCellIdx->get(i(0),writtenUnknown)];
			break;
		}
	}
	
	DataSet table = target->createDataSet(component_name.str(), PredType::NATIVE_FLOAT, patch_space);
	table.write(componentPatch, PredType::NATIVE_DOUBLE);
	
	delete[] componentPatch;
	
	// write all meta information about the table
	
	offsetOfPatch

	Attribute origin = table.createAttribute("origin", PredType::NATIVE_DOUBLE, dtuple);
	origin.write(PredType::NATIVE_DOUBLE, offsetOfPatch.data());
	
	
	int level_data = 0; // TODO: Read out real (AMR) cell level (needs to be passed from ADERDG2CarpetHDF5 class, for instance)
	Attribute level = table.createAttribute("level", PredType::NATIVE_INT, H5S_SCALAR);
	level.write(PredType::NATIVE_INT, &level_data);

	// This is not a Carpet metadata but something we add from the ExaHyPE side.
	// It tells about the limiter status and typically reads like 0-O,1..2-DG,3..4-FV,5-T
	Attribute limiterStatus = table.createAttribute("limiterStatus", PredType::NATIVE_INT, H5S_SCALAR);
	limiterStatus.write(PredType::NATIVE_INT, &limiterStatus_data);
	
	Attribute timestep = table.createAttribute("timestep", PredType::NATIVE_INT, H5S_SCALAR);
	timestep.write(PredType::NATIVE_INT, &iteration);

	Attribute time = table.createAttribute("time", PredType::NATIVE_DOUBLE, H5S_SCALAR);
	time.write(PredType::NATIVE_DOUBLE, &timeStamp);

	// dx in terms of Cactus: Real seperation from each value
	Attribute delta = table.createAttribute("delta", PredType::NATIVE_FLOAT, dtuple);
	delta.write(PredType::NATIVE_DOUBLE, dx.data()); // issue: conversion from double to float
	
	StrType t_str = H5::StrType(H5::PredType::C_S1, name.size()+1); // Todo: not sure about +1 for \0
	Attribute aname = table.createAttribute("name", t_str, H5S_SCALAR);
	aname.write(t_str, name.c_str());
	
	int cctk_nghostzones_data = 0;  // no ghostzones for ADERDG; we won't save ghostzones for FV solvers.
	Attribute nghostzones = table.createAttribute("cctk_nghostzones", PredType::NATIVE_INT, H5S_SCALAR);
	nghostzones.write(PredType::NATIVE_INT, &cctk_nghostzones_data);
}


#endif /* HDF5 */
