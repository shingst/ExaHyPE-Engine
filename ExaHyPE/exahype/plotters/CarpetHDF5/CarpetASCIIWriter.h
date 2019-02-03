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
 *
 * @authors: Sven Koeppel
 **/

#ifndef _EXAHYPE_PLOTTERS_CARPET_ASCII_WRITER_
#define _EXAHYPE_PLOTTERS_CARPET_ASCII_WRITER_

#include "exahype/plotters/Plotter.h"
#include "exahype/plotters/CarpetHDF5/CarpetWriter.h"
#include "exahype/plotters/slicing/CartesianSlicer.h"
#include "kernels/KernelUtils.h" // idx::kernels

namespace exahype {
  namespace plotters {
    class CarpetASCIIWriter;

    namespace CarpetASCIIDatasets {
        struct CompatFull;
	struct Coord1D;
	struct Coord2D;
	struct Coord3D;
    }

    namespace ascii { class CSVStackWriter; } // external forward decl, #include exahype/plotters/ascii/CSVWriter.h
  }
}

/**
 * <h2>Writing CarpetHDF5 files which are compatible to Cactus/EinsteinToolkit</h2>
 * 
 * This writer allows you to write Carpet ASCII files, you don't need HDF5 to use it.
 * 
 * @author Sven Köppel
 *
 **/
class exahype::plotters::CarpetASCIIWriter : public exahype::plotters::CarpetWriter {
  tarch::logging::Log _log;

public:
  typedef tarch::la::Vector<DIMENSIONS, double> dvec;
  typedef tarch::la::Vector<DIMENSIONS, bool> boolvec;
  typedef tarch::la::Vector<DIMENSIONS, int> ivec;

	
  bool fullCarpetCompatibility; ///< Generate all columns of carpet


  /**
   * cf. also the documentation in the ADERDG2CarpetHDF5.h
   * 
   * oneFilePerTimestep: You might want to have this to view the result during computation
   *     as HDF5 is very lazy at writing. Note that writing out data this form is not compilant
   *     with most CarpetHDF5 readers (ie. the visit reader). You must join seperate HDF5 files afterwards
   *     manually.
   * 
   * allUnknownsInOneFile: Write different fields in a single H5 combined file. Typically for Cactus
   *     as structure of arrays is to write each unknown in its own file (ie. one file per physical field).
   *
   **/
  CarpetASCIIWriter(const std::string& _filename, int _basisSize, int _solverUnknowns, int _writtenUnknowns, exahype::parser::ParserView _plotterParameters,
		   char** writtenQuantitiesNames, bool oneFilePerTimestep_=false, bool allUnknownsInOneFile_=false);

  virtual void openFile(); ///< Opens or switchs the currently active file or the list of files. Closes if neccessary.
  virtual void flushFile(); ///< Flushs all file output buffers. Always flushs before.
  virtual void closeFile(); ///< Closes all files. Closes, deletes and nulls the file objects.

  typedef exahype::plotters::ascii::CSVStackWriter  CSVWriterType;
  std::vector<CSVWriterType> files; ///< List of pointers to H5Files. Has length 1 if allUnknownsInOneFile.

  void startPlotting(double time);
  void finishPlotting();

  // Default values for limiterStatus for plotPatch* functions, used for instance from
  // a pure FV solver which has no limiter status flag.
  constexpr static int nonLimitingLimiterStatus = -1;

  void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      double* mappedCell, double timeStamp, int limiterStatus=nonLimitingLimiterStatus);
  
  void plotPatchForSingleUnknown(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp, int limiterStatus_data,
      int writtenUnknown, CSVWriterType& target);
}; // class

#endif /* _EXAHYPE_PLOTTERS_CARPET_ASCII_WRITER_ */
