/*
 * Inline Time Series reductions as helpers for UserOnTheFlyPostProcessing classes.
 * MPI version, 2017-02-20, SvenK.
 */
#ifndef __TIMESERIES_EXA_MPI__
#define __TIMESERIES_EXA_MPI__

#include <vector>
#include <map>
#include <limits>
#include <cmath>

#include <stdio.h>
#include <string.h>

#ifdef Parallel
#include "tarch/parallel/Node.h" // for MPI
#include "tarch/parallel/NodePool.h" 
#endif

#include "tarch/logging/Log.h"
#include "kernels/GaussLegendreBasis.h" // for gaussLegendreWeights

namespace reductions {

/**
 * A helper function to determine the volume element (dV=dx*dy*dz) for an integration
 * in the ADER scheme where the Gauss Legendre interpolation points are not even equally
 * distributed.
 * 
 * This function could easily be  constexpr if gaussLegendreWeights would be a constexpr.
 **/
inline double ADERDGVolume(const int order, const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, const tarch::la::Vector<DIMENSIONS, int>& pos) {
	// Gauss-Legendre weights from pos argument
	double wx = kernels::legendre::weights[order][pos[0]];
	double wy = kernels::legendre::weights[order][pos[1]];
	double wz = 1;
	#ifdef Dim3
	wz = kernels::legendre::weights[order][pos[2]];
	#endif

	// volume form for integration
	double scaling = tarch::la::volume(sizeOfPatch);
	
	return scaling*wx*wy*wz;
}

/**
 * Non-MPI aware single-core reductor.
 * 
 * TimeSeriesReductions is a C++ class to be used in UserOnTheFlyPostProcessing classes
 * to compute reductions of fields and write them out as ASCII files. The reductions are
 * 
 *   l1: L^1 norm of field (l1 = sum of all cells)
 *   l2: L^2 norm of field (l2 = sqrt(sum of cell[i]^2))
 *   max: Maximum value of field (equal to L^\inf norm)
 *   min: Minimum value of field
 *   avg: Average value of field
 *
 * 
 * @author SvenK
 * 
 **/
class TimeSeriesReductions {
public:
    enum index                          { tidx=0    , time=1, l1=2   , l2=3   , max=4, min=5, avg=6, sum=7 , nelem=8     , LEN=9 };
    const char * const colnames[LEN]  = {"plotindex","time" ,"l1norm","l2norm","max" ,"min" ,"avg" , "sum" , "numelements", };
    const char * const colformat[LEN] = {"%.0f\t"   ,"%e\t" ,"%e\t"  ,"%e\t"  ,"%e\t","%e\t","%e\t", "%e\t", "%.0f\t"     , };
    double data[LEN];
    //int avgcnt; // now moved to nelem for MPI_Send primitive. If you don't want it, add another
    // const bool const colmask[LEN] = { PRINT, SKIP, PRINT, PRINT, ...} with PRINT=true, SKIP=false or so.
    // (better use colprint[LEN]={true,false,...}

    TimeSeriesReductions() {
	data[tidx] = -1.0; // to start with 0 at first startRow
	// to debug whether strange things happen
	data[l1] = 42.123;
    }

    void startRow(double current_time) {
	data[tidx] += 1.0;
	data[time] = current_time;
	data[l1] = 0.0;
	data[l2] = 0.0;
	data[max] = 0.0;
	data[min] = std::numeric_limits<double>::max();
	data[avg] = 0.0;
	data[sum] = 0.0;
	data[nelem] = 0.0;
    }

    void addValue(double val, double dx) {
	// dx is the scaling (volume form) as computed by peano::la::volume,
        // compare the calculation in MyEulerSolver_Plotter1.cpp 
        data[l1] += std::abs(val) * dx;
        data[l2] += val * val * dx;
        data[max] = std::max( data[max], val);
        data[min] = std::min( data[min], val);
        data[avg] += val;
	data[sum] += val * dx;
        data[nelem] += 1.0; // double, unfortunately
    }

    void addValues(double input[LEN]) {
        data[l1] += input[l1];
        data[l2] += input[l2];
        data[max] = std::max( data[max], input[max] );
        data[min] = std::min( data[min], input[min] );
        data[avg] += input[avg];
	data[sum] += input[sum];
	data[nelem] += input[nelem];
    }

    void addValues(TimeSeriesReductions& input) {
	addValues(input.data);
    }

    void finishRow() {
	data[l2] = std::sqrt(data[l2]);
	data[avg] = data[avg] / data[nelem];
    }
};


#define master tarch::parallel::Node::getInstance().isGlobalMaster()

/**
 * Actual MPI-aware writer of reductions over all ranks.
 * 
 * @author SvenK
 * 
 **/
class ReductionsWriter : public TimeSeriesReductions {
    FILE* asc;
    int reductionTag;
public:
    const std::string filename;

    ReductionsWriter(const std::string _filename) : filename(_filename) {
	if(master)
		openFile();

	#ifdef Parallel
	reductionTag = tarch::parallel::Node::getInstance().reserveFreeTag(
		std::string("TimeSeriesReductions(")+filename+")::finishRow()");
	#endif

    }
    
    ~ReductionsWriter() {
	if(master && asc)
		fclose(asc);

	#ifdef Parallel
	tarch::parallel::Node::getInstance().releaseTag(reductionTag);
	#endif
    }

    void openFile() {
	static tarch::logging::Log _log("ReductionsWriter");
	logInfo("openFile()", "Trying to open file '" << filename << "' in the ReductionsWriter"); // this is strange. When not touching the string, it gets removed...
	const char* fn = filename.c_str();
	asc = fopen(fn, "w");
	if(asc == NULL) {
		//fprintf(stderr, "ASCII writer: Could not open output file '%s'\n", fn);//, strerror(errno));
		// in MPI environment:
		logError("openFile()", "Cannot open output ASCII file at '" << filename << "'.");
		exit(-1);
	}
	
	// print the file header line
	// fprintf(asc, "# ");
	for(int i=0; i<LEN; i++) {
		fputs(colnames[i], asc);
		fputs(" ", asc);
	}
	fprintf(asc, "\n");
    }

    void writeRow() {
	if(asc == NULL) {
		static tarch::logging::Log _log("ReductionsWriter");
		logError("writeRow()", "Oh no.");
		exit(-1);
	}
	writeRow(asc);
    }

    void writeRow(FILE* stream) {
	for(int i=0; i<LEN; i++) {
            fprintf(stream, colformat[i], data[i]);
        }
        fprintf(stream, "\n");
	fflush(stream); // write out this line immediately
    }

    void finishRow() {
	#ifdef Parallel
		if(master) {
			double recieved[LEN];
			for (int rank=1; rank<tarch::parallel::Node::getInstance().getNumberOfNodes(); rank++) {
				if(!tarch::parallel::NodePool::getInstance().isIdleNode(rank)) {
					MPI_Recv( &recieved[0], LEN, MPI_DOUBLE, rank, reductionTag, tarch::parallel::Node::getInstance().getCommunicator(), MPI_STATUS_IGNORE );
					addValues(recieved);
				}
			}
		} else {
			MPI_Send( &data[0], LEN, MPI_DOUBLE, tarch::parallel::Node::getGlobalMasterRank(), reductionTag, tarch::parallel::Node::getInstance().getCommunicator());
		}
	#endif
	if(master) {
		TimeSeriesReductions::finishRow();
		writeRow();
	}
    }
};


/**
 * Convenience class to organize reducing multiple fields in a single call.
 * This class basically manages a list of ReductionWriters, each with names and associated
 * field indices in the vector of conserved variables Q. This allows quick mapping of
 * a vector onto reduction files.
 * 
 * 
 * @author SvenK
 * 
 **/
class MultipleReductionsWriter {
	std::vector<ReductionsWriter*> reductions;
	std::vector<int> positionsQ; ///< the positions in reductions corresponding to the reductions

	/// sth like "output/" or "output/reductions-" to put before each filename
	std::string reductionsPrefix;
	
	/// sth. like ".asc" to put afer each filename
	std::string reductionsSuffix;

public:	
	MultipleReductionsWriter(const std::string& _reductionsPrefix="", const std::string& _reductionsSuffix=".asc") :
		reductionsPrefix(_reductionsPrefix),
		reductionsSuffix(_reductionsSuffix)
	{}
	
	std::string composeFilename(const std::string& reductionName) {
		return reductionsPrefix + reductionName + reductionsSuffix;
	}
	
	int size() {
		return reductions.size();
	}

	/**
	 * Make sure your vector Q is always larger than any posq you pass here.
	 * 
	 **/
	void add(int posq, const std::string& redname) {
		reductions.push_back(new ReductionsWriter(composeFilename(redname)));
		positionsQ.push_back(posq);
	}
	
	void startRow(double current_time) {
		for (auto& r : reductions)
			r->startRow(current_time);
	}
	
	void finishRow() {
		for (auto& r : reductions)
			r->finishRow();
	}
	
	void addValue(double val, double dx) {
		for (auto& r : reductions)
			r->addValue(val, dx);
	}
	
	void addValue(double Q[], double dx) {
		for (size_t i = 0; i < reductions.size(); i++)
			reductions[i]->addValue(Q[ positionsQ[i] ], dx);
	}
};


} // namespace reductions

#endif /* __TIMESERIES_EXA_MPI__ */
