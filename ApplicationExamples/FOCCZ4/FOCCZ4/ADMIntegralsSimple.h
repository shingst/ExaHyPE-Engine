#ifndef __EXAHYPE_CCZ4_POSTPROCESSING_SIMPLE_ADM_INTEGRALS__
#define __EXAHYPE_CCZ4_POSTPROCESSING_SIMPLE_ADM_INTEGRALS__

#include "exahype/plotters/ascii/CSVWriter.h"

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"
#include <vector>

/**
 * The "simple" ADM Integrals classes implement the surface integrals on a sphere, implementing
 * rather simple equations (without conformal decomposition). This is described on
 * https://arxiv.org/pdf/gr-qc/0701123.pdf (or in any textbook), and I call it simple compared to
 * what they suggest in terms of volume integrals and quasi-excision around punctures.
 *
 * Following good ExaHyPE practice, the implementation here is merely independent of ExaHyPE and
 * serves as a library which can be used in an ExaHyPE plotter (here: ADMIntegralsWriter).
 **/
namespace ADMIntegralsSimple {
	struct Sphere {
		static constexpr int sphereDim = 2;
		
		const double radius;             ///< Radius (in code units, i.e. M_Sol) where to do the surface integral
		const int num_points[sphereDim]; ///< Number of equidistant grid in the 2-manifold where to integrate on
		int measured_points[sphereDim];  ///< How many points where measured by a grid sweep in each dimension
		
		static constexpr int numIntegrands = 4; /// ADM Mass scalar and Angular Momentum vector
		static constexpr int idxM = 0; /// ADM Mass
		static constexpr int idxJ = 1; /// ADM Angular Momentum with three quantities
		double I[numIntegrands]; ///< The Integrand values determined in a grid sweep

		Sphere(double radius, int ntheta=50, int nphi=50); /// < Initialize the Spherical Integration
		void reset(); ///< Set Quantities and measured points to zero

		
		void plotPatch(const double* const offsetOfPatch, const double* const sizeOfPatch, double* uPatch, double timeStamp);
		void Integrand(const double* const pos, const double* const dx, const double* const Q, double* I);
		
		double coverage() const; ///< Give a number between [0-1] which indicates the fraction of points which have been integrated
		
		std::vector<double> times; ///< Internal bookkeeping for times in previous patches
		void assertTimes(double timeStamp); ///< Fails with exception if there is local time stepping
		
		std::string toString() const; /// Short identifier of constant data hold by this Sphere() instance
	}; // struct Sphere
	
	class MultipleSpheresWriter {
		exahype::plotters::ascii::CSVWriter writer;
		
		struct Line {
			int plotindex;
			double time;
			double radius;
			double M;
			double J0, J1, J2;
			double coverage;
		};
		
		Line cur_line;
	public:
		std::vector<Sphere> spheres;

		MultipleSpheresWriter(
			const std::string filenamebase,
			const std::string reductionsSuffix=".asc");
		
		// Syntactic sugar
		void addSphere(double radius, int ntheta=50, int nphi=50) {
            spheres.push_back(Sphere(radius,ntheta,nphi));
        }

		void startPlotting(double current_time);
		void plotPatch(const double* const offsetOfPatch, const double* const sizeOfPatch, double* uPatch, double timeStamp);
		void finishPlotting();
	}; // class MultipleSpheresWriter
} // ns ADMIntegralsSimple
#endif /* __EXAHYPE_CCZ4_POSTPROCESSING_SIMPLE_ADM_INTEGRALS__ */
