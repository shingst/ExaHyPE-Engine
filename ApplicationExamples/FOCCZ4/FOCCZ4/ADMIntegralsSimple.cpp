#include "ADMIntegralsSimple.h"

// kernels::interpolate for DG evaluation on an arbitrary point in a patch
#include "kernels/GaussLegendreBasis.h"

// Importing compile-time constants from the local ExaHyPE application:
#include "AbstractFOCCZ4Solver_ADERDG.h"
#include "FOCCZ4Solver_ADERDG_Variables.h"
namespace idx = FOCCZ4::FOCCZ4Solver_ADERDG_Variables::shortcuts;
constexpr int nVar  = FOCCZ4::AbstractFOCCZ4Solver_ADERDG::NumberOfVariables;
constexpr int order = FOCCZ4::AbstractFOCCZ4Solver_ADERDG::Order;

#define TDIM DIMENSIONS
#include "SVEC/tensish.cpph"     // We use tensish for simplicity
//#include "InitialData/ADMBase.h" // We use ADMBase for nonconformal computation

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <sstream>

namespace ADMIntegralsSimple { // just for keeping everyhing constrained

ADMIntegralsSimple::Sphere::Sphere(double radius, int ntheta, int nphi) :
	radius(radius), num_points{ntheta,nphi}
{
	reset();
}

void ADMIntegralsSimple::Sphere::reset() {
	TDO(i, sphereDim) measured_points[i] = 0;
	TDO(i, numIntegrands) I[i] = 0;
	times.clear();
}

constexpr double pi = M_PI;

using namespace ADMIntegralsSimple;
using namespace std;
using namespace tensish;
//using namespace FOCCZ4_InitialData; // ADMBase

using array3 = std::array<double,3>; // instead of using dvec = tarch::la::Vector<DIMENSIONS, double>;

// helpers:

inline array3 cartesian(const double R, const double theta, const double phi) {
	array3 pos;
	pos[0] = R * sin(theta) * cos(phi);
	pos[1] = R * sin(theta) * sin(phi);
	pos[2] = R * cos(theta);
	return pos;
}

namespace SphericalCoordinates { ///< to avoid magic numbers
	constexpr int r     = 0;
	constexpr int theta = 1;
	constexpr int phi   = 2;
}


#define eps  tensish::sym::levi  // epsilon function alias

// These are stolen from ADMBase.cpp
using metric_t = tensish::metric::generic<square::const_shadow_D>;
using auxGamma_t = tensish::generic::stored<square::const_shadow<sym::row_major, DIMENSIONS>, DIMENSIONS>;

void ADMIntegralsSimple::Sphere::Integrand(const double* const xpos, const double* const dxpos, const double* const Q, double* I) {
	using namespace SphericalCoordinates;

	enum { GTXX= 0, GTXY, GTXZ, GTYY, GTYZ, GTZZ};
	enum { ATXX= idx::K, ATXY, ATXZ, ATYY, ATYZ, ATZZ};
	
	double ginv[3][3];
	ginv[0][0] = -Q[idx::G+GTYZ] * Q[idx::G+GTYZ] + Q[idx::G+GTYY] * Q[idx::G+GTZZ];
	ginv[0][1] =  Q[idx::G+GTYZ] * Q[idx::G+GTXZ] - Q[idx::G+GTXY] * Q[idx::G+GTZZ];
	ginv[0][2] = -Q[idx::G+GTXZ] * Q[idx::G+GTYY] + Q[idx::G+GTXY] * Q[idx::G+GTYZ];
	ginv[1][1] = -Q[idx::G+GTXZ] * Q[idx::G+GTXZ] + Q[idx::G+GTXX] * Q[idx::G+GTZZ];
	ginv[1][2] =  Q[idx::G+GTXY] * Q[idx::G+GTXZ] - Q[idx::G+GTXX] * Q[idx::G+GTYZ];
	ginv[2][2] = -Q[idx::G+GTXY] * Q[idx::G+GTXY] + Q[idx::G+GTXX] * Q[idx::G+GTYY];

	ginv[1][0] = ginv[0][1];
	ginv[2][0] = ginv[0][2];
	ginv[2][1] = ginv[1][2];

	auto const dxG =&(Q[idx::dxG]);
	auto const dyG =&(Q[idx::dyG]);
	auto const dzG =&(Q[idx::dzG]);

	double Dx[3][3];
	Dx[0][0] = dxG[GTXX];
	Dx[0][1] = dxG[GTXY];
	Dx[0][2] = dxG[GTXZ];
	Dx[1][0] = dyG[GTXX];
	Dx[1][1] = dyG[GTXY];
	Dx[1][2] = dyG[GTXZ];
	Dx[2][0] = dzG[GTXX];
	Dx[2][1] = dzG[GTXY];
	Dx[2][2] = dzG[GTXZ];

	double Dy[3][3];
	Dy[0][0] = dxG[GTXY];
	Dy[0][1] = dxG[GTYY];
	Dy[0][2] = dxG[GTYZ];
	Dy[1][0] = dyG[GTXY];
	Dy[1][1] = dyG[GTYY];
	Dy[1][2] = dyG[GTYZ];
	Dy[2][0] = dzG[GTXY];
	Dy[2][1] = dzG[GTYY];
	Dy[2][2] = dzG[GTYZ];

	double Dz[3][3];
	Dz[0][0] = dxG[GTXZ];
	Dz[0][1] = dxG[GTYZ];
	Dz[0][2] = dxG[GTZZ];
	Dz[1][0] = dyG[GTXZ];
	Dz[1][1] = dyG[GTYZ];
	Dz[1][2] = dyG[GTZZ];
	Dz[2][0] = dzG[GTXZ];
	Dz[2][1] = dzG[GTYZ];
	Dz[2][2] = dzG[GTZZ];

	std::array<double,3> Ghat_d {{0., 0., 0.}};
	
	for(int i=0;i<3; ++i)
  	  for(int j=0;j<3; ++j){
	    Ghat_d[0] += 2.*Dx[i][j] * ginv[i][j];
	    Ghat_d[1] += 2.*Dy[i][j] * ginv[i][j];
	    Ghat_d[2] += 2.*Dz[i][j] * ginv[i][j];
	  }

	std::array<double,3> Ghat {{
	  Ghat_d[0] * ginv[0][0] + Ghat_d[1] * ginv[0][1] + Ghat_d[2] * ginv[0][2],
	  Ghat_d[0] * ginv[0][1] + Ghat_d[1] * ginv[1][1] + Ghat_d[2] * ginv[1][2],
	  Ghat_d[0] * ginv[0][2] + Ghat_d[1] * ginv[2][1] + Ghat_d[2] * ginv[2][2]
	}};

	auto const psi = std::exp(-0.5*Q[idx::phi]);
	auto const psim2 = std::exp(Q[idx::phi]);
	auto const psi2 = 1./psim2;
	auto const psi6 = psi2 * psi2 * psi2;

	std::array<double,3> dpsi {{
	  -0.5*psi*Q[idx::P+0],
	  -0.5*psi*Q[idx::P+1],
	  -0.5*psi*Q[idx::P+2],
	}};

	std::array<double,3> dpsi_u {{
	  dpsi[0] * ginv[0][0] + dpsi[1] * ginv[0][1] + dpsi[2] * ginv[0][2],
	  dpsi[0] * ginv[0][1] + dpsi[1] * ginv[1][1] + dpsi[2] * ginv[1][2],
	  dpsi[0] * ginv[0][2] + dpsi[1] * ginv[2][1] + dpsi[2] * ginv[2][2]
	}};

	auto const M4PI = 1./(4.*M_PI);

	std::array<double,3> integrandM {{
	  2.*M4PI*(1./8. * Ghat[0] - dpsi_u[0]),	  
	  2.*M4PI*(1./8. * Ghat[1] - dpsi_u[1]),	  
	  2.*M4PI*(1./8. * Ghat[2] - dpsi_u[2]),	  
	}};

	auto const r = std::sqrt(xpos[0]*xpos[0] + xpos[1]*xpos[1] + xpos[2]*xpos[2]);
	
	auto costheta = xpos[2]/r;
	if(!std::isnormal(costheta)) costheta = std::copysign(1., xpos[2]);
  	auto const theta = std::acos(costheta);
	auto const sintheta = std::sin(theta);	
	
	std::array<double,3> dS {{
	  xpos[0] * r *sintheta * dxpos[1]*dxpos[2],
	  xpos[1] * r *sintheta * dxpos[1]*dxpos[2],
	  xpos[2] * r *sintheta * dxpos[1]*dxpos[2],
	}};

	I[idxM] = dS[0] * integrandM[0] + dS[1] *integrandM[1] + dS[2] * integrandM[2];
	I[idxJ+0] = 0;
	I[idxJ+1] = 0;


	std::array<double,3> Ax_u {{
	  Q[ATXX] * ginv[0][0] + Q[ATXY] * ginv[0][1] + Q[ATXZ] * ginv[0][2],
	  Q[ATXX] * ginv[0][1] + Q[ATXY] * ginv[1][1] + Q[ATXZ] * ginv[1][2],
	  Q[ATXX] * ginv[0][2] + Q[ATXY] * ginv[2][1] + Q[ATXZ] * ginv[2][2]
	}};

	std::array<double,3> Ay_u {{
	  Q[ATXY] * ginv[0][0] + Q[ATYY] * ginv[0][1] + Q[ATYZ] * ginv[0][2],
	  Q[ATXY] * ginv[0][1] + Q[ATYY] * ginv[1][1] + Q[ATYZ] * ginv[1][2],
	  Q[ATXY] * ginv[0][2] + Q[ATYY] * ginv[2][1] + Q[ATYZ] * ginv[2][2]
	}};

	std::array<double,3> integrandJz {{
	  0.5*M4PI*psi6*(xpos[0] * Ay_u[0] - xpos[1] * Ax_u[0] ),
	  0.5*M4PI*psi6*(xpos[0] * Ay_u[1] - xpos[1] * Ax_u[1] ),
	  0.5*M4PI*psi6*(xpos[0] * Ay_u[2] - xpos[1] * Ax_u[2] ),
	}};


	I[idxJ+2] = dS[0] * integrandJz[0] + dS[1] *integrandJz[1] + dS[2] * integrandJz[2];

//	std::cout << "I[idxJ+2] \t"
//		  << "x\t" << xpos[0] << "\t"
//		  << "y\t" << xpos[1] << "\t"
//		  << "z\t" << xpos[2] << "\t"
//	          << "Value\t"<< I[idxJ+2] << std::endl;

}

void ADMIntegralsSimple::Sphere::plotPatch(const double* const offsetOfPatch,
        const double* const sizeOfPatch, double* uPatch, double timeStamp) {
	using namespace SphericalCoordinates;
	assert(DIMENSIONS == 3);
	assertTimes(timeStamp);
	
	// This code comes from the GRMHD/Writers/SphereIntegrals.cpp, a dumb an inefficient
	// way to sample the 2d sphere with a uniform rectangular grid with ntheta x nphi points
	// and a loop which just checks for each one whether it is contained in the current patch or not.

	const double dx[3] = { /* dr */ 0, /* dtheta */ pi/num_points[0], /* dphi */ 2*pi/num_points[1] };
	int ind[3]; // running indices for rho (not used), theta, phi
	
	array3 pos; // integration point
	for(ind[theta]=0; ind[theta]<num_points[0]; ind[theta]++) {
        for(ind[phi]=0; ind[phi]<num_points[1]; ind[phi]++) {
		// compute the cartesian position of this integration point
		pos = cartesian(radius, ind[theta]*dx[theta], ind[phi]*dx[phi]);
		
		// check if it is inside the current cell
		bool isinside = true;
		DFOR(d) 
            isinside = isinside && offsetOfPatch[d] < pos[d] && pos[d] < (offsetOfPatch[d]+sizeOfPatch[d]);
		if(!isinside) continue;
		
		double Q[nVar], dI[numIntegrands];
		
		// for the time being, evaluate all quantities on the integration point
		TDO(k,nVar) Q[k] = kernels::legendre::interpolate(offsetOfPatch, sizeOfPatch, pos.data(), nVar, k, order, uPatch);

		// Evaluate Integral contribution at point pos
		Integrand(pos.data(), (double*)dx, (double*)Q, (double*)dI);
		
		// add to Integrand value
		for(int nn=idxM; nn < idxJ+3; ++nn){
            		I[nn] += dI[nn];
        	}
        } // for 2d sphere
    }
} // plotPatch

void ADMIntegralsSimple::Sphere::assertTimes(double timeStamp) {
	times.push_back(timeStamp);
	double globalTime = times[0];
	for(auto&& time : times) {
		if(std::abs((time - globalTime)/globalTime) > 1e-15) {
			std::stringstream ss;
			ss << "ADMIntegrals::Sphere does not support local time stepping, but I have patches with times " << time << " and " << globalTime << " which differ.";
			throw std::out_of_range(ss.str());
		}
	}
}


double ADMIntegralsSimple::Sphere::coverage() const {
	double coverage[sphereDim];
	TDO(d,sphereDim) coverage[d] = (double)measured_points[d] / num_points[d];
	
	double totalCoverage = 0;
	TDO(d,sphereDim) totalCoverage += coverage[d] / sphereDim;
	
	return totalCoverage;
}

std::string ADMIntegralsSimple::Sphere::toString() const {
	std::stringstream ss;
	ss << "ADMIntegralsSimple::Sphere(radius="<< radius << ",num_points={"<<num_points[0]<<","<<num_points[1]<<"})";
	return ss.str();
}

//////////////////////// 2ND PART, Administration /////////////////////////////

ADMIntegralsSimple::MultipleSpheresWriter::MultipleSpheresWriter(const std::string filenamebase, const std::string reductionsSuffix) {
	writer.columns = {
		CSVWRITER_INTEGER_COLUMN(Line, plotindex, "The number of plots where this belongs to"),
		CSVWRITER_DOUBLE_COLUMN (Line, time,      "Coordinate time [M]"),
		CSVWRITER_DOUBLE_COLUMN (Line, radius,    "Radius at which the Surface Integral was evaluated"),
		CSVWRITER_DOUBLE_COLUMN (Line, M,         "ADM Mass determined by surface integral (May require a global reduction)"),
		CSVWRITER_DOUBLE_COLUMN (Line, J0,        "ADM Momentum (x) determined by surface integral (May require a global reduction)"),
		CSVWRITER_DOUBLE_COLUMN (Line, J1,        "ADM Momentum (y) determined by surface integral (May require a global reduction)"),
		CSVWRITER_DOUBLE_COLUMN (Line, J2,        "ADM Momentum (z) determined by surface integral (May require a global reduction)"),
		CSVWRITER_DOUBLE_COLUMN (Line, coverage,  "Point coverage by this MPI process. Summing all reductions should give 1."),
	};
	
	writer.openFile(filenamebase, reductionsSuffix);
	writer.writeCommentLine("FOCCZ4::ADMIntegralsSimple::MultipleSpheresWriter output");
	cur_line.plotindex = 0;
}

void ADMIntegralsSimple::MultipleSpheresWriter::startPlotting(double current_time) {
	if(cur_line.plotindex == 0) {
		// Make a list of registered spheres and finish the header
		writer.writeCommentLine("Requested Surfaces for Integration:");
		for(auto&& sphere : spheres) writer.writeCommentLine( sphere.toString() );
		
		writer.writeHeader();
	}
	
	for(auto&& sphere : spheres) sphere.reset();
	cur_line.plotindex += 1.0;
	cur_line.time = current_time;
}

void ADMIntegralsSimple::MultipleSpheresWriter::plotPatch(const double* const offsetOfPatch, const double* const sizeOfPatch, double* uPatch, double timeStamp) {
	for(auto&& sphere : spheres) sphere.plotPatch(offsetOfPatch,sizeOfPatch,uPatch,timeStamp);
}

void ADMIntegralsSimple::MultipleSpheresWriter::finishPlotting() {
	for(auto&& sphere : spheres) {
		cur_line.radius = sphere.radius;
		cur_line.M = sphere.I[sphere.idxM];
		cur_line.J0 = sphere.I[sphere.idxJ+0];
		cur_line.J1 = sphere.I[sphere.idxJ+1];
		cur_line.J2 = sphere.I[sphere.idxJ+2];
		cur_line.coverage = sphere.coverage();
		
		CSVWRITER_WRITE_ROW(writer, cur_line);
	}
}


} // ns ADMIntegralsSimple
