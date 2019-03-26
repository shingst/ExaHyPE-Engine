#include "GridDemonstrator_ADERDG.h"

#include "GridDemonstrator_ADERDG_Variables.h"

#include <cmath>
#include <random>


tarch::logging::Log GridDemonstrator::GridDemonstrator_ADERDG::_log( "GridDemonstrator::GridDemonstrator_ADERDG" );

struct RandomInitialDataPerTimestep {
	constexpr static int nx = 20; // grid points we want to store in advance in each direction
	constexpr static int nt = 100; // Time levels we want to store in advance. Is treated periodically.
	constexpr static double dx = 0.1; // grid distance in each direction
	constexpr static double dt = 0.5; // time distance
	constexpr static int numpoints = nt*nx*nx*nx;  // for simplicity, always 3d
	double out[numpoints];
	static double rand01() { return ((double) rand() / (RAND_MAX)); } ///< random \in [0,1]
	RandomInitialDataPerTimestep() {
		for(int i=0; i<numpoints; i++) out[i] = rand01();		
	}
	double interpolate(const double* const x, double t) {
		int ti, xi[3], index;
		ti = int(std::floor(t / dt)) % nt;
		for(int i=0; i<DIMENSIONS; i++) xi[i] = int(std::floor(x[i] / dx)) % nx;
		if(DIMENSIONS<3) xi[2] = 0;
		index = nt*ti;
		for(int i=0; i<DIMENSIONS; i++) index += xi[i]*nx;
		return out[index];
	}
};

RandomInitialDataPerTimestep* randId;

void GridDemonstrator::GridDemonstrator_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
	// Create random data for all time steps well in advance.
	// This avoids any problems with MPI and TBB.
	randId = new RandomInitialDataPerTimestep;s
}

bool GridDemonstrator::GridDemonstrator_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {

	// admissible means the solution is acceptable as being a physical solution.
	bool isAdmissible;
	
	typedef tarch::la::Vector <DIMENSIONS,double> vec;
	
	vec abscenter = tarch::la::abs(center);
	
	// lower left, upper right radius of cell
	double l = tarch::la::norm2(abscenter - dx/2.0);
	double r = tarch::la::norm2(abscenter + dx/2.0);
	// center radius of cell
	double c = tarch::la::norm2(center);

	// works well with maximum-mesh-depth = 3
	// and grid size = [-20,20]
	bool singleNS = true;
	if(singleNS) {
		constexpr double radius_shell = 13; // limit around this shell
		isAdmissible = (l > radius_shell) || (r <= radius_shell);
		//printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");
	}
	
	// a single NS with some atmosphere thickness
	bool atmoNS = true;
	if(atmoNS) {
		double rnd = solution[0]; // randomness moved into adjustPointSolution, isPhysicallyAdmissible shall be deterministic
		
		// limit around this shell:
		double radius_shell = 10 + rnd * 4;
		double atmo_thickness = 2;
		isAdmissible = (l > radius_shell+atmo_thickness/2) || (r <= radius_shell-atmo_thickness/2);
		//printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");
	}
	
	// works well with maximum-mesh-depth = 5
	// and grid size = 20
	bool singleBH = false;
	if(singleBH) {
		constexpr double horizon_size = .5;
		isAdmissible = c > horizon_size;
	}
	
	// two points rotating around the center, roughly representing
	// the inspiral phase of two NS.
	// Status: Does not work in the moment. There is a bug somewhere.
	bool inspiralNS = false;
	constexpr double pi = M_PI;
	constexpr double radius = 6; // size of each NS
	if(inspiralNS) {
		vec pol[2];  // centers in polar coordinates {r,phi}

		double rstart = 15, rziel = 2, tstart = 0, tziel = 10*dt;
		double r = rstart + (rstart - rziel) / (tstart - tziel) * t;
		
		double runden = 5, T = tziel / runden, omega = 2*pi / T, offset = pi;
		double phi = omega * t;

		pol[0] = r, phi;
		pol[1] = r, phi + offset;
		
		bool Admissible[2];
		for(int i=0; i<2; i++) {
			double r=pol[i][0], phi = pol[i][1];
			
			// determine center of star in cartesian coordinates
			vec star;
			star[0] = r * std::cos(phi); // x
			star[1] = r * std::sin(phi); // y

			// determine radius positions of current cell center relative to star center
			
			vec relcenter = tarch::la::abs(center - star);
			// lower left, upper right radius in these coordinates
			l = tarch::la::norm2(relcenter - dx/2.0);
			r = tarch::la::norm2(relcenter + dx/2.0);
			
			Admissible[i] = (l > radius) || (r <= radius);
		}
		
		isAdmissible = (Admissible[0] && Admissible[1]);		
	}

	return isAdmissible;

}



void GridDemonstrator::GridDemonstrator_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 1 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    Q[0] = 0.0;
  }
  
  // to introduce per timestep randomness:
  Q[0] = randId->interpolate(x,t);
}

void GridDemonstrator::GridDemonstrator_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 1 + 0

  // @todo Please implement/augment if required
  stateOut[0] = 0.0;

  fluxOut[0] = 0.0;
}

exahype::solvers::Solver::RefinementControl GridDemonstrator::GridDemonstrator_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void GridDemonstrator::GridDemonstrator_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 1 + 0
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
}




//You can either implement this method or modify fusedSource
void GridDemonstrator::GridDemonstrator_ADERDG::algebraicSource(const double* const Q,double* S) {
  // @todo Please implement/augment if required
  S[0] = 0.0;
}



