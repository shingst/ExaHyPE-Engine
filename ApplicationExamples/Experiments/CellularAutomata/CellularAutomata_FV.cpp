#include "CellularAutomata_FV.h"
#include "CellularAutomata_FV_Variables.h"

/**
 * Game Of Life / Cellular Automata in ExaHyPE
 * an Experiment.
 * 
 * Mainly based on Omohundro 1984: Modelling Cellular Automata
 * with Partial Differential Equations
 * http://fab.cba.mit.edu/classes/862.16/notes/computation/Omohundro-1984.pdf
 * 
 * Also interesting:
 * https://arxiv.org/pdf/1111.1567.pdf
 * https://arxiv.org/pdf/1003.1983.pdf
 * 
 * This application implements a purely nonconservative scheme which shall
 * mimic cellular automata in continuum limit.
 * 
 * Current status: It's not really working.
 * 
 **/

#include <iostream>

#include <algorithm> // std::fill
#include "kernels/KernelUtils.h" // kernels::idx2

#include "bfuncs.h"

#include "tarch/la/Vector.h"
typedef tarch::la::Vector<DIMENSIONS,int> vint;
typedef tarch::la::Vector<DIMENSIONS,double> vdouble;

// define all the aux stuff
const double c1=10; // Factor of Hysteresis
const double a1=0.00198742; // Integrate[b4[t], {t, -0.11, 0.11}] for shift speed
const double c2=0.5; // blow up the neighbourhood
const double c3=0.1;

// This assumes timestep size 1.0 and breaks a timestep into
// four intervals 0.20 each.
double t1 = 0.11;
double t2 = 0.36;
double t3 = 0.61;
double t4 = 0.86;

int T(const double* const _Q) {
	// Transition law T(N, S_i)
	const CA::CellularAutomata_FV::ReadOnlyVariables Q(_Q);
	using namespace tarch::la;
	
	// implement the rules of Conways Game of Life
	double neighbors = std::round(sum(Q.S()));
	if((0 <= neighbors && neighbors <= 1) || 4 <= neighbors) {
		return 0; // die
	} else if(equals(neighbors,3)) {
		return 1; // get born
	} else if(equals(neighbors,2)) {
		return Q.N(); // Nothing happens
	}
	// This should never happen
	return Q.N();
}

double H(const double x, const int c) {
	// Hysteresis to drive x towards c
	// Simple for GoL where only two states n=0 and n=1 are possible.
	// So we can drive with a harmonic oscillator
	return c1*SQ(x-c);
}

// round for the cell
vint cell(const double x, const double y) {
	vint ret;
	ret(0) = std::floor(x);
	ret(1) = std::floor(y);
	return ret;
}


// convenience
vint cell(vdouble& pos) {
	return cell(pos(0),pos(1));
}

// convenience
vint cell(const double* const _Q) {
	const CA::CellularAutomata_FV::ReadOnlyVariables Q(_Q);
	return cell(Q.x(),Q.y());
}



int id_for_neighbors(int id[3][3], int i, int j) {
	if(i>=0 && i<3 && j>=0 && j<3) return id[i][j];
	else return 0; // Boundary value: 0
}

tarch::logging::Log CA::CellularAutomata_FV::_log( "CA::CellularAutomata_FV" );


void CA::CellularAutomata_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::Parser::ParserView& constants) {
  // Draw some gnuplot files to proof the correctness of b6.
}

bool CA::CellularAutomata_FV::useAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) const {
  return true; // tarch::la::equals(t,0.0);
}

void CA::CellularAutomata_FV::adjustSolution(const double* const x,const double w,const double t,const double dt, double* _Q) {
  // Dimensions             = 2
  // Number of variables    = 10 + #parameters
  using namespace tarch::la;
  Variables Q(_Q);

  if(tarch::la::equals(t,0.0)) {
	// Initial conditions: random points
	vint pos = cell(x[0],x[1]);
	std::fill_n(_Q,NumberOfVariables+NumberOfParameters, 0);
	
	// Set the simplest possible data: a blinker
	// as it oscillates with frequency 2, we know the old=future data.
	int id[3][3] = {
		{0,1,0},
		{0,1,0},
		{0,1,0}
	};
	int old[3][3] = {
		{0,0,0},
		{1,1,1},
		{0,0,0}
	};

	int i=pos(0), j=pos(1);
	if(i<3 && j<3) {
		Q.N() = id[i][j];
		Q.F() = old[i][j];
		Q.S(0) = id_for_neighbors(id, i-1, j);
		Q.S(1) = id_for_neighbors(id, i-1, j+1);
		Q.S(2) = id_for_neighbors(id, i,   j+1);
		Q.S(3) = id_for_neighbors(id, i+1, j+1);
		Q.S(4) = id_for_neighbors(id, i+1, j);
		Q.S(5) = id_for_neighbors(id, i+1, j-1);
		Q.S(6) = id_for_neighbors(id, i,   j-1);
		Q.S(7) = id_for_neighbors(id, i-1, j-1);
	} else {
		// keep vacuum
	}

	// Set automata cell coordinates
	Q.x() = pos(0);
	Q.y() = pos(1);
	
	Q.dt() = 0;
  }

  // inject time
  Q.t() = t;
  Q.dt() += 1.0; // play around with timesteps
  // std::cout << "timestep: " << Q.dt() << std::endl;
}

struct stepguard {
	int s;
	stepguard(const double* const _Q) {
		const CA::CellularAutomata_FV::ReadOnlyVariables Q(_Q);
		s = (int)std::fmod(Q.dt(), 4.0);
	}
	double operator()(int wanted_step) {
		return (s == wanted_step ? 1.0 : 0.0);
	}
};


void CA::CellularAutomata_FV::coefficientMatrix(const double* const _Q,const int d,double* Bn) {
	const ReadOnlyVariables Q(_Q);
	kernels::idx2 idx_Bn(NumberOfVariables+NumberOfParameters, NumberOfVariables+NumberOfParameters);
	std::fill_n(Bn, idx_Bn.size, 0.0);
	
	const double coord = (d==0 ? Q.x() : Q.y());
	const double t = Q.t();
	
	// positive sign shifts f one space left.
	// Counting the S(0)..S(8) again from left around the center, like
	/********** 
	 *2  3  4 *
	 *1     5 *
	 *8  7  6 *
	 **********/
	static const int derivSign[DIMENSIONS][8] = {
		{+1,+1, 0,-1,-1,-1, 0,+1}, // x direciton
		{ 0,-1,-1,-1, 0,+1,+1,+1}, // y direction
	};

	/*// Literature
	Bn[idx_Bn(0,0)] = -b6(t-t1)*c2*db8(coord);
	for(int i=0; i<8; i++) { // wrong indexing
		Bn[idx_Bn(i,i)] = +b6(t-t3)/a1*derivSign[d][i];
	}
	*/
	
	// adhoc try
	stepguard tstep(_Q);
	Bn[idx_Bn(0,0)] = -tstep(1)*c2*db8(coord);
	for(int i=0; i<8; i++) {
		Bn[idx_Bn(2+i,2+i)] = +tstep(3)/a1*derivSign[d][i];
	}
}

void CA::CellularAutomata_FV::source(const double* const _Q,double* _S) {
	const ReadOnlyVariables Q(_Q);
	Variables Source(_S);
	
	const double t=Q.t(), x=Q.x(), y=Q.y();
	
	/* // Literature
	Source.F()  = -b6(t-t1) * c3 * Q.F() + b6(t-t4) * H(Q.F(), T(_Q));
	Source.N()  = -b6(t-t1) * c3 * Q.N() + b6(t-t2) * c3 * (b6(x)*b6(y)*Q.F() - Q.N());
	for(int i=0; i<8; i++) {
		Source.S(i) = -b6(t-t1) * c3 * Q.S(i) + b6(t-t2)*c3*(b6(x)*b6(y)*Q.F()-Q.S(i));
	}
	*/
	
	// adhoc try
	stepguard tstep(_Q);
	Source.F()  = -tstep(1) * c3 * Q.F() + tstep(4) * H(Q.F(), T(_Q));
	Source.N()  = -tstep(1) * c3 * Q.N() + tstep(2) * c3 * (b6(x)*b6(y)*Q.F() - Q.N());
	for(int i=0; i<8; i++) {
		Source.S(i) = -tstep(1) * c3 * Q.S(i) + tstep(2)*c3*(b6(x)*b6(y)*Q.F()-Q.S(i));
	}
}

exahype::solvers::Solver::RefinementControl CA::CellularAutomata_FV::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void CA::CellularAutomata_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 10 + #parameters
  std::fill_n(lambda, NumberOfVariables, 0.0);

  // we must make sure that simulation dt is much smaller than 0.2
  lambda[0] = 15.0;
}

void CA::CellularAutomata_FV::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 10 + #parameters
  
  for(int i=0; i<DIMENSIONS; i++)
    std::fill_n(F[i], NumberOfVariables, 0.0);
}


void CA::CellularAutomata_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 10 + #parameters

  // zero boundary conditions (death BC)
  std::fill_n(stateOutside, NumberOfVariables, 0.0);
}

