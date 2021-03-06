// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ErrorPlotter.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iomanip>

using std::cout; using std::cerr;
using std::endl; using std::string;
using std::ifstream; using std::ostringstream;
using std::istringstream;

ifstream input_file;

double maxRelError;

NavierStokes::ErrorPlotter::ErrorPlotter(NavierStokes::NavierStokesSolver_ADERDG& solver) {
  //reference file (generated with the same plotter) to compare numerical solution against
  input_file.open("/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/ApplicationExamples/CompressibleNavierStokes/ref.csv", std::ifstream::in);
}

NavierStokes::ErrorPlotter::~ErrorPlotter() {
  input_file.close();
}

void NavierStokes::ErrorPlotter::startPlotting( double time) {
  maxRelError = 0;
}


void NavierStokes::ErrorPlotter::finishPlotting() {
  std::cout<<"MAX REL ERROR "<<maxRelError<<std::endl;
}

void NavierStokes::ErrorPlotter::dissectOutputLine(std::string line, double *t, double *x, double *Q) {
  size_t pos = 0;
  std::string tokens[1+DIMENSIONS+5];
  
  //std::cout<<"parsing "<<line<<std::endl;
   
  for(int i=0; i<DIMENSIONS+6;i++) {
    pos = line.find("\t");
    tokens[i] = line.substr(0, pos);
    line.erase(0, pos + 1);
  }

  *t = std::stod(tokens[0]);

  for(int i=0; i<DIMENSIONS; i++) {
    x[i] = std::stod(tokens[i+1]);
  }
  for(int i=0; i<5; i++) {
    Q[i] = std::stod(tokens[i+DIMENSIONS]);
  }

}

/**
* Todo: this error checking may not be threadsafe, if ExaHyPE/Peano calls mapQuantities from multiple threads
**/
void NavierStokes::ErrorPlotter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 5;
  std::string line = "";  
  getline(input_file,line);
  //std::cout<<line<<std::endl;
  
  double t_in = 0;
  double Q_in[5];
  double x_in[DIMENSIONS];
   
  dissectOutputLine(line, &t_in, x_in, Q_in);

  // the conversion to str in the reference file must be applied to current state, too 
  //(such that input data and current cell state data are rounded similarly)
  std::string separator = "\t";
  std::stringstream stateString;
  stateString<<std::setprecision(16)<<timeStamp;
  for(int d = 0;d<DIMENSIONS;d++) {
    stateString<<separator<<x[d];
  }
  for (int i=0; i<writtenUnknowns; i++){ 
    stateString<<separator<<Q[i];
  }
  stateString<<std::endl;

  double t_out = 0;
  double Q_out[5];
  double x_out[DIMENSIONS];

  dissectOutputLine(stateString.str(), &t_out, x_out, Q_out);

  for(int d = 0;d<DIMENSIONS;d++) {
    assert(x_out[d]==x_in[d]);
  } 
  
  for (int i=0; i<writtenUnknowns; i++){ 
    if(Q_out[i]!=Q_in[i]) {
       double relError = std::abs(Q_out[i]-Q_in[i])/std::abs(Q_in[i]);
       maxRelError = std::max(relError, maxRelError);
    }
  }
   
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i]=Q[i];
  }
}
