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
  input_file.open("/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/ApplicationExamples/CompressibleNavierStokes/ref.csv", std::ifstream::in);
}

NavierStokes::ErrorPlotter::~ErrorPlotter() {
  input_file.close();
}

void NavierStokes::ErrorPlotter::startPlotting( double time) {
  // @TODO Please insert your code here.
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

  std::string separator = "\t";
  std::stringstream output;
  output<<std::setprecision(16)<<timeStamp;
  for(int d = 0;d<DIMENSIONS;d++) {
    output<<separator<<x[d];
  }
  for (int i=0; i<writtenUnknowns; i++){ 
    output<<separator<<Q[i];
  }
  output<<std::endl;

  double t_out = 0;
  double Q_out[5];
  double x_out[DIMENSIONS];

  dissectOutputLine(output.str(), &t_out, x_out, Q_out);

  //if(t_in!=t_out) {
  //if(output.str()!=line) {
  //  std::cout<<"output "<<output.str()<<std::endl;
  //  std::cout<<"line "<<line<<std::endl;
  //}
  //}

  //assert(t_in==t_out);
  for(int d = 0;d<DIMENSIONS;d++) {
    assert(x_out[d]==x_in[d]);
  } 
  //if(output.str()!=line) {
  //  std::cout<<"output "<<output.str()<<std::endl;
  //  std::cout<<"line "<<line<<std::endl;
  //}

  
  for (int i=0; i<writtenUnknowns; i++){ 
    if(Q_out[i]!=Q_in[i]) {
       double relError = std::abs(Q_out[i]-Q_in[i])/std::abs(Q_in[i]);
       maxRelError = std::max(relError, maxRelError);
    }
  }

  //assert(output.str()==line);
   
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i]=Q[i];
  }
}
