#include "FOCCZ4Solver_FV.h"

#include "FOCCZ4Solver_FV_Variables.h"


tarch::logging::Log FOCCZ4::FOCCZ4Solver_FV::_log( "FOCCZ4::FOCCZ4Solver_FV" );

void FOCCZ4::FOCCZ4Solver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  
  // @todo Please implement/augment if required
}

void FOCCZ4::FOCCZ4Solver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
  Q[7] = 0.0;
  Q[8] = 0.0;
  Q[9] = 0.0;
  Q[10] = 0.0;
  Q[11] = 0.0;
  Q[12] = 0.0;
  Q[13] = 0.0;
  Q[14] = 0.0;
  Q[15] = 0.0;
  Q[16] = 0.0;
  Q[17] = 0.0;
  Q[18] = 0.0;
  Q[19] = 0.0;
  Q[20] = 0.0;
  Q[21] = 0.0;
  Q[22] = 0.0;
  Q[23] = 0.0;
  Q[24] = 0.0;
  Q[25] = 0.0;
  Q[26] = 0.0;
  Q[27] = 0.0;
  Q[28] = 0.0;
  Q[29] = 0.0;
  Q[30] = 0.0;
  Q[31] = 0.0;
  Q[32] = 0.0;
  Q[33] = 0.0;
  Q[34] = 0.0;
  Q[35] = 0.0;
  Q[36] = 0.0;
  Q[37] = 0.0;
  Q[38] = 0.0;
  Q[39] = 0.0;
  Q[40] = 0.0;
  Q[41] = 0.0;
  Q[42] = 0.0;
  Q[43] = 0.0;
  Q[44] = 0.0;
  Q[45] = 0.0;
  Q[46] = 0.0;
  Q[47] = 0.0;
  Q[48] = 0.0;
  Q[49] = 0.0;
  Q[50] = 0.0;
  Q[51] = 0.0;
  Q[52] = 0.0;
  Q[53] = 0.0;
  Q[54] = 0.0;
  Q[55] = 0.0;
  Q[56] = 0.0;
  Q[57] = 0.0;
  Q[58] = 0.0;
  Q[59] = 0.0;
  Q[60] = 0.0;
  Q[61] = 0.0;
  Q[62] = 0.0;
  Q[63] = 0.0;
  Q[64] = 0.0;
  Q[65] = 0.0;
  Q[66] = 0.0;
  Q[67] = 0.0;
  Q[68] = 0.0;
  Q[69] = 0.0;
  Q[70] = 0.0;
  Q[71] = 0.0;
  Q[72] = 0.0;
  Q[73] = 0.0;
  Q[74] = 0.0;
  Q[75] = 0.0;
  Q[76] = 0.0;
  Q[77] = 0.0;
  Q[78] = 0.0;
  Q[79] = 0.0;
  Q[80] = 0.0;
  Q[81] = 0.0;
  Q[82] = 0.0;
  Q[83] = 0.0;
  Q[84] = 0.0;
  Q[85] = 0.0;
  Q[86] = 0.0;
  Q[87] = 0.0;
  Q[88] = 0.0;
  Q[89] = 0.0;
  Q[90] = 0.0;
  Q[91] = 0.0;
  Q[92] = 0.0;
  Q[93] = 0.0;
  Q[94] = 0.0;
  Q[95] = 0.0;
}

void FOCCZ4::FOCCZ4Solver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
  lambda[5] = 1.0;
  lambda[6] = 1.0;
  lambda[7] = 1.0;
  lambda[8] = 1.0;
  lambda[9] = 1.0;
  lambda[10] = 1.0;
  lambda[11] = 1.0;
  lambda[12] = 1.0;
  lambda[13] = 1.0;
  lambda[14] = 1.0;
  lambda[15] = 1.0;
  lambda[16] = 1.0;
  lambda[17] = 1.0;
  lambda[18] = 1.0;
  lambda[19] = 1.0;
  lambda[20] = 1.0;
  lambda[21] = 1.0;
  lambda[22] = 1.0;
  lambda[23] = 1.0;
  lambda[24] = 1.0;
  lambda[25] = 1.0;
  lambda[26] = 1.0;
  lambda[27] = 1.0;
  lambda[28] = 1.0;
  lambda[29] = 1.0;
  lambda[30] = 1.0;
  lambda[31] = 1.0;
  lambda[32] = 1.0;
  lambda[33] = 1.0;
  lambda[34] = 1.0;
  lambda[35] = 1.0;
  lambda[36] = 1.0;
  lambda[37] = 1.0;
  lambda[38] = 1.0;
  lambda[39] = 1.0;
  lambda[40] = 1.0;
  lambda[41] = 1.0;
  lambda[42] = 1.0;
  lambda[43] = 1.0;
  lambda[44] = 1.0;
  lambda[45] = 1.0;
  lambda[46] = 1.0;
  lambda[47] = 1.0;
  lambda[48] = 1.0;
  lambda[49] = 1.0;
  lambda[50] = 1.0;
  lambda[51] = 1.0;
  lambda[52] = 1.0;
  lambda[53] = 1.0;
  lambda[54] = 1.0;
  lambda[55] = 1.0;
  lambda[56] = 1.0;
  lambda[57] = 1.0;
  lambda[58] = 1.0;
  lambda[59] = 1.0;
  lambda[60] = 1.0;
  lambda[61] = 1.0;
  lambda[62] = 1.0;
  lambda[63] = 1.0;
  lambda[64] = 1.0;
  lambda[65] = 1.0;
  lambda[66] = 1.0;
  lambda[67] = 1.0;
  lambda[68] = 1.0;
  lambda[69] = 1.0;
  lambda[70] = 1.0;
  lambda[71] = 1.0;
  lambda[72] = 1.0;
  lambda[73] = 1.0;
  lambda[74] = 1.0;
  lambda[75] = 1.0;
  lambda[76] = 1.0;
  lambda[77] = 1.0;
  lambda[78] = 1.0;
  lambda[79] = 1.0;
  lambda[80] = 1.0;
  lambda[81] = 1.0;
  lambda[82] = 1.0;
  lambda[83] = 1.0;
  lambda[84] = 1.0;
  lambda[85] = 1.0;
  lambda[86] = 1.0;
  lambda[87] = 1.0;
  lambda[88] = 1.0;
  lambda[89] = 1.0;
  lambda[90] = 1.0;
  lambda[91] = 1.0;
  lambda[92] = 1.0;
  lambda[93] = 1.0;
  lambda[94] = 1.0;
  lambda[95] = 1.0;
}

void FOCCZ4::FOCCZ4Solver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // @todo Please implement/augment if required
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];
  stateOutside[4] = stateInside[4];
  stateOutside[5] = stateInside[5];
  stateOutside[6] = stateInside[6];
  stateOutside[7] = stateInside[7];
  stateOutside[8] = stateInside[8];
  stateOutside[9] = stateInside[9];
  stateOutside[10] = stateInside[10];
  stateOutside[11] = stateInside[11];
  stateOutside[12] = stateInside[12];
  stateOutside[13] = stateInside[13];
  stateOutside[14] = stateInside[14];
  stateOutside[15] = stateInside[15];
  stateOutside[16] = stateInside[16];
  stateOutside[17] = stateInside[17];
  stateOutside[18] = stateInside[18];
  stateOutside[19] = stateInside[19];
  stateOutside[20] = stateInside[20];
  stateOutside[21] = stateInside[21];
  stateOutside[22] = stateInside[22];
  stateOutside[23] = stateInside[23];
  stateOutside[24] = stateInside[24];
  stateOutside[25] = stateInside[25];
  stateOutside[26] = stateInside[26];
  stateOutside[27] = stateInside[27];
  stateOutside[28] = stateInside[28];
  stateOutside[29] = stateInside[29];
  stateOutside[30] = stateInside[30];
  stateOutside[31] = stateInside[31];
  stateOutside[32] = stateInside[32];
  stateOutside[33] = stateInside[33];
  stateOutside[34] = stateInside[34];
  stateOutside[35] = stateInside[35];
  stateOutside[36] = stateInside[36];
  stateOutside[37] = stateInside[37];
  stateOutside[38] = stateInside[38];
  stateOutside[39] = stateInside[39];
  stateOutside[40] = stateInside[40];
  stateOutside[41] = stateInside[41];
  stateOutside[42] = stateInside[42];
  stateOutside[43] = stateInside[43];
  stateOutside[44] = stateInside[44];
  stateOutside[45] = stateInside[45];
  stateOutside[46] = stateInside[46];
  stateOutside[47] = stateInside[47];
  stateOutside[48] = stateInside[48];
  stateOutside[49] = stateInside[49];
  stateOutside[50] = stateInside[50];
  stateOutside[51] = stateInside[51];
  stateOutside[52] = stateInside[52];
  stateOutside[53] = stateInside[53];
  stateOutside[54] = stateInside[54];
  stateOutside[55] = stateInside[55];
  stateOutside[56] = stateInside[56];
  stateOutside[57] = stateInside[57];
  stateOutside[58] = stateInside[58];
  stateOutside[59] = stateInside[59];
  stateOutside[60] = stateInside[60];
  stateOutside[61] = stateInside[61];
  stateOutside[62] = stateInside[62];
  stateOutside[63] = stateInside[63];
  stateOutside[64] = stateInside[64];
  stateOutside[65] = stateInside[65];
  stateOutside[66] = stateInside[66];
  stateOutside[67] = stateInside[67];
  stateOutside[68] = stateInside[68];
  stateOutside[69] = stateInside[69];
  stateOutside[70] = stateInside[70];
  stateOutside[71] = stateInside[71];
  stateOutside[72] = stateInside[72];
  stateOutside[73] = stateInside[73];
  stateOutside[74] = stateInside[74];
  stateOutside[75] = stateInside[75];
  stateOutside[76] = stateInside[76];
  stateOutside[77] = stateInside[77];
  stateOutside[78] = stateInside[78];
  stateOutside[79] = stateInside[79];
  stateOutside[80] = stateInside[80];
  stateOutside[81] = stateInside[81];
  stateOutside[82] = stateInside[82];
  stateOutside[83] = stateInside[83];
  stateOutside[84] = stateInside[84];
  stateOutside[85] = stateInside[85];
  stateOutside[86] = stateInside[86];
  stateOutside[87] = stateInside[87];
  stateOutside[88] = stateInside[88];
  stateOutside[89] = stateInside[89];
  stateOutside[90] = stateInside[90];
  stateOutside[91] = stateInside[91];
  stateOutside[92] = stateInside[92];
  stateOutside[93] = stateInside[93];
  stateOutside[94] = stateInside[94];
  stateOutside[95] = stateInside[95];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void FOCCZ4::FOCCZ4Solver_FV::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;
  F[0][5] = 0.0;
  F[0][6] = 0.0;
  F[0][7] = 0.0;
  F[0][8] = 0.0;
  F[0][9] = 0.0;
  F[0][10] = 0.0;
  F[0][11] = 0.0;
  F[0][12] = 0.0;
  F[0][13] = 0.0;
  F[0][14] = 0.0;
  F[0][15] = 0.0;
  F[0][16] = 0.0;
  F[0][17] = 0.0;
  F[0][18] = 0.0;
  F[0][19] = 0.0;
  F[0][20] = 0.0;
  F[0][21] = 0.0;
  F[0][22] = 0.0;
  F[0][23] = 0.0;
  F[0][24] = 0.0;
  F[0][25] = 0.0;
  F[0][26] = 0.0;
  F[0][27] = 0.0;
  F[0][28] = 0.0;
  F[0][29] = 0.0;
  F[0][30] = 0.0;
  F[0][31] = 0.0;
  F[0][32] = 0.0;
  F[0][33] = 0.0;
  F[0][34] = 0.0;
  F[0][35] = 0.0;
  F[0][36] = 0.0;
  F[0][37] = 0.0;
  F[0][38] = 0.0;
  F[0][39] = 0.0;
  F[0][40] = 0.0;
  F[0][41] = 0.0;
  F[0][42] = 0.0;
  F[0][43] = 0.0;
  F[0][44] = 0.0;
  F[0][45] = 0.0;
  F[0][46] = 0.0;
  F[0][47] = 0.0;
  F[0][48] = 0.0;
  F[0][49] = 0.0;
  F[0][50] = 0.0;
  F[0][51] = 0.0;
  F[0][52] = 0.0;
  F[0][53] = 0.0;
  F[0][54] = 0.0;
  F[0][55] = 0.0;
  F[0][56] = 0.0;
  F[0][57] = 0.0;
  F[0][58] = 0.0;
  F[0][59] = 0.0;
  F[0][60] = 0.0;
  F[0][61] = 0.0;
  F[0][62] = 0.0;
  F[0][63] = 0.0;
  F[0][64] = 0.0;
  F[0][65] = 0.0;
  F[0][66] = 0.0;
  F[0][67] = 0.0;
  F[0][68] = 0.0;
  F[0][69] = 0.0;
  F[0][70] = 0.0;
  F[0][71] = 0.0;
  F[0][72] = 0.0;
  F[0][73] = 0.0;
  F[0][74] = 0.0;
  F[0][75] = 0.0;
  F[0][76] = 0.0;
  F[0][77] = 0.0;
  F[0][78] = 0.0;
  F[0][79] = 0.0;
  F[0][80] = 0.0;
  F[0][81] = 0.0;
  F[0][82] = 0.0;
  F[0][83] = 0.0;
  F[0][84] = 0.0;
  F[0][85] = 0.0;
  F[0][86] = 0.0;
  F[0][87] = 0.0;
  F[0][88] = 0.0;
  F[0][89] = 0.0;
  F[0][90] = 0.0;
  F[0][91] = 0.0;
  F[0][92] = 0.0;
  F[0][93] = 0.0;
  F[0][94] = 0.0;
  F[0][95] = 0.0;
  
  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
  F[1][5] = 0.0;
  F[1][6] = 0.0;
  F[1][7] = 0.0;
  F[1][8] = 0.0;
  F[1][9] = 0.0;
  F[1][10] = 0.0;
  F[1][11] = 0.0;
  F[1][12] = 0.0;
  F[1][13] = 0.0;
  F[1][14] = 0.0;
  F[1][15] = 0.0;
  F[1][16] = 0.0;
  F[1][17] = 0.0;
  F[1][18] = 0.0;
  F[1][19] = 0.0;
  F[1][20] = 0.0;
  F[1][21] = 0.0;
  F[1][22] = 0.0;
  F[1][23] = 0.0;
  F[1][24] = 0.0;
  F[1][25] = 0.0;
  F[1][26] = 0.0;
  F[1][27] = 0.0;
  F[1][28] = 0.0;
  F[1][29] = 0.0;
  F[1][30] = 0.0;
  F[1][31] = 0.0;
  F[1][32] = 0.0;
  F[1][33] = 0.0;
  F[1][34] = 0.0;
  F[1][35] = 0.0;
  F[1][36] = 0.0;
  F[1][37] = 0.0;
  F[1][38] = 0.0;
  F[1][39] = 0.0;
  F[1][40] = 0.0;
  F[1][41] = 0.0;
  F[1][42] = 0.0;
  F[1][43] = 0.0;
  F[1][44] = 0.0;
  F[1][45] = 0.0;
  F[1][46] = 0.0;
  F[1][47] = 0.0;
  F[1][48] = 0.0;
  F[1][49] = 0.0;
  F[1][50] = 0.0;
  F[1][51] = 0.0;
  F[1][52] = 0.0;
  F[1][53] = 0.0;
  F[1][54] = 0.0;
  F[1][55] = 0.0;
  F[1][56] = 0.0;
  F[1][57] = 0.0;
  F[1][58] = 0.0;
  F[1][59] = 0.0;
  F[1][60] = 0.0;
  F[1][61] = 0.0;
  F[1][62] = 0.0;
  F[1][63] = 0.0;
  F[1][64] = 0.0;
  F[1][65] = 0.0;
  F[1][66] = 0.0;
  F[1][67] = 0.0;
  F[1][68] = 0.0;
  F[1][69] = 0.0;
  F[1][70] = 0.0;
  F[1][71] = 0.0;
  F[1][72] = 0.0;
  F[1][73] = 0.0;
  F[1][74] = 0.0;
  F[1][75] = 0.0;
  F[1][76] = 0.0;
  F[1][77] = 0.0;
  F[1][78] = 0.0;
  F[1][79] = 0.0;
  F[1][80] = 0.0;
  F[1][81] = 0.0;
  F[1][82] = 0.0;
  F[1][83] = 0.0;
  F[1][84] = 0.0;
  F[1][85] = 0.0;
  F[1][86] = 0.0;
  F[1][87] = 0.0;
  F[1][88] = 0.0;
  F[1][89] = 0.0;
  F[1][90] = 0.0;
  F[1][91] = 0.0;
  F[1][92] = 0.0;
  F[1][93] = 0.0;
  F[1][94] = 0.0;
  F[1][95] = 0.0;
  
}




//You can either implement this method or modify fusedSource
void FOCCZ4::FOCCZ4Solver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  // @todo Please implement/augment if required
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = 0.0;
  S[6] = 0.0;
  S[7] = 0.0;
  S[8] = 0.0;
  S[9] = 0.0;
  S[10] = 0.0;
  S[11] = 0.0;
  S[12] = 0.0;
  S[13] = 0.0;
  S[14] = 0.0;
  S[15] = 0.0;
  S[16] = 0.0;
  S[17] = 0.0;
  S[18] = 0.0;
  S[19] = 0.0;
  S[20] = 0.0;
  S[21] = 0.0;
  S[22] = 0.0;
  S[23] = 0.0;
  S[24] = 0.0;
  S[25] = 0.0;
  S[26] = 0.0;
  S[27] = 0.0;
  S[28] = 0.0;
  S[29] = 0.0;
  S[30] = 0.0;
  S[31] = 0.0;
  S[32] = 0.0;
  S[33] = 0.0;
  S[34] = 0.0;
  S[35] = 0.0;
  S[36] = 0.0;
  S[37] = 0.0;
  S[38] = 0.0;
  S[39] = 0.0;
  S[40] = 0.0;
  S[41] = 0.0;
  S[42] = 0.0;
  S[43] = 0.0;
  S[44] = 0.0;
  S[45] = 0.0;
  S[46] = 0.0;
  S[47] = 0.0;
  S[48] = 0.0;
  S[49] = 0.0;
  S[50] = 0.0;
  S[51] = 0.0;
  S[52] = 0.0;
  S[53] = 0.0;
  S[54] = 0.0;
  S[55] = 0.0;
  S[56] = 0.0;
  S[57] = 0.0;
  S[58] = 0.0;
  S[59] = 0.0;
  S[60] = 0.0;
  S[61] = 0.0;
  S[62] = 0.0;
  S[63] = 0.0;
  S[64] = 0.0;
  S[65] = 0.0;
  S[66] = 0.0;
  S[67] = 0.0;
  S[68] = 0.0;
  S[69] = 0.0;
  S[70] = 0.0;
  S[71] = 0.0;
  S[72] = 0.0;
  S[73] = 0.0;
  S[74] = 0.0;
  S[75] = 0.0;
  S[76] = 0.0;
  S[77] = 0.0;
  S[78] = 0.0;
  S[79] = 0.0;
  S[80] = 0.0;
  S[81] = 0.0;
  S[82] = 0.0;
  S[83] = 0.0;
  S[84] = 0.0;
  S[85] = 0.0;
  S[86] = 0.0;
  S[87] = 0.0;
  S[88] = 0.0;
  S[89] = 0.0;
  S[90] = 0.0;
  S[91] = 0.0;
  S[92] = 0.0;
  S[93] = 0.0;
  S[94] = 0.0;
  S[95] = 0.0;
}

void  FOCCZ4::FOCCZ4Solver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver_FV.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  BgradQ[4] = 0.0;
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = 0.0;
  BgradQ[9] = 0.0;
  BgradQ[10] = 0.0;
  BgradQ[11] = 0.0;
  BgradQ[12] = 0.0;
  BgradQ[13] = 0.0;
  BgradQ[14] = 0.0;
  BgradQ[15] = 0.0;
  BgradQ[16] = 0.0;
  BgradQ[17] = 0.0;
  BgradQ[18] = 0.0;
  BgradQ[19] = 0.0;
  BgradQ[20] = 0.0;
  BgradQ[21] = 0.0;
  BgradQ[22] = 0.0;
  BgradQ[23] = 0.0;
  BgradQ[24] = 0.0;
  BgradQ[25] = 0.0;
  BgradQ[26] = 0.0;
  BgradQ[27] = 0.0;
  BgradQ[28] = 0.0;
  BgradQ[29] = 0.0;
  BgradQ[30] = 0.0;
  BgradQ[31] = 0.0;
  BgradQ[32] = 0.0;
  BgradQ[33] = 0.0;
  BgradQ[34] = 0.0;
  BgradQ[35] = 0.0;
  BgradQ[36] = 0.0;
  BgradQ[37] = 0.0;
  BgradQ[38] = 0.0;
  BgradQ[39] = 0.0;
  BgradQ[40] = 0.0;
  BgradQ[41] = 0.0;
  BgradQ[42] = 0.0;
  BgradQ[43] = 0.0;
  BgradQ[44] = 0.0;
  BgradQ[45] = 0.0;
  BgradQ[46] = 0.0;
  BgradQ[47] = 0.0;
  BgradQ[48] = 0.0;
  BgradQ[49] = 0.0;
  BgradQ[50] = 0.0;
  BgradQ[51] = 0.0;
  BgradQ[52] = 0.0;
  BgradQ[53] = 0.0;
  BgradQ[54] = 0.0;
  BgradQ[55] = 0.0;
  BgradQ[56] = 0.0;
  BgradQ[57] = 0.0;
  BgradQ[58] = 0.0;
  BgradQ[59] = 0.0;
  BgradQ[60] = 0.0;
  BgradQ[61] = 0.0;
  BgradQ[62] = 0.0;
  BgradQ[63] = 0.0;
  BgradQ[64] = 0.0;
  BgradQ[65] = 0.0;
  BgradQ[66] = 0.0;
  BgradQ[67] = 0.0;
  BgradQ[68] = 0.0;
  BgradQ[69] = 0.0;
  BgradQ[70] = 0.0;
  BgradQ[71] = 0.0;
  BgradQ[72] = 0.0;
  BgradQ[73] = 0.0;
  BgradQ[74] = 0.0;
  BgradQ[75] = 0.0;
  BgradQ[76] = 0.0;
  BgradQ[77] = 0.0;
  BgradQ[78] = 0.0;
  BgradQ[79] = 0.0;
  BgradQ[80] = 0.0;
  BgradQ[81] = 0.0;
  BgradQ[82] = 0.0;
  BgradQ[83] = 0.0;
  BgradQ[84] = 0.0;
  BgradQ[85] = 0.0;
  BgradQ[86] = 0.0;
  BgradQ[87] = 0.0;
  BgradQ[88] = 0.0;
  BgradQ[89] = 0.0;
  BgradQ[90] = 0.0;
  BgradQ[91] = 0.0;
  BgradQ[92] = 0.0;
  BgradQ[93] = 0.0;
  BgradQ[94] = 0.0;
  BgradQ[95] = 0.0;
}

