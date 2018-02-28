#include "Demonstrator_FV.h"

#include "Demonstrator_FV_Variables.h"


tarch::logging::Log Qt::Demonstrator_FV::_log( "Qt::Demonstrator_FV" );

void Qt::Demonstrator_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // @todo Please implement/augment if required
}

void Qt::Demonstrator_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  // Dimensions             = 3
  // Number of variables    = 60 + #parameters
  
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
}

void Qt::Demonstrator_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 60 + #parameters
  
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
}

void Qt::Demonstrator_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 3
  // Number of variables    = 60 + #parameters

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
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit




//You can either implement this method or modify fusedSource
void Qt::Demonstrator_FV::algebraicSource(const double* const Q,double* S) {
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
}

void  Qt::Demonstrator_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
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
}

