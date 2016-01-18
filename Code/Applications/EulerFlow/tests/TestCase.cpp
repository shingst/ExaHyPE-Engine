#include "EulerFlow/tests/TestCase.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"
#include "EulerFlow/dg/ADERDG.h"
#include "EulerFlow/Constants.h"

#include <iostream>
using std::cout;
using std::endl;

registerTest(exahype::tests::TestCase)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif


exahype::tests::TestCase::TestCase():
  tarch::tests::TestCase( "exahype::tests::TestCase" ) {
}


exahype::tests::TestCase::~TestCase() {
}


void exahype::tests::TestCase::run() {
  // @todo If you have further tests, add them here
  //testMethod( test1 );

  testMethod( testSpaceTimePredictor );
  testMethod( testVolumeIntegral     );
  testMethod( testRiemannSolver      );
  testMethod( testSurfaceIntegral    );
  testMethod( testUpdateSolution     );

  exit(0);
}

/*
void exahype::tests::TestCase::test1() {
  // @todo Add your test here
  validateEquals(1,1);




  int trianglesA = 1;
  double xA[]    = {0.0, 0.0, 0.0};
  double yA[]    = {0.0, 1.0, 0.0};
  double zA[]    = {0.0, 0.0, 1.0};

  int trianglesB = 1;
  double xB[]    = {0.2, 0.2, 0.2};
  double yB[]    = {0.0, 1.0, 0.0};
  double zB[]    = {0.0, 0.0, 1.0};

  int numberOfContactPoints;

  double xC[]    = {-1.0, -1.0, -1.0};
  double yC[]    = {-1.0, -1.0, -1.0};
  double zC[]    = {-1.0, -1.0, -1.0};

  double xN[]    = {-1.0, -1.0, -1.0};
  double yN[]    = {-1.0, -1.0, -1.0};
  double zN[]    = {-1.0, -1.0, -1.0};

  ttd::ContactPoints::computeCollisionPoints(
    trianglesA, xA, yA, zA,
    trianglesB, xB, yB, zB,
    0.1,
    numberOfContactPoints,
    xC, yC, zC,
    xN, yN, zN
  );

  validateEquals(numberOfContactPoints,0);

  ttd::ContactPoints::computeCollisionPoints(
    trianglesA, xA, yA, zA,
    trianglesB, xB, yB, zB,
    0.25,
    numberOfContactPoints,
    xC, yC, zC,
    xN, yN, zN
  );

  validateEquals(numberOfContactPoints,1);

  validateNumericalEqualsWithParams3(xC[0],0.1,xC[0],yC[0],zC[0]);
  validateNumericalEqualsWithParams3(yC[0],0.0,xC[0],yC[0],zC[0]);
  validateNumericalEqualsWithParams3(zC[0],0.0,xC[0],yC[0],zC[0]);

  validateNumericalEqualsWithParams3(xN[0],0.1,xN[0],yN[0],zN[0]);
  validateNumericalEqualsWithParams3(yN[0],0.0,xN[0],yN[0],zN[0]);
  validateNumericalEqualsWithParams3(zN[0],0.0,xN[0],yN[0],zN[0]);


}
*/

void exahype::tests::TestCase::testSpaceTimePredictor() {

  cout << "Test space time predictor, ORDER=3, DIM=2" << endl;

  // input:
  double *luh = new double[80]();  // space DOF
  for(int i=0; i<16; i++) {
    luh[5*i+0] = 1.00000000000000000000e+00;
    luh[5*i+4] = 2.50000000000000044409e+00;
  }
  const double dx[2]        = {3.70370370370370349811e-02, 3.70370370370370349811e-02};
  const double timeStepSize =  1.40831757919882352703e-03;

  // local:
  double *lQi = new double[320]; // space-time DOF
  double *lFi = new double[640];
  double *rhs0 = new double[320];
  double *rhs  = new double[320];
  double *tmp  = new double[20];

  // output:
  double *lQhi = new double[80];
  double *lFhi = new double[160];
  double *lQhbnd = new double[80];
  double *lFhbnd = new double[80];

  dg::spaceTimePredictor<2>(lQi, lFi, luh, lQhi, lFhi, lQhbnd, lFhbnd, rhs0, rhs, tmp, dx, timeStepSize);

  //validateNumericalEquals(lQhi[i],<referencesolution>);
  validateNumericalEqualsWithEps(lQhi[0], 1, eps);
  validateNumericalEqualsWithEps(lQhi[1], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[2], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[3], 0, eps);
  validateNumericalEqualsWithEps(lQhi[4], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[5], 1, eps);
  validateNumericalEqualsWithEps(lQhi[6], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[7], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[8], 0, eps);
  validateNumericalEqualsWithEps(lQhi[9], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[10], 1, eps);
  validateNumericalEqualsWithEps(lQhi[11], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[12], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[13], 0, eps);
  validateNumericalEqualsWithEps(lQhi[14], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[15], 1, eps);
  validateNumericalEqualsWithEps(lQhi[16], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[17], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[18], 0, eps);
  validateNumericalEqualsWithEps(lQhi[19], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[20], 1, eps);
  validateNumericalEqualsWithEps(lQhi[21], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[22], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[23], 0, eps);
  validateNumericalEqualsWithEps(lQhi[24], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[25], 1, eps);
  validateNumericalEqualsWithEps(lQhi[26], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[27], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[28], 0, eps);
  validateNumericalEqualsWithEps(lQhi[29], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[30], 1, eps);
  validateNumericalEqualsWithEps(lQhi[31], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[32], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[33], 0, eps);
  validateNumericalEqualsWithEps(lQhi[34], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[35], 1, eps);
  validateNumericalEqualsWithEps(lQhi[36], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[37], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[38], 0, eps);
  validateNumericalEqualsWithEps(lQhi[39], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[40], 1, eps);
  validateNumericalEqualsWithEps(lQhi[41], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lQhi[42], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[43], 0, eps);
  validateNumericalEqualsWithEps(lQhi[44], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[45], 1, eps);
  validateNumericalEqualsWithEps(lQhi[46], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lQhi[47], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[48], 0, eps);
  validateNumericalEqualsWithEps(lQhi[49], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[50], 1, eps);
  validateNumericalEqualsWithEps(lQhi[51], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[52], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lQhi[53], 0, eps);
  validateNumericalEqualsWithEps(lQhi[54], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[55], 1, eps);
  validateNumericalEqualsWithEps(lQhi[56], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[57], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[58], 0, eps);
  validateNumericalEqualsWithEps(lQhi[59], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[60], 1, eps);
  validateNumericalEqualsWithEps(lQhi[61], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lQhi[62], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[63], 0, eps);
  validateNumericalEqualsWithEps(lQhi[64], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[65], 1, eps);
  validateNumericalEqualsWithEps(lQhi[66], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lQhi[67], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[68], 0, eps);
  validateNumericalEqualsWithEps(lQhi[69], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[70], 1, eps);
  validateNumericalEqualsWithEps(lQhi[71], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lQhi[72], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lQhi[73], 0, eps);
  validateNumericalEqualsWithEps(lQhi[74], 2.5, eps);
  validateNumericalEqualsWithEps(lQhi[75], 1, eps);
  validateNumericalEqualsWithEps(lQhi[76], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[77], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lQhi[78], 0, eps);
  validateNumericalEqualsWithEps(lQhi[79], 2.5, eps);

  // lFhi_x
  validateNumericalEqualsWithEps(lFhi[0], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[1], 1, eps);
  validateNumericalEqualsWithEps(lFhi[2], 3.90655760696e-35, eps);
  validateNumericalEqualsWithEps(lFhi[3], 0, eps);
  validateNumericalEqualsWithEps(lFhi[4], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[5], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[6], 1, eps);
  validateNumericalEqualsWithEps(lFhi[7], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[8], 0, eps);
  validateNumericalEqualsWithEps(lFhi[9], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[10], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[11], 1, eps);
  validateNumericalEqualsWithEps(lFhi[12], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[13], 0, eps);
  validateNumericalEqualsWithEps(lFhi[14], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[15], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[16], 1, eps);
  validateNumericalEqualsWithEps(lFhi[17], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[18], 0, eps);
  validateNumericalEqualsWithEps(lFhi[19], 2.68302860327e-17, eps);
  validateNumericalEqualsWithEps(lFhi[20], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[21], 1, eps);
  validateNumericalEqualsWithEps(lFhi[22], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[23], 0, eps);
  validateNumericalEqualsWithEps(lFhi[24], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[25], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[26], 1, eps);
  validateNumericalEqualsWithEps(lFhi[27], 1.09880295341e-34, eps);
  validateNumericalEqualsWithEps(lFhi[28], 0, eps);
  validateNumericalEqualsWithEps(lFhi[29], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[30], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[31], 1, eps);
  validateNumericalEqualsWithEps(lFhi[32], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[33], 0, eps);
  validateNumericalEqualsWithEps(lFhi[34], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[35], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[36], 1, eps);
  validateNumericalEqualsWithEps(lFhi[37], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[38], 0, eps);
  validateNumericalEqualsWithEps(lFhi[39], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[40], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[41], 1, eps);
  validateNumericalEqualsWithEps(lFhi[42], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[43], 0, eps);
  validateNumericalEqualsWithEps(lFhi[44], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[45], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[46], 1, eps);
  validateNumericalEqualsWithEps(lFhi[47], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[48], 0, eps);
  validateNumericalEqualsWithEps(lFhi[49], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[50], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[51], 1, eps);
  validateNumericalEqualsWithEps(lFhi[52], 4.8794232613e-35, eps);
  validateNumericalEqualsWithEps(lFhi[53], 0, eps);
  validateNumericalEqualsWithEps(lFhi[54], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[55], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[56], 1, eps);
  validateNumericalEqualsWithEps(lFhi[57], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[58], 0, eps);
  validateNumericalEqualsWithEps(lFhi[59], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[60], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[61], 1, eps);
  validateNumericalEqualsWithEps(lFhi[62], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[63], 0, eps);
  validateNumericalEqualsWithEps(lFhi[64], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[65], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[66], 1, eps);
  validateNumericalEqualsWithEps(lFhi[67], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[68], 0, eps);
  validateNumericalEqualsWithEps(lFhi[69], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[70], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[71], 1, eps);
  validateNumericalEqualsWithEps(lFhi[72], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[73], 0, eps);
  validateNumericalEqualsWithEps(lFhi[74], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[75], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[76], 1, eps);
  validateNumericalEqualsWithEps(lFhi[77], 7.37687416339e-35, eps);
  validateNumericalEqualsWithEps(lFhi[78], 0, eps);
  validateNumericalEqualsWithEps(lFhi[79], 2.68302860327e-17, eps);
  // lFhi_y
  validateNumericalEqualsWithEps(lFhi[80], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[81], 3.90655760696e-35, eps);
  validateNumericalEqualsWithEps(lFhi[82], 1, eps);
  validateNumericalEqualsWithEps(lFhi[83], 0, eps);
  validateNumericalEqualsWithEps(lFhi[84], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[85], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[86], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[87], 1, eps);
  validateNumericalEqualsWithEps(lFhi[88], 0, eps);
  validateNumericalEqualsWithEps(lFhi[89], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[90], -1.48132747725e-17, eps);
  validateNumericalEqualsWithEps(lFhi[91], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[92], 1, eps);
  validateNumericalEqualsWithEps(lFhi[93], 0, eps);
  validateNumericalEqualsWithEps(lFhi[94], -5.18464617038e-17, eps);
  validateNumericalEqualsWithEps(lFhi[95], -5.78725778411e-18, eps);
  validateNumericalEqualsWithEps(lFhi[96], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[97], 1, eps);
  validateNumericalEqualsWithEps(lFhi[98], 0, eps);
  validateNumericalEqualsWithEps(lFhi[99], -2.02554022444e-17, eps);
  validateNumericalEqualsWithEps(lFhi[100], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[101], -6.28505465062e-35, eps);
  validateNumericalEqualsWithEps(lFhi[102], 1, eps);
  validateNumericalEqualsWithEps(lFhi[103], 0, eps);
  validateNumericalEqualsWithEps(lFhi[104], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[105], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[106], 1.09880295341e-34, eps);
  validateNumericalEqualsWithEps(lFhi[107], 1, eps);
  validateNumericalEqualsWithEps(lFhi[108], 0, eps);
  validateNumericalEqualsWithEps(lFhi[109], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[110], 7.63123721751e-18, eps);
  validateNumericalEqualsWithEps(lFhi[111], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[112], 1, eps);
  validateNumericalEqualsWithEps(lFhi[113], 0, eps);
  validateNumericalEqualsWithEps(lFhi[114], 2.67093302613e-17, eps);
  validateNumericalEqualsWithEps(lFhi[115], 3.32798659798e-18, eps);
  validateNumericalEqualsWithEps(lFhi[116], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[117], 1, eps);
  validateNumericalEqualsWithEps(lFhi[118], 0, eps);
  validateNumericalEqualsWithEps(lFhi[119], 1.16479530929e-17, eps);
  validateNumericalEqualsWithEps(lFhi[120], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[121], 6.02226203114e-35, eps);
  validateNumericalEqualsWithEps(lFhi[122], 1, eps);
  validateNumericalEqualsWithEps(lFhi[123], 0, eps);
  validateNumericalEqualsWithEps(lFhi[124], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[125], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[126], -7.10324058069e-35, eps);
  validateNumericalEqualsWithEps(lFhi[127], 1, eps);
  validateNumericalEqualsWithEps(lFhi[128], 0, eps);
  validateNumericalEqualsWithEps(lFhi[129], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[130], -4.17981679914e-18, eps);
  validateNumericalEqualsWithEps(lFhi[131], 4.8794232613e-35, eps);
  validateNumericalEqualsWithEps(lFhi[132], 1, eps);
  validateNumericalEqualsWithEps(lFhi[133], 0, eps);
  validateNumericalEqualsWithEps(lFhi[134], -1.4629358797e-17, eps);
  validateNumericalEqualsWithEps(lFhi[135], -3.07890876627e-18, eps);
  validateNumericalEqualsWithEps(lFhi[136], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[137], 1, eps);
  validateNumericalEqualsWithEps(lFhi[138], 0, eps);
  validateNumericalEqualsWithEps(lFhi[139], -1.07761806819e-17, eps);
  validateNumericalEqualsWithEps(lFhi[140], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[141], -5.27718213855e-35, eps);
  validateNumericalEqualsWithEps(lFhi[142], 1, eps);
  validateNumericalEqualsWithEps(lFhi[143], 0, eps);
  validateNumericalEqualsWithEps(lFhi[144], 2.68302860327e-17, eps);
  validateNumericalEqualsWithEps(lFhi[145], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[146], 5.00229224796e-35, eps);
  validateNumericalEqualsWithEps(lFhi[147], 1, eps);
  validateNumericalEqualsWithEps(lFhi[148], 0, eps);
  validateNumericalEqualsWithEps(lFhi[149], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[150], 1.23156350651e-17, eps);
  validateNumericalEqualsWithEps(lFhi[151], -4.73179353928e-35, eps);
  validateNumericalEqualsWithEps(lFhi[152], 1, eps);
  validateNumericalEqualsWithEps(lFhi[153], 0, eps);
  validateNumericalEqualsWithEps(lFhi[154], 4.31047227278e-17, eps);
  validateNumericalEqualsWithEps(lFhi[155], 7.66579600935e-18, eps);
  validateNumericalEqualsWithEps(lFhi[156], 7.37687416339e-35, eps);
  validateNumericalEqualsWithEps(lFhi[157], 1, eps);
  validateNumericalEqualsWithEps(lFhi[158], 0, eps);
  validateNumericalEqualsWithEps(lFhi[159], 2.68302860327e-17, eps);


  // TODO
  // lQhbnd, lQhbnd

  delete[] luh;
  delete[] lQi;
  delete[] lFi;
  delete[] rhs0;
  delete[] rhs;
  delete[] tmp;
  delete[] lQhi;
  delete[] lFhi;
  delete[] lQhbnd;
  delete[] lFhbnd;

} // testSpaceTimePredictor


void exahype::tests::TestCase::testVolumeIntegral() {


} // testVolumeIntegral


void exahype::tests::TestCase::testRiemannSolver() {

} // testRiemannSolver


void exahype::tests::TestCase::testSurfaceIntegral() {

} // testSurfaceIntegral



void exahype::tests::TestCase::testUpdateSolution() {

} // testUpdateSolution


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif

