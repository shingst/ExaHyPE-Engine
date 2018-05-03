#include "IntegralsWriter.h"

#include "Fortran/PDE.h"
#include "InitialData/InitialData.h"
#include "CCZ4Solver_ADERDG.h"
#include "CCZ4Solver_FV.h"

#include "kernels/aderdg/generic/c/sizes.cpph"

#include <cmath>

CCZ4::IntegralsWriter::IntegralsWriter(exahype::solvers::LimitingADERDGSolver&  solver)
	: IntegralsWriter() { plotForADERSolver = true; }

CCZ4::IntegralsWriter::IntegralsWriter(CCZ4Solver_ADERDG&  solver)
	: IntegralsWriter() { plotForADERSolver = true; }

CCZ4::IntegralsWriter::IntegralsWriter(CCZ4Solver_FV&  solver)
	: IntegralsWriter() { plotForADERSolver = false; }

CCZ4::IntegralsWriter::IntegralsWriter() :
	fields("output/fields-", ".asc", exahype::plotters::ascii::parallel::postprocess),
	errors("output/error-", ".asc", exahype::plotters::ascii::parallel::postprocess)
{
	
	fields.add( 0, "00-gxx");
	fields.add( 1, "01-gxy");
	fields.add( 2, "02-gxz");
	fields.add( 3, "03-gyy");
	fields.add( 4, "04-gyz");
	fields.add( 5, "05-gzz");
	fields.add( 6, "06-kxx");
	fields.add( 7, "07-kxy");
	fields.add( 8, "08-kxz");
	fields.add( 9, "09-kyy");
	fields.add(10, "10-kyz");
	fields.add(11, "11-kzz");
	fields.add(12, "12-Z1");
	fields.add(13, "13-Z2");
	fields.add(14, "14-Z3");
	fields.add(15, "15-Theta");
	fields.add(16, "16-lapse");
	fields.add(17, "17-shift1");
	fields.add(18, "18-shift2");
	fields.add(19, "19-shift3");
	fields.add(20, "20-b1");
	fields.add(21, "21-b2");
	fields.add(22, "22-b3");
	fields.add(23, "23-A1");
	fields.add(24, "24-A2");
	fields.add(25, "25-A3");
	fields.add(26, "26-B11");
	fields.add(27, "27-B21");
	fields.add(28, "28-B31");
	fields.add(29, "29-B12");
	fields.add(30, "30-B22");
	fields.add(31, "31-B32");
	fields.add(32, "32-B13");
	fields.add(33, "33-B23");
	fields.add(34, "34-B33");
	fields.add(35, "35-D111");
	fields.add(36, "36-D112");
	fields.add(37, "37-D113");
	fields.add(38, "38-D122");
	fields.add(39, "39-D123");
	fields.add(40, "40-D133");
	fields.add(41, "41-D211");
	fields.add(42, "42-D212");
	fields.add(43, "43-D213");
	fields.add(44, "44-D222");
	fields.add(45, "45-D223");
	fields.add(46, "46-D233");
	fields.add(47, "47-D311");
	fields.add(48, "48-D312");
	fields.add(49, "49-D313");
	fields.add(50, "50-D322");
	fields.add(51, "51-D323");
	fields.add(52, "52-D333");
	fields.add(53, "53-traceK");
	fields.add(54, "54-phi");
	fields.add(55, "55-P1");
	fields.add(56, "56-P2");
	fields.add(57, "57-P3");
	fields.add(58, "58-K0");

	errors.add( 0, "gxx");
	errors.add( 1, "gxy");
	errors.add( 2, "gxz");
	errors.add( 3, "gyy");
	errors.add( 4, "gyz");
	errors.add( 5, "gzz");
	errors.add( 6, "kxx");
	errors.add( 7, "kxy");
	errors.add( 8, "kxz");
	errors.add( 9, "kyy");
	errors.add(10, "kyz");
	errors.add(11, "kzz");
	
	errors.add(16, "lapse");
	errors.add(54, "phi");
	errors.add(53, "traceK");
}


CCZ4::IntegralsWriter::~IntegralsWriter() {
  // @todo Please insert your code here
}


void CCZ4::IntegralsWriter::startPlotting(double time) {
  	fields.startRow(time);
	errors.startRow(time);
}


void CCZ4::IntegralsWriter::finishPlotting() {
  	fields.finishRow();
	errors.finishRow();
}


void CCZ4::IntegralsWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  	// make sure this plotter has no output associated
	assertion( outputQuantities == nullptr );

	const int nVar  = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	
	double dV;
	if(plotForADERSolver) {
		const int order = CCZ4::AbstractCCZ4Solver_ADERDG::Order;
		dV = kernels::ADERDGVolume(order, sizeOfPatch, pos);
	} else {
		const int patchSize = CCZ4::AbstractCCZ4Solver_FV::PatchSize;
		dV = tarch::la::volume(sizeOfPatch)/patchSize; // correct is probably (patchSize+1)
	}
	
	double V[nVar], ExactCons[nVar], ExactPrims[nVar];
	CCZ4Fortran::Cons2Prim(V, Q);

	// reduce the primitive quantities
	fields.addValue(V, dV);

	InitialData(x.data(),timeStamp,ExactCons);
	CCZ4Fortran::Cons2Prim(ExactPrims, ExactCons);
	
	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = std::fabs(ExactPrims[i] - V[i]);
		/*printf("ExactPrims[i] - V[i] = %f - %f = %.20e\n",
			ExactPrims[i],V[i], ExactPrims[i] - V[i]
		);*/
	}
	
	errors.addValue(localError, dV);
}


