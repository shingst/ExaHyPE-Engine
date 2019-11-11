#include "ADMIntegralsWriter.h"

FOCCZ4::ADMIntegralsWriter::ADMIntegralsWriter() : exahype::plotters::ADERDG2UserDefined::ADERDG2UserDefined(),
	simple("output/adm-integrals")
{
	simple.addSphere(5);
}

void FOCCZ4::ADMIntegralsWriter::startPlotting( double time) {
	simple.startPlotting(time);
}


void FOCCZ4::ADMIntegralsWriter::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
    double timeStamp) {

	simple.plotPatch(offsetOfPatch.data(), sizeOfPatch.data(), u, timeStamp);
}

void FOCCZ4::ADMIntegralsWriter::finishPlotting() {
	simple.finishPlotting();
}
