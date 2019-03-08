Output at 2019-03-02 is:

sven@nils:~/numrel/exahype/Engine-ExaHyPE/ApplicationExamples/Experiments/HybridInterfaceDemonstration$ ./ExaHyPE-HybridInterfaceDemonstration --built-in-specfile | grep '>>'
>> >> called HybridInterfaceDemonstration::StandardWriter::StandardWriter
>> >> called HybridInterfaceDemonstration::DummySolver_ADERDG::init
>> >> called HybridInterfaceDemonstration::DummySolver_FV::init
>> >> called HybridInterfaceDemonstration::DummySolver_ADERDG::adjustPointSolution
>> Setting DG ID since t=0
>> >> called HybridInterfaceDemonstration::DummySolver_ADERDG::isPhysicallyAdmissible
>> >> called HybridInterfaceDemonstration::DummySolver_FV::adjustSolution
>> Setting FV ID since t=0
>> >> called HybridInterfaceDemonstration::DummySolver_ADERDG::eigenvalues
>> >> called HybridInterfaceDemonstration::DummySolver_ADERDG::beginTimeStep
>> >> called HybridInterfaceDemonstration::StandardWriter::startPlotting
>> >> called HybridInterfaceDemonstration::StandardWriter::mapQuantities
>> >> called HybridInterfaceDemonstration::StandardWriter::finishPlotting
>> >> called HybridInterfaceDemonstration::DummySolver_FV::boundaryValues
>> >> called HybridInterfaceDemonstration::DummySolver_FV::algebraicSource
>> >> called HybridInterfaceDemonstration::DummySolver_FV::nonConservativeProduct
>> >> called HybridInterfaceDemonstration::DummySolver_FV::flux
>> >> called HybridInterfaceDemonstration::DummySolver_FV::eigenvalues
>> >> called HybridInterfaceDemonstration::DummySolver_ADERDG::endTimeStep

