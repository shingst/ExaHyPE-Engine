package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class UserDefinedFiniteVolumesinC implements Solver {
  public static final String Identifier = "user::defined";

  private String _projectName;
  private String _solverName;
  private int _numberOfVariables;
  private int _numberOfParameters;
  private int _patchSize;
  private boolean _hasConstants;

  public UserDefinedFiniteVolumesinC(String projectName, String solverName, int dimensions, int numberOfVariables, int numberOfParameters, int patchSize, boolean enableProfiler, boolean hasConstants) {
    _projectName        = projectName;
    _solverName         = solverName;
    _numberOfVariables  = numberOfVariables;
    _numberOfParameters = numberOfParameters;
    _patchSize          = patchSize;
    _hasConstants       = hasConstants;
  }

  @Override
  public String getSolverName() {
    return _solverName;
  }
  
  @Override
  public void writeAbstractHeader(BufferedWriter writer) throws IOException {
    // TODO Auto-generated method stub
    
  }
  
  @Override
  public void writeAbstractImplementation(BufferedWriter writer) throws IOException {
    // TODO Auto-generated method stub
    
  }
  
  public void writeHeader(java.io.BufferedWriter writer)
      throws java.io.IOException {
     Helpers.writeMinimalFiniteVolumesSolverHeader(_solverName, writer, _projectName, _hasConstants);

    writer.write("};\n\n\n");
  }

  // @todo THIS MUST BE REVSISED!!!
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException {
    writer.write("#include \"" + _solverName + ".h\"\n");
    writer.write("\n\n\n");
    
    //    if (_hasConstants) {
    //      writer.write(_projectName + "::" + _solverName + "::" + _solverName + "(int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler, exahype::Parser::ParserView parser):\n");
    //    }
    //    else {
    //      
    //    }
    writer.write(_projectName + "::" + _solverName + "::" + _solverName + "(int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):\n");
    writer.write("  exahype::solvers::FiniteVolumesSolver("
            + "\""+_solverName+"\", "+_numberOfVariables+"/* numberOfVariables */, "+_numberOfParameters+" /* numberOfParameters */, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {\n");
    writer.write("  // @todo Please implement/augment if required\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + _projectName + "::" + _solverName
        + "::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("double " + _projectName + "::" + _solverName
        + "::stableTimeStepSize(double* luh, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return 1.0;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
        + "::solutionAdjustment(double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("bool " + _projectName + "::" + _solverName
            + "::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return false;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("exahype::solvers::Solver::RefinementControl " + _projectName + "::" + _solverName
            + "::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return exahype::solvers::Solver::RefinementControl::Keep;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
            + "::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
            + "::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
            + "::volumeUnknownsProlongation(  double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    writer.write("  //@todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + _solverName
            + "::volumeUnknownsRestriction(  double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }
  public void writeUserPDE(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a PDF.f90.\n");
  }
  public void writeTypesDef(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a typesDef.f90.\n");
  }
  
  @Override
  public boolean supportsVariables() {
    return false;
  }
}
