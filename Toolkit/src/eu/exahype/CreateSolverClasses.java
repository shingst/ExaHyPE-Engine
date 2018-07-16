package eu.exahype;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.Set;
import java.util.Collections;
import java.util.Collection;
import java.util.Arrays;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.kernel.ADERDGKernel;
import eu.exahype.kernel.FiniteVolumesKernel;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AComputationalDomain;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;
import eu.exahype.node.PSolver;
import eu.exahype.node.AArchitecture;
import eu.exahype.node.PArchitecture;
import eu.exahype.solvers.Solver;
import eu.exahype.solvers.SolverFactory;
import eu.exahype.variables.Variables;
import eu.exahype.io.FileSearch;
import eu.exahype.io.IOUtils;

public class CreateSolverClasses extends DepthFirstAdapter {
  
  private static final Collection<String> SUPPORTED_MICROARCHITECTURES = Collections.unmodifiableCollection(Arrays.asList("wsm", "snb", "hsw", "knc", "knl", "noarch"));
  
  public Boolean valid = true;
  
  private String _outputPath;
  private String _projectName;
  private String _microarchitecture;
  private Collection<String> _definedSolvers;
  private int _dimensions;
  private boolean _enableProfiler;

  public CreateSolverClasses(DirectoryAndPathChecker directoryAndPathChecker) {
    _outputPath = directoryAndPathChecker.outputDirectory.getAbsolutePath();
    _enableProfiler = false;
  }

  @Override
  public void inAProject(AProject node) {
    _projectName     = node.getName().getText();
    _definedSolvers  = new java.util.HashSet<String>();

    if (node.getSolver().size() == 0) {
      System.out.println("there are no solvers in the specification file ... nothing to be done");
    }

    if (node.getArchitecture()!=null) {
      _microarchitecture = ( (eu.exahype.node.AArchitecture)( node.getArchitecture() )).getIdentifier().getText().toLowerCase();
    }
    else {
      _microarchitecture = "noarch";
    }
    if (!SUPPORTED_MICROARCHITECTURES.contains(_microarchitecture)) {
      System.out.println("Unknown architecture specified ... fallback solution \"noarch\" taken");
      _microarchitecture = "noarch";
    }
  }

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    _dimensions = Integer.parseInt( node.getDimension().getText() );
    if (_dimensions!=2 && _dimensions!=3) {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
    }
  }


  @Override
  public void inAProfiling(AProfiling node) {
    _enableProfiler = !node.getProfiler().getText().equals("NoOpProfiler");
  };

  // @todo This function should be a member of Solver.java.
  private boolean validate(Variables variables, int order, String kernel, String language,
      String solverName,eu.exahype.solvers.Solver solver) {
    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiple definitions." );
      return false;
    }
    
    if (!language.equals("C") && !language.equals("Fortran")) {
      System.err.println("ERROR: unknown language for solver " + solverName
          + ". Supported languages are C and Fortran");
      return false;
    }
    
    // @todo Does not work with Finite Volumes
/*    if (order < 1 || order > 9) {
      System.err.println("ERROR: Only polynomial degrees of 1..9 are supported.");
      return false;
    }
*/
    
    if (solver == null) {
      System.err.println("ERROR: creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + language + " is not supported");
      return false;
    }
    return true;
  }

  
  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    String  solverName   = node.getName().getText();
    String  language     = node.getLanguage().getText();
    int     order        = Integer.parseInt(node.getOrder().getText());
    boolean hasConstants = node.getConstants()!=null;
    Variables variables  = new Variables(solverName, node);
    boolean isFortran    = language.equals("Fortran");

    try {
      ADERDGKernel kernel = new ADERDGKernel(node);
    
      SolverFactory solverFactory = new SolverFactory(_projectName, _dimensions, _enableProfiler, _microarchitecture);
      eu.exahype.solvers.Solver solver = solverFactory.createADERDGSolver(
          solverName, kernel, isFortran, variables.getNumberOfVariables(), variables.getNumberOfParameters(),variables.getNamingSchemeNames(), order, hasConstants);
      valid = validate(variables,order,kernel.toString(),language,solverName,solver);
      
      if (valid) {
        _definedSolvers.add(solverName);

        // write the files
      
        tryWriteSolverHeader(solver);
        tryWriteSolverUserImplementation(solver);
        
        tryWriteAbstractSolverHeader(solver);
        tryWriteAbstractSolverImplementation(solver);

        if (solver.supportsVariables()) {
          tryWriteVariablesHeader(variables);
        }
      }
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }

  @Override
  public void inAFiniteVolumesSolver(AFiniteVolumesSolver node) {
    String solverName    = node.getName().getText();
    String  language     = node.getLanguage().getText();
    int     patchSize    = Integer.parseInt(node.getPatchSize().getText());
    boolean hasConstants = node.getConstants()!=null;
    Variables variables  = new Variables(solverName, node);
    boolean isFortran    = language.equals("Fortran");
    
    SolverFactory solverFactory = new SolverFactory(_projectName, _dimensions, _enableProfiler, _microarchitecture);
    FiniteVolumesKernel kernel  = new FiniteVolumesKernel(node);
    eu.exahype.solvers.Solver solver = solverFactory.createFiniteVolumesSolver(
        solverName, kernel, isFortran, variables.getNumberOfVariables(), variables.getNumberOfParameters(),variables.getNamingSchemeNames(), patchSize, hasConstants);

    valid = validate(variables,patchSize,kernel.toString(),language,solverName,solver);

    if (valid) {
      _definedSolvers.add(solverName);

      // write the files
      try {
        tryWriteSolverHeader(solver);
        tryWriteSolverUserImplementation(solver);
        
        tryWriteAbstractSolverHeader(solver);
        tryWriteAbstractSolverImplementation(solver);

        if (solver.supportsVariables()) {
          tryWriteVariablesHeader(variables);
        }
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        exc.printStackTrace();
        valid = false;
      }
    }
  }
  
  @Override
  public void inALimitingAderdgSolver(ALimitingAderdgSolver node) {
    String solverName        = node.getName().getText();
    String  language         = node.getLanguage().getText();
    int     order            = Integer.parseInt(node.getOrder().getText());
    int     DmpObservables   = Integer.parseInt(node.getDmpObservables().getText());
    boolean hasConstants     = node.getConstants()!=null;
    
    boolean isFortran        = language.equals("Fortran");
    
    int     patchSize        = 2*order+1;
    String  limiterLanguage  = node.getLanguageLimiter().getText();
    
    String solverNameADERDG  = solverName+"_ADERDG";
    String solverNameFV      = solverName+"_FV";
    String solverNameLimiter = solverName+"_Limiter";
    
    Variables variablesSolver  = new Variables(solverNameADERDG, node);
    Variables variablesLimiter = new Variables(solverNameFV, node);
    try {
      ADERDGKernel         aderdgKernel  = new ADERDGKernel(node);
      FiniteVolumesKernel  FVKernel      = new FiniteVolumesKernel(node);
      aderdgKernel.setNumberOfObservables(DmpObservables);
      
      SolverFactory solverFactory = new SolverFactory(_projectName, _dimensions, _enableProfiler, _microarchitecture);
      Solver solverAderdg  = solverFactory.createADERDGSolver(
          solverNameADERDG, aderdgKernel,isFortran,variablesSolver.getNumberOfVariables(),variablesSolver.getNumberOfParameters(),variablesSolver.getNamingSchemeNames(),order,hasConstants);
      Solver solverFV = solverFactory.createFiniteVolumesSolver(
          solverNameFV, FVKernel,isFortran,variablesLimiter.getNumberOfVariables(),variablesLimiter.getNumberOfParameters(),variablesLimiter.getNamingSchemeNames(),patchSize,hasConstants);
      Solver solverLimiter = solverFactory.createLimiterSolver(solverNameLimiter, aderdgKernel, FVKernel, solverAderdg, solverFV);

      valid  = validate(variablesSolver,order,aderdgKernel.toString(),language,solverName,solverAderdg);
      valid &= validate(variablesLimiter,patchSize,FVKernel.toString(),limiterLanguage,solverName,solverFV);
      
      if (valid) {        
        _definedSolvers.add(solverName);

        // write the files
        tryWriteSolverHeader(solverAderdg);
        tryWriteSolverHeader(solverFV);

        tryWriteSolverUserImplementation(solverAderdg);
        tryWriteSolverUserImplementation(solverFV);
        
        tryWriteAbstractSolverHeader(solverAderdg);
        tryWriteAbstractSolverHeader(solverFV);
        tryWriteAbstractSolverHeader(solverLimiter);
        
        tryWriteAbstractSolverImplementation(solverAderdg);
        tryWriteAbstractSolverImplementation(solverFV);
        tryWriteAbstractSolverImplementation(solverLimiter);

        if (solverAderdg.supportsVariables()) {
          tryWriteVariablesHeader(variablesSolver);
        }
        if (solverFV.supportsVariables()) {
          tryWriteVariablesHeader(variablesLimiter);
        }
      }  
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }
  
  private void tryWriteSolverHeader(Solver solver) throws IOException,IllegalArgumentException {
    java.io.File solverHeaderFile = FileSearch.relocatableFile(
        _outputPath + "/" + solver.getSolverName() + ".h");
    
    if (solverHeaderFile.exists()) {
      System.out.println("create header of solver " + solver.getSolverName() + " ... header "
          + solverHeaderFile.getAbsoluteFile()
          + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)");
    } else {
      BufferedWriter headerWriter =
          new BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      solver.writeHeader(headerWriter);
      System.out.println("create header of solver " + solver.getSolverName() + " ... ok");
      headerWriter.close();
    }
  }
  
  private void tryWriteSolverUserImplementation(Solver solver) throws IOException,IllegalArgumentException {
    java.io.File solverUserImplementationFile = FileSearch.relocatableFile(
        _outputPath + "/" + solver.getSolverName() + ".cpp");
    
    if (solverUserImplementationFile.exists()) {
      System.out.println("user's implementation file of solver " + solver.getSolverName()
          + " ... does exist already. Is not overwritten");
    } else {
      BufferedWriter userImplementationWriter =
          new BufferedWriter(new java.io.FileWriter(solverUserImplementationFile));
      solver.writeUserImplementation(userImplementationWriter);
      System.out.println(
          "create user implementation template of solver " + solver.getSolverName() + " ... please complete");
      userImplementationWriter.close();
    }
  }

  private void tryWriteAbstractSolverHeader(Solver solver) throws IOException,IllegalArgumentException {
    java.io.File abstractSolverHeaderFile = FileSearch.relocatableFile(
        _outputPath + "/Abstract" + solver.getSolverName() + ".h");
    
    if (abstractSolverHeaderFile.exists()) {
      System.out.println("implementation file for abstract solver superclass Abstract" + solver.getSolverName()
          + " ... does exist already. Is overwritten");
    }

    BufferedWriter writer =
        new BufferedWriter(new java.io.FileWriter(abstractSolverHeaderFile));
    solver.writeAbstractHeader(writer);
    System.out.println("create header file for abstract solver superclass Abstract" + solver.getSolverName() + " ... ok");
    writer.close();
  }
  
  private void tryWriteAbstractSolverImplementation(Solver solver) throws IOException,IllegalArgumentException {
    java.io.File abstractSolverImplementationFile = FileSearch.relocatableFile(
        _outputPath + "/Abstract" + solver.getSolverName() + ".cpp");
    
    if (abstractSolverImplementationFile.exists()) {
      System.out.println("implementation file for abstract solver superclass Abstract" + solver.getSolverName()
          + " ... does exist already. Is overwritten");
    }

    BufferedWriter writer =
        new BufferedWriter(new java.io.FileWriter(abstractSolverImplementationFile));
    solver.writeAbstractImplementation(writer);
    System.out.println("create implementation file for abstract solver superclass Abstract" + solver.getSolverName() + " ... ok");
    writer.close();
  }

  
  private void tryWriteVariablesHeader(Variables variables) throws IOException {
    java.io.File solverHeaderFile = FileSearch.relocatableFile(
        _outputPath + "/" + variables.getSolverName() + "_Variables.h");
    
    if (solverHeaderFile.exists()) {
      BufferedWriter headerWriter =
          new BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      variables.writeHeader(headerWriter, variables.getSolverName(), _projectName);
      System.out.println("create header of variables for solver " + variables.getSolverName() + " ... header "
          + solverHeaderFile.getAbsoluteFile()
          + " does exist already and will be overwritten");
      headerWriter.close();
    } else {
      BufferedWriter headerWriter =
          new BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      variables.writeHeader(headerWriter, variables.getSolverName(), _projectName);
      System.out.println("create header of variables for solver " + variables.getSolverName() + " ... ok");
      headerWriter.close();
    }
  }
}
