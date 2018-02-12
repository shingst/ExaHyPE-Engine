package eu.exahype.solvers;

import java.util.Set;
import java.util.List;
import java.util.Arrays;

import eu.exahype.kernel.ADERDGKernel;
import eu.exahype.kernel.FiniteVolumesKernel;

/**
 * Creates a solver
 * 
 * The idea is that you pass in a kernel object (either FiniteVolumesKernel or 
 * ADERDGKernel) which returns an enum describing the solver type. This 
 * information together with additional properties then is used to instantiate
 * the correct solver class.
 * 
 * <h2> Add support for a new solver type </h2>
 *
 * Switch to either the FiniteVolumesKernel or ADERDGKernel. Add your new solver 
 * type to the enum called KernelType. Next, add a new branch to the routine
 * getKernelType() which returns your brand new solver. Finally, extend the 
 * switch statements in this class to instantiate this new solver.
 *
 */
public class SolverFactory {
  private String _projectName;
  private int _dimensions;
  private boolean _enableProfiler;
  private boolean _enableDeepProfiler;
  private String _microarchitecture;
  private String _pathToLibxsmm;

  
  public SolverFactory(
      String projectName,
      int dimensions,
      boolean enableProfiler,
      boolean enableDeepProfiler,
      String microarchitecture) {
    _projectName = projectName;  
    _dimensions = dimensions;
    _enableProfiler = enableProfiler;
    _enableDeepProfiler = enableDeepProfiler;
    _microarchitecture = microarchitecture;
  }
  
  /**
   * Generates the writer for an ADER-DG solver
   * 
   * Consult class documentation for further details.
   */
  public Solver createADERDGSolver(String solverName, ADERDGKernel kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int order,boolean hasConstants) {
    try { //some solver initialisation can throw IllegalArgumentException if the options are wrong or IOException
      switch (kernel.getKernelType()) {
        case GenericADERDG: 
          return new eu.exahype.solvers.GenericADERDG(_projectName, solverName, _dimensions,
            numberOfVariables, numberOfParameters, namingSchemeNames, order, _enableProfiler, hasConstants, isFortran, kernel );
        case OptimisedADERDG:
          return new eu.exahype.solvers.OptimisedADERDG(_projectName, solverName, _dimensions,
            numberOfVariables, numberOfParameters, namingSchemeNames, order, _microarchitecture,
            _enableProfiler, _enableDeepProfiler, hasConstants, kernel);
      }
      System.err.println("ERROR: solver configuration is not supported: "+kernel.toString() );
      return null;
    } catch(Exception e) {
      System.err.println("ERROR: can't create the solver. Error: "+e );
      return null;
    }
  }
  
  /**
   * Generates the writer for an Finite Volumes solver
   * 
   * Consult class documentation for further details.
   */
  public Solver createFiniteVolumesSolver(String solverName, FiniteVolumesKernel kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int patchSize,boolean hasConstants) {
    try { //some solver initialisation can throw IllegalArgumentException if the options are wrong or IOException
      switch (kernel.getKernelType()) {
        case GenericMUSCLHancock: 
          if (isFortran) {
        	// @todo Does not exist yet
            //return new eu.exahype.solvers.GenericFiniteVolumesMUSCLHancockInFortran(_projectName, solverName,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
          }
          else {
            return new eu.exahype.solvers.GenericFiniteVolumesMUSCLHancockInC(_projectName, solverName,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants, kernel);
          }
          break;
        case GenericGodunov: 
            if (isFortran) {
              // @todo Does not exist yet
              //return new eu.exahype.solvers.GenericFiniteVolumesGodunovInFortran(_projectName, solverName,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
            }
            else {
              return new eu.exahype.solvers.GenericFiniteVolumesGodunovInC(_projectName, solverName,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants, kernel);
            }
            break;
        case UserDefined: 
            if (isFortran) {
              return new eu.exahype.solvers.UserDefinedFiniteVolumesinFortran(_projectName, solverName,_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
            }
            else {
              return new eu.exahype.solvers.UserDefinedFiniteVolumesinC(_projectName, solverName,_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
            }
        }
      System.err.println("ERROR: solver configuration is not supported: "+kernel.toString() );
      return null;
    } catch(Exception e) {
      System.err.println("ERROR: can't create the solver. Error: "+e );
      return null;
    }
  }
  
  public Solver createLimiterSolver(String solverName, ADERDGKernel aderdgKernel, FiniteVolumesKernel FVKernel, Solver ADERDGsolver, Solver FVsolver) { //take the kernel to know if generic or optimised
    try {
      return new eu.exahype.solvers.Limiter(_projectName, solverName, ADERDGsolver, FVsolver); //TODO JMG optimized vs generic
    } catch(Exception e) {
      System.err.println("ERROR: can't create the solver. Error: "+e );
      return null;
    }
  }
}
