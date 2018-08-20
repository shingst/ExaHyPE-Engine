package eu.exahype;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;

import eu.exahype.kernel.ADERDGKernel;

public class CodeGeneratorHelper {
  
  //configuration parameters
  //------------------------
  private static String OPT_KERNELS_PATH_PREFIX = "kernels";       //Desired relative path to the generated code, starts from application root
  private static String CODEGENERATOR_PATH      = "CodeGenerator/codegenerator"; //Relative path to the CodeGenerator, starts from exahype root (ExaHyPE-Engine)
  private static String defineNamespace(String projectName, String solverName) {return projectName+"::"+solverName+"_kernels::aderdg";}  //build the generated code's namespace
  
  //options flags
  private static String useFluxOptionFlag            = "--useFlux";
  private static String useFluxVectOptionFlag        = "--useFluxVect";
  private static String useNCPOptionFlag             = "--useNCP";
  private static String useSourceOptionFlag          = "--useSource";
  private static String useFusedSourceOptionFlag     = "--useFusedSource";
  private static String useFusedSourceVectOptionFlag = "--useFusedSourceVect";
  private static String useMaterialParamOptionFlag   = "--useMaterialParam";
  private static String usePointSourcesOptionFlag    = "--usePointSources";
  private static String useLimiterOptionFlag         = "--useLimiter";
  private static String useGaussLobattoOptionFlag    = "--useGaussLobatto";
  private static String ghostLayerWidthOptionFlag    = "--ghostLayerWidth";
  private static String useCERKGuessOptionFlag       = "--useCERKGuess";
  
  
  //Internal states
  //---------------
  private Map<String,String> _optKernelsPaths;      //stores the paths to the generated code (used for imports in the KernelRegistration and in the Makefile)
  private Map<String,String> _optKernelsNamespaces; //stores the namespace used. The specific namespace depend on the solvername (assume projectname is constant)
  private static String _pathToApplication = null;  //static to not initialize the whole CodeGeneratorHelper when not required
  
  
  //Singleton pattern (to be able to access the instance everywhere)
  //-----------------
  private static volatile CodeGeneratorHelper instance = null;

  private CodeGeneratorHelper() {
    _optKernelsPaths = new HashMap<String,String>();
    _optKernelsNamespaces = new HashMap<String,String>();
  }

  public static CodeGeneratorHelper getInstance() {
      if (instance == null) {
          synchronized(CodeGeneratorHelper.class) {
              if (instance == null) {
                  instance = new CodeGeneratorHelper();
              }
          }
      }
      return instance;
  }
  
  // shortcut method to build a unique key out of projectname + solver name
  private static String getKey(String projectName, String solverName) {
    return projectName + "::" + solverName;
  }
  
  
  //Setter
  //------
  public static void setPaths(DirectoryAndPathChecker directoryAndPathChecker) {
    _pathToApplication = directoryAndPathChecker.outputDirectory.getPath();
  }
  
  
  //Getter
  //------
  
  public String getIncludePath(String projectName, String solverName) {
    return _optKernelsPaths.get(getKey(projectName,solverName));
  }
  
  public Collection<String> getIncludePaths() {
    return _optKernelsPaths.values();
  }
  
  public String getNamespace(String projectName, String solverName) {
    return _optKernelsNamespaces.get(getKey(projectName,solverName));
  }
  
  public Collection<String> getNamespaces() {
    return _optKernelsNamespaces.values();
  }
  
  
  //Generate code
  //-------------
  public String invokeCodeGenerator(String projectName, String solverName, int numberOfUnknowns, int numberOfParameters, int order, int dimensions, String microarchitecture, ADERDGKernel kernel)
      throws IOException {
    
    //check and defines paths       
    if(_pathToApplication == null) {
      System.err.println("ERROR: Path to the application for the CodeGenerator not found");
      throw new IOException();
    }
              
    java.io.File pathToCodeGenerator_File =
        new java.io.File(CodeGeneratorHelper.CODEGENERATOR_PATH);
    if (!pathToCodeGenerator_File.exists()) {
      System.err.println("ERROR: CodeGenerator not found. Can't generate optimised kernels. Path: " + pathToCodeGenerator_File.getCanonicalPath());
      throw new IOException();
    }
    String pathToCodeGenerator = pathToCodeGenerator_File.getCanonicalPath();
    String optKernelPath = (new java.io.File(OPT_KERNELS_PATH_PREFIX,solverName)).getPath();
    
    //define the CodeGenerator arguments
    String namespace = defineNamespace(projectName, solverName);    
    String numericsParameter = kernel.isLinear() ? "linear" : "nonlinear";
    String options =  (kernel.useFlux() ? (kernel.useFluxVect() ? useFluxVectOptionFlag : useFluxOptionFlag)+" " : "")
                    + (kernel.useSource() ? (kernel.useFusedSourceVect() ? useFusedSourceVectOptionFlag : (kernel.useFusedSource() ? useFusedSourceOptionFlag : useSourceOptionFlag))+" " : "") 
                    + (kernel.useNCP() ?  useNCPOptionFlag+" " : "") 
                    + (kernel.usePointSources() ?  usePointSourcesOptionFlag+" "+kernel.getNumberOfPointSources()+" " : "") 
                    + (kernel.useCERKGuess() ? useCERKGuessOptionFlag+" " : "") 
                    + (kernel.useMaterialParameterMatrix() ? useMaterialParamOptionFlag+" " : "")
                    + (kernel.useGaussLobatto() ? useGaussLobattoOptionFlag+" " : "")
                    + (kernel.useLimiter() ?  useLimiterOptionFlag+" "+kernel.getNumberOfObservables()+" " : "")
                    ;

    // set up the command to execute the code generator
    String args =   " " + _pathToApplication 
                  + " " + optKernelPath 
                  + " " + namespace
                  + " " + projectName + "::" + solverName
                  + " " + numberOfUnknowns 
                  + " " + numberOfParameters 
                  + " " + order 
                  + " " + dimensions 
                  + " " + numericsParameter 
                  + " " + microarchitecture 
                  + " " + options; 
    
    String bashCommand = "env python3 "  + pathToCodeGenerator +  args;

    // execute the command line program
    Runtime runtime = Runtime.getRuntime();
    System.out.println("CodeGenerator command line: "+bashCommand);
    Process codeGenerator = runtime.exec(bashCommand);

    // capture any output that is produced by the code generator and print it line-by-line
    java.io.InputStream stdout = codeGenerator.getInputStream();
    java.io.BufferedReader stdoutReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stdout));
    String line = "";
    while ((line = stdoutReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }
    java.io.InputStream stderr = codeGenerator.getErrorStream();
    java.io.BufferedReader stderrReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stderr));
    while ((line = stderrReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }

    // in order to stop further toolkit execution if the code generator fails,
    // explicitly wait for the process
    try {
        int exitValue = codeGenerator.waitFor();
        if(exitValue != 0) {
            System.err.println("ERROR: CodeGenerator failed with exit value " + exitValue);
            throw new IOException();
        }
    } catch(InterruptedException e) {
        System.err.println("This is very bad. I don't know what's going on.");
        throw new IOException();
    }
    
    _optKernelsPaths.put(getKey(projectName,solverName),optKernelPath);
    _optKernelsNamespaces.put(getKey(projectName,solverName), namespace);
    
    return optKernelPath;
    
  } // invokeCodeGenerator
  
}
