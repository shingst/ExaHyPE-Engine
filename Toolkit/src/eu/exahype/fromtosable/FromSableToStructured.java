package eu.exahype.fromtosable;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.*;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.*;

import eu.exahype.variables.Variables;


/**
 * A SableCC visitor to transform a specfile into a structured
 * container.
 *
 * Which subsequently can be exported to a format such as JSON.
 *
 * This can be helpful to process specfiles without rough
 * regexpes or similiar.
 *
 **/
public class FromSableToStructured extends DepthFirstAdapter {
  /// This maps an [hierarchical path key] -> [value]
  // private java.util.Map<String, String> _data;
  Structured _data;
  
  public FromSableToStructured() {
    //_data = new java.util.HashMap<String,String>();
    _data = new StructuredDumper();
  }
  
  public void dump() {
  // where does entry come from
  /*   for(Entry<String,String> entry : _data.entrySet()) {
        System.out.println(entry.getKey() + "=" + entry.getValue());
     }
  */
  }
  
  public void mandatory(String target, Token token) {
    _data.put(target, token.getText());
  }
  
  public boolean optional(String target, Token token) {
    if(token != null)
      _data.put(target, token.getText());
    return (token != null);
  }

  @Override
  public void inAPaths(APaths node) {
    // see below for an implementation which traverses nicer.
    System.out.println("APath: " + node.toString());
  }

  @Override
  public void inAProject(AProject node) {
    // for debugging:
    System.out.println("AProject: " + node.toString());
  
    mandatory("project/name", node.getName());
    
    APaths paths = (APaths) node.getPaths();
    mandatory("project/paths/peano-path", paths.getPeanoKernelPath());
    mandatory("project/paths/exahype-path", paths.getExahypePath());
    mandatory("project/paths/output-directory",  paths.getOutputDirectory());
    
    AArchitecture architecture = (AArchitecture) node.getArchitecture();
    mandatory("project/get-architecture", architecture.getIdentifier());
    
    ALogfile logfile = (ALogfile) node.getLogfile();
    optional("project/logfile", logfile.getFilename());
    
    AComputationalDomain computationalDomain = (AComputationalDomain) node.getComputationalDomain();
    // *probably* todo: As integer. Better not here but detect
    // somewhere centrally.
    int dimensions = Integer.parseInt(computationalDomain.getDimension().getText());
    _data.put("project/computational_domain/dimensions", dimensions);
    
    mandatory("project/computational_domain/width/x", computationalDomain.getWidthX());
    mandatory("project/computational_domain/width/y", computationalDomain.getWidthY());
    if(dimensions>2)
    mandatory("project/computational_domain/width/z", computationalDomain.getWidthZ());

    mandatory("project/computational_domain/offset/x", computationalDomain.getOffsetX());
    mandatory("project/computational_domain/offset/y", computationalDomain.getOffsetY());
    if(dimensions>2)
    mandatory("project/computational_domain/offset/z", computationalDomain.getOffsetZ());

    optional("project/computational_domain/endTime", computationalDomain.getEndTime());
    optional("project/computational_domain/timeSteps", computationalDomain.getTimeSteps());
    
    ASharedMemory sharedMemory = (ASharedMemory) node.getSharedMemory();
    if(sharedMemory != null) {
        mandatory("project/shared_memory/identifier", sharedMemory.getIdentifier());
        mandatory("project/shared_memory/cores", sharedMemory.getCores());
        mandatory("project/shared_memory/properties_file", sharedMemory.getPropertiesFile());
    }
  }


  /// Child of AProject
  @Override
  public void inAProfiling(AProfiling node) {
    /*
        @SuppressWarnings("hiding") TIdentifier _profiler_,
        @SuppressWarnings("hiding") List<PItem> _metrics_,
        @SuppressWarnings("hiding") TFilename _likwidInc_,
        @SuppressWarnings("hiding") TFilename _likwidLib_,
        @SuppressWarnings("hiding") TFilename _ipcmInc_,
        @SuppressWarnings("hiding") TFilename _ipcmLib_,
        @SuppressWarnings("hiding") TTokenOnOff _deepProfiling_)
    */
  };
  
  /// Child of AProject
  @Override
  public void inADistributedMemory(ADistributedMemory node) {
        mandatory("project/distributed_memory/identifier", node.getIdentifier());
        mandatory("project/distributed_memory/configure", node.getConfigure());
        mandatory("project/distributed_memory/buffer_size", node.getBuffersize());
        mandatory("project/distributed_memory/timeout", node.getTimeout());
  }
  
  
  /// Child of AProject
  @Override
  public void inAOptimisation(AOptimisation node) {
        String prefix = "project/distributed_memory";
        optional(prefix+"fuseAlgorithmSteps", node.getFuseAlgorithmSteps());
        optional(prefix+"fuseAlgorithmStepsFactor", node.getFuseAlgorithmStepsFactor());
        optional(prefix+"spawnPredictor", node.getSpawnPredictor());
        optional(prefix+"batchTimesteps", node.getBatchTimesteps());
        optional(prefix+"skipReduction", node.getSkipReduction());
        optional(prefix+"disableAmr", node.getDisableAmr());
        // usw., for
        //  TFloatNumber _doubleCompression_,
        // TTokenOnOff _spawnDoubleCompression_)
  }
  
  public void registerVariables(Variables quantities) {
  }

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    String solverName    = node.getName().getText();
    _data.put("domain/solver/name", node.getName().getText());
    _data.put("domain/solver/language", node.getLanguage().getText());
    int     order        = Integer.parseInt(node.getOrder().getText());
    // interestingly, here the constants start:
    boolean hasConstants = node.getConstants()!=null;
    // make use of the structured variables class
    Variables quantities  = new Variables(solverName, node);

    // can feed this into the data.
    Map<String,Integer> variables = quantities.getVariablesMap();
    Map<String,Integer> parameters = quantities.getParametersMap();    
  }

  @Override
  public void inAFiniteVolumesSolver(AFiniteVolumesSolver node) {
    String solverName    = node.getName().getText();
    String  language     = node.getLanguage().getText();
    int     patchSize    = Integer.parseInt(node.getPatchSize().getText());
    boolean hasConstants = node.getConstants()!=null;
    Variables variables  = new Variables(solverName, node);
    boolean isFortran    = language.equals("Fortran");
  }
  
  @Override
  public void inALimitingAderdgSolver(ALimitingAderdgSolver node) {
    String solverName    = node.getName().getText();
    String  language     = node.getLanguage().getText();
    int     order        = Integer.parseInt(node.getOrder().getText());
    boolean hasConstants = node.getConstants()!=null;
    String  limiterLanguage  = node.getLanguageLimiter().getText();
  }
}
