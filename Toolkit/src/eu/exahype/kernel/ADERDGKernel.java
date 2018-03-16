package eu.exahype.kernel;

import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import eu.exahype.node.PIds;
import eu.exahype.node.AIds;
import eu.exahype.node.AIdentifierId;
import eu.exahype.node.AIdentifierWithMultId;

import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.PSolver;

/**
 * Circumscribes an ADERDG Kernel
 * 
 * For a description how to add new variants consult SolverFactory.
 */
public class ADERDGKernel {
  
   
  /**
   * Configuration parameter: id of the options
   */
  public static final String GENERIC_OPTION_ID            = "generic";
  public static final String OPTIMISED_OPTION_ID          = "optimised";

  public static final String LINEAR_OPTION_ID             = "linear";
  public static final String NONLINEAR_OPTION_ID          = "nonlinear";
  public static final String USER_DEFINED_OPTION_ID       = "user";
  public static final String LEGENDRE_OPTION_ID           = "gausslegendre";
  public static final String LOBATTO_OPTION_ID            = "gausslobatto";
    
  public static final String FLUX_OPTION_ID               = "flux";
  public static final String SOURCE_OPTION_ID             = "source";
  public static final String NCP_OPTION_ID                = "ncp";
  public static final String POINTSOURCES_OPTION_ID       = "pointsources";
  public static final String MATERIALPARAMETER_OPTION_ID  = "materialparameters";

  public static final String NO_TIME_AVG_OPTION_ID        = "notimeavg";
  public static final String PATCHWISE_ADJUST_OPTION_ID   = "patchwiseadjust";
  public static final String TEMP_VARS_ON_STACK_OPTION_ID = "usestack";
  public static final String MAX_PICARD_ITER_ID           = "maxpicarditer";
  public static final String CONVERTER_OPTION_ID          = "converter"; //for debug only, not in guidebook
  public static final String FLOPS_OPTION_ID              = "flops"; //for debug only, not in guidebook
  
  private Set<String> type;
  private Map<String, Integer> terms;
  private Map<String, Integer> optimisation;
  
  private int ghostLayerWidth = -1;
  private int numberOfObservables = -1;
  
  public ADERDGKernel(PSolver solver) throws IllegalArgumentException {
    if(solver instanceof AAderdgSolver) {
      type = parseIds(((AAderdgSolver) solver).getKernelType());
      terms = parseIdsToMap(((AAderdgSolver) solver).getKernelTerms());
      optimisation =  parseIdsToMap(((AAderdgSolver) solver).getKernelOpt());
    } else if(solver instanceof ALimitingAderdgSolver) {
      type = parseIds(((ALimitingAderdgSolver) solver).getKernelType());
      terms = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelTerms());
      optimisation =  parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelOpt());
    } else {
      throw new IllegalArgumentException("No kernel definition found");
    }
    
    validate();
  }
  
  //return null on error, use only after the program should already have failed with invalid kernel
  public static ADERDGKernel noExceptionContructor(PSolver solver) {
    try {
      return new ADERDGKernel(solver);
    } catch(IllegalArgumentException e) {
      return null;
    }
  }
  
  private static Set<String> parseIds(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().map(e -> (e instanceof AIdentifierId) ? ((AIdentifierId)e).getValue().getText() : ((AIdentifierWithMultId)e).getValue().getText()).collect(Collectors.toSet());
  }
  
  private static Map<String, Integer> parseIdsToMap(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().collect(Collectors.toMap(
      (e -> (e instanceof AIdentifierId) ? ((AIdentifierId)e).getValue().getText() : ((AIdentifierWithMultId)e).getValue().getText()),
      (e -> (e instanceof AIdentifierId) ? Integer.valueOf(-1) : Integer.valueOf(((AIdentifierWithMultId)e).getMultiplicity().getText()))
    ));
  }
  
  private void validate() throws IllegalArgumentException {
    if(!type.contains(LINEAR_OPTION_ID) ^ type.contains(NONLINEAR_OPTION_ID)) {//should be only one
      throw new IllegalArgumentException("nonlinear or linear not specified or both specified in the kernel type");
    }
    if(usePointSources() && getNumberOfPointSources() < 0) {
      throw new IllegalArgumentException("point sources used but number not specified! In the specification file, use "+POINTSOURCES_OPTION_ID+":X, with X the number of point sources.");
    }
  }

  public enum KernelType {
    GenericADERDG,
    OptimisedADERDG,
    UserDefined,
    Unknown
  }

  public boolean isLinear() throws IllegalArgumentException {
    return type.contains(LINEAR_OPTION_ID);
  }

  public KernelType getKernelType() {
    if (optimisation.containsKey(OPTIMISED_OPTION_ID)) {
     return KernelType.OptimisedADERDG;
   }
      
   // default kernel - must be last   
   if(optimisation.containsKey(GENERIC_OPTION_ID)) {
      return KernelType.GenericADERDG;
  }
  
    return  KernelType.Unknown;
  }

  public boolean usesOptimisedKernels() {
    // assert: !optimisation.containsKey(GENERIC_OPTION_ID)
    return optimisation.containsKey(OPTIMISED_OPTION_ID);
  }

  public boolean useFlux() {
    return terms.containsKey(FLUX_OPTION_ID);
  }
  
  public boolean useSource() {
    return terms.containsKey(SOURCE_OPTION_ID);
  }
  
  public boolean useNCP() {
    return terms.containsKey(NCP_OPTION_ID);
  }
  
  public boolean usePointSources() {
    return terms.containsKey(POINTSOURCES_OPTION_ID);
  }
  
  public boolean useMaterialParameterMatrix() {
    return terms.containsKey(MATERIALPARAMETER_OPTION_ID);
  }
  
  public boolean noTimeAveraging() {
    return optimisation.containsKey(NO_TIME_AVG_OPTION_ID);
  }
  
  public boolean patchwiseAdjust() {
    return optimisation.containsKey(PATCHWISE_ADJUST_OPTION_ID);
  }
  
  public boolean tempVarsOnStack() {
    return optimisation.containsKey(TEMP_VARS_ON_STACK_OPTION_ID);
  }
  
  public int maxPicardIterations() {
    if(useMaxPicardIterations()) {
      return optimisation.get(MAX_PICARD_ITER_ID);
    }
    return -1;
  }
  
  public boolean useMaxPicardIterations() {
    return optimisation.containsKey(MAX_PICARD_ITER_ID) &&
           optimisation.get(MAX_PICARD_ITER_ID)!=-1;
  }
  
  public boolean useConverterDebug() {
    return optimisation.containsKey(CONVERTER_OPTION_ID);
  }
  
  public boolean useFlopsDebug() {
    return optimisation.containsKey(FLOPS_OPTION_ID);
  }
  
  public int getNumberOfPointSources() {
    if(usePointSources()) {
      return terms.get(POINTSOURCES_OPTION_ID);
    }
    return -1;
  }
  
  //Used set the GhostLayerWidth for LimitingSolver
  public void setGhostLayerWidth(int glw) {
    this.ghostLayerWidth = glw;
  }
  
  public int getGhostLayerWidth() {
    return ghostLayerWidth; // -1 by default
  }
  
  //Used set the numberOfObservable for LimitingSolver
  public void setNumberOfObservables(int obs) {
    this.numberOfObservables = obs;
  }
  
  public int getNumberOfObservables() {
    return numberOfObservables; // -1 by default
  }
  
  // useLimiter only if the two limiter parameter have been set
  public boolean useLimiter() {
    return ghostLayerWidth > -1 && numberOfObservables > -1;
  }

  //(type: [...], terms: [...], opt: [...])
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("(type: [");
    for(String s : type) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], terms: [");
    for(String s : terms.keySet()) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], opt: [");
    for(String s : optimisation.keySet()) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("])");
    
    return sb.toString();
  }
  
}
