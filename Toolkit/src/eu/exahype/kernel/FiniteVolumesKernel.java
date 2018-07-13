package eu.exahype.kernel;

import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.stream.Collectors;

import eu.exahype.node.PIds;
import eu.exahype.node.AIds;
import eu.exahype.node.AIdentifierId;
import eu.exahype.node.AIdentifierWithMultId;

import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.PSolver;

/**
 * Circumscribes a Finite Volume Kernel
 * 
 * For a description how to add new variants consult SolverFactory.
 */
public class FiniteVolumesKernel {
  
  /**
   * Configuration parameter: id of the options
   */
  private static final Map<String, String> TYPE_OPTION_IDS = new HashMap<String, String>();
  static {
    TYPE_OPTION_IDS.put("GODUNOV_OPTION_ID",      "godunov");
    TYPE_OPTION_IDS.put("MUSCL_OPTION_ID",        "musclhancock");
    //TYPE_OPTION_IDS.put("USER_DEFINED_OPTION_ID", "user");
    TYPE_OPTION_IDS.put("LEGENDRE_OPTION_ID",     "Legendre");
    TYPE_OPTION_IDS.put("LOBATTO_OPTION_ID",      "Lobatto");
  }
  
  private static final Map<String, String> TERM_OPTION_IDS = new HashMap<String, String>();
  static {
    TERM_OPTION_IDS.put("FLUX_OPTION_ID",              "flux");
    TERM_OPTION_IDS.put("SOURCE_OPTION_ID",            "source");
    TERM_OPTION_IDS.put("NCP_OPTION_ID",               "ncp");
    TERM_OPTION_IDS.put("POINTSOURCES_OPTION_ID",      "pointsources");
    //TERM_OPTION_IDS.put("MATERIALPARAMETER_OPTION_ID", "materialparameters"); //TODO
  }
  
  private static final Map<String, String> OPTIMIZATION_OPTION_IDS = new HashMap<String, String>();
  static {
    OPTIMIZATION_OPTION_IDS.put("GENERIC_OPTION_ID",            "generic");
    OPTIMIZATION_OPTION_IDS.put("OPTIMIZED_OPTION_ID",          "optimised");
    //OPTIMIZATION_OPTION_IDS.put("PATCHWISE_ADJUST_OPTION_ID",   "patchwiseadjust"); //TODO
    OPTIMIZATION_OPTION_IDS.put("TEMP_VARS_ON_STACK_OPTION_ID", "usestack");
  }
  
  /** 
   * Maps of parsed option
   * associated int is -1 when not present in the spec file 
   */
  private Map<String, Integer> type;
  private Map<String, Integer> terms;
  private Map<String, Integer> optimization;
  
  public FiniteVolumesKernel(PSolver solver) throws IllegalArgumentException {
    if(solver instanceof AFiniteVolumesSolver) {
      type = parseIdsToMap(((AFiniteVolumesSolver) solver).getKernelType());
      terms = parseIdsToMap(((AFiniteVolumesSolver) solver).getKernelTerms());
      optimization = parseIdsToMap(((AFiniteVolumesSolver) solver).getKernelOpt());
    } else if(solver instanceof ALimitingAderdgSolver) {
      type = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelLimiterType()); 
      // Does not differ between ADER-DG solver and FV limiter 
      terms = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelTerms());
      optimization = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelLimiterOpt());
    } else {
      throw new IllegalArgumentException("No kernel definition found");
    }
	  
    validate();
  }
  
  //return null on error, use only after the program should already have failed with invalid kernel
  public static FiniteVolumesKernel noExceptionContructor(PSolver solver) {
    try {
      return new FiniteVolumesKernel(solver);
    } catch(IllegalArgumentException e) {
      return null;
    }
  }
  
  private static Map<String, Integer> parseIdsToMap(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().collect(Collectors.toMap(
      (e -> (e instanceof AIdentifierId) ? ((AIdentifierId)e).getValue().getText() : ((AIdentifierWithMultId)e).getValue().getText()),
      (e -> (e instanceof AIdentifierId) ? Integer.valueOf(-1) : Integer.valueOf(((AIdentifierWithMultId)e).getMultiplicity().getText()))
    ));
  }
  
  private void validate() throws IllegalArgumentException {
     //check if all parsed arguments are recognized
    for(String parsedId : type.keySet()) {
      if(!TYPE_OPTION_IDS.containsValue(parsedId)) {
        throw new IllegalArgumentException("Type key \""+parsedId+"\" not recognized");
      }
    }
    for(String parsedId : terms.keySet()) {
      if(!TERM_OPTION_IDS.containsValue(parsedId)) {
        throw new IllegalArgumentException("Terms key \""+parsedId+"\" not recognized");
      }
    }
    for(String parsedId : optimization.keySet()) {
      if(!OPTIMIZATION_OPTION_IDS.containsValue(parsedId)) {
        throw new IllegalArgumentException("Optimisation key \""+parsedId+"\" not recognized");
      }
    }
  }
  
  public enum KernelType {
    GenericMUSCLHancock,
    GenericGodunov,
    UserDefined,
    Unknown
  }
  
  public KernelType getKernelType() {
    if (type.containsKey(TYPE_OPTION_IDS.get("MUSCL_OPTION_ID")) && optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("OPTIMIZED_OPTION_ID"))) {
      return  KernelType.Unknown;
    }
    if (type.containsKey(TYPE_OPTION_IDS.get("MUSCL_OPTION_ID")) && optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("GENERIC_OPTION_ID"))) {
      return  KernelType.GenericMUSCLHancock;
    }
    if (type.containsKey(TYPE_OPTION_IDS.get("GODUNOV_OPTION_ID")) &&  optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("OPTIMIZED_OPTION_ID"))) {
      return  KernelType.Unknown;
    }
    if (type.containsKey(TYPE_OPTION_IDS.get("GODUNOV_OPTION_ID")) && optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("GENERIC_OPTION_ID"))) {
      return  KernelType.GenericGodunov;
    }
    //if ( type.contains(USER_DEFINED_OPTION_ID) ) {
    //  return  KernelType.UserDefined; //Not supported anymore
    //}
		
    return KernelType.Unknown;
  }
  
  public int getGhostLayerWidth() throws UnsupportedOperationException, IllegalArgumentException {
    switch (getKernelType()) {
      case GenericMUSCLHancock:
        return 2;
      case GenericGodunov:
        return 1;

    }
    throw new IllegalArgumentException("Kerneltype not recognized");
  }
  
  public boolean useGaussLobatto() {
    return type.containsKey(TYPE_OPTION_IDS.get("LOBATTO_OPTION_ID"));
  }

  public boolean usesOptimisedKernels() {
    return false;
  }
  
   public boolean useFlux() {
    return terms.containsKey(TERM_OPTION_IDS.get("FLUX_OPTION_ID"));
  }
  
  public boolean useSource() {
    return terms.containsKey(TERM_OPTION_IDS.get("SOURCE_OPTION_ID"));
  }
  
  public boolean useNCP() {
    return terms.containsKey(TERM_OPTION_IDS.get("NCP_OPTION_ID"));
  }
  
  public boolean usePointSources() {
    return terms.containsKey(TERM_OPTION_IDS.get("POINTSOURCES_OPTION_ID"));
  }
  
  public boolean useMaterialParameterMatrix() {
    return terms.containsKey(TERM_OPTION_IDS.get("MATERIALPARAMETER_OPTION_ID"));
  }
  
  public boolean tempVarsOnStack() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("TEMP_VARS_ON_STACK_OPTION_ID"));
  }
    
  //(type: [...], terms: [...], opt: [...])
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("(type: [");
    for(String s : type.keySet()) {
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
    for(String s : optimization.keySet()) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("])");
    
    return sb.toString();
  }
  
}
