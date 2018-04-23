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
  private static final Map<String, String> TYPE_OPTION_ID = new HashMap<String, String>();
  static {
    TYPE_OPTION_ID.put("GODUNOV_OPTION_ID",      "godunov");
    TYPE_OPTION_ID.put("MUSCL_OPTION_ID",        "musclhancock");
    //TYPE_OPTION_ID.put("USER_DEFINED_OPTION_ID", "user");
    TYPE_OPTION_ID.put("LEGENDRE_OPTION_ID",     "Legendre");
    TYPE_OPTION_ID.put("LOBATTO_OPTION_ID",      "Lobatto");
  }
  
  private static final Map<String, String> TERMS_OPTION_ID = new HashMap<String, String>();
  static {
    TERMS_OPTION_ID.put("FLUX_OPTION_ID",              "flux");
    TERMS_OPTION_ID.put("SOURCE_OPTION_ID",            "source");
    TERMS_OPTION_ID.put("NCP_OPTION_ID",               "ncp");
    TERMS_OPTION_ID.put("POINTSOURCES_OPTION_ID",      "pointsources");
    //TERMS_OPTION_ID.put("MATERIALPARAMETER_OPTION_ID", "materialparameters"); //TODO
  }
  
  private static final Map<String, String> OPTIMIZATION_OPTION_ID = new HashMap<String, String>();
  static {
    OPTIMIZATION_OPTION_ID.put("GENERIC_OPTION_ID",            "generic");
    OPTIMIZATION_OPTION_ID.put("OPTIMIZED_OPTION_ID",          "optimised");
    //OPTIMIZATION_OPTION_ID.put("PATCHWISE_ADJUST_OPTION_ID",   "patchwiseadjust"); //TODO
    OPTIMIZATION_OPTION_ID.put("TEMP_VARS_ON_STACK_OPTION_ID", "usestack");
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
      if(!TYPE_OPTION_ID.containsValue(parsedId)) {
        throw new IllegalArgumentException("Type key \""+parsedId+"\" not recognized");
      }
    }
    for(String parsedId : terms.keySet()) {
      if(!TERMS_OPTION_ID.containsValue(parsedId)) {
        throw new IllegalArgumentException("Terms key \""+parsedId+"\" not recognized");
      }
    }
    for(String parsedId : optimization.keySet()) {
      if(!OPTIMIZATION_OPTION_ID.containsValue(parsedId)) {
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
    if (type.containsKey(TYPE_OPTION_ID.get("MUSCL_OPTION_ID")) && optimization.containsKey(OPTIMIZATION_OPTION_ID.get("OPTIMIZED_OPTION_ID"))) {
      return  KernelType.Unknown;
    }
    if (type.containsKey(TYPE_OPTION_ID.get("MUSCL_OPTION_ID")) && optimization.containsKey(OPTIMIZATION_OPTION_ID.get("GENERIC_OPTION_ID"))) {
      return  KernelType.GenericMUSCLHancock;
    }
    if (type.containsKey(TYPE_OPTION_ID.get("GODUNOV_OPTION_ID")) &&  optimization.containsKey(OPTIMIZATION_OPTION_ID.get("OPTIMIZED_OPTION_ID"))) {
      return  KernelType.Unknown;
    }
    if (type.containsKey(TYPE_OPTION_ID.get("GODUNOV_OPTION_ID")) && optimization.containsKey(OPTIMIZATION_OPTION_ID.get("GENERIC_OPTION_ID"))) {
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

  public boolean usesOptimisedKernels() {
    return false;
  }
  
   public boolean useFlux() {
    return terms.containsKey(TERMS_OPTION_ID.get("FLUX_OPTION_ID"));
  }
  
  public boolean useSource() {
    return terms.containsKey(TERMS_OPTION_ID.get("SOURCE_OPTION_ID"));
  }
  
  public boolean useNCP() {
    return terms.containsKey(TERMS_OPTION_ID.get("NCP_OPTION_ID"));
  }
  
  public boolean usePointSources() {
    return terms.containsKey(TERMS_OPTION_ID.get("POINTSOURCES_OPTION_ID"));
  }
  
  public boolean useMaterialParameterMatrix() {
    return terms.containsKey(TERMS_OPTION_ID.get("MATERIALPARAMETER_OPTION_ID"));
  }
  
  public boolean tempVarsOnStack() {
    return optimization.containsKey(OPTIMIZATION_OPTION_ID.get("TEMP_VARS_ON_STACK_OPTION_ID"));
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
