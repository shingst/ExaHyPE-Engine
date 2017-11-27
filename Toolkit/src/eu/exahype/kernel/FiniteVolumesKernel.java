package eu.exahype.kernel;

import java.util.Set;
import java.util.stream.Collectors;

import eu.exahype.node.PIds;
import eu.exahype.node.AIds;
import eu.exahype.node.AIdentifierId;

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
  public static final String GENERIC_OPTION_ID      = "generic";
  public static final String OPTIMISED_OPTION_ID    = "optimised";
  
  public static final String GODUNOV_OPTION_ID      = "godunov";
  public static final String MUSCL_OPTION_ID        = "musclhancock";
  public static final String USER_DEFINED_OPTION_ID = "user";
  
  public static final String FLUX_OPTION_ID         = "flux";
  public static final String SOURCE_OPTION_ID       = "source";
  public static final String NCP_OPTION_ID          = "ncp";
  public static final String POINTSOURCES_OPTION_ID = "pointsources";
  
  private Set<String> type;
  private Set<String> terms;
  private Set<String> optimization;
  
  public FiniteVolumesKernel(PSolver solver) throws IllegalArgumentException {
    if(solver instanceof AFiniteVolumesSolver) {
      type = parseIds(((AFiniteVolumesSolver) solver).getKernelType());
      terms = parseIds(((AFiniteVolumesSolver) solver).getKernelTerms());
      optimization = parseIds(((AFiniteVolumesSolver) solver).getKernelOpt());
    } else if(solver instanceof ALimitingAderdgSolver) {
      type = parseIds(((ALimitingAderdgSolver) solver).getKernelLimiterType()); 
      // Does not differ between ADER-DG solver and FV limiter 
      terms = parseIds(((ALimitingAderdgSolver) solver).getKernelTerms());
      optimization = parseIds(((ALimitingAderdgSolver) solver).getKernelLimiterOpt());
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
  
  private static Set<String> parseIds(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().map(e -> ((AIdentifierId)e).getValue().getText()).collect(Collectors.toSet());
  }
  
  private void validate() throws IllegalArgumentException {
  }
  
  public enum KernelType {
    GenericMUSCLHancock,
    GenericGodunov,
    UserDefined,
    Unknown
  }
  
  public KernelType getKernelType() {
	if ( type.contains(MUSCL_OPTION_ID) && optimization.contains(OPTIMISED_OPTION_ID) ) {
      return  KernelType.Unknown;
	}
	if ( 
      type.contains(MUSCL_OPTION_ID) && optimization.contains(GENERIC_OPTION_ID)
      || 
      type.contains(MUSCL_OPTION_ID)
    ) {
	  return  KernelType.GenericMUSCLHancock;
	}
	if ( type.contains(GODUNOV_OPTION_ID) && optimization.contains(OPTIMISED_OPTION_ID) ) {
      return  KernelType.Unknown;
	}
	if ( 
      type.contains(GODUNOV_OPTION_ID) && optimization.contains(GENERIC_OPTION_ID)
      || 
      type.contains(GODUNOV_OPTION_ID)
    ) {
	  return  KernelType.GenericGodunov;
	}
	if ( type.contains(USER_DEFINED_OPTION_ID) ) {
      return  KernelType.UserDefined;
    }
		
    return KernelType.Unknown;
  }

  public boolean usesOptimisedKernels() {
	return false;
  }
  
  public boolean useFlux() {
    return terms.contains(FLUX_OPTION_ID);
  }
  
  
  public boolean useSource() {
    return terms.contains(SOURCE_OPTION_ID);
  }
  
  public boolean useNCP() {
    return terms.contains(NCP_OPTION_ID);
  }
  
  public boolean usePointSources() {
    return terms.contains(POINTSOURCES_OPTION_ID);
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
    for(String s : terms) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], opt: [");
    for(String s : optimization) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("])");
    
    return sb.toString();
  }
  
}