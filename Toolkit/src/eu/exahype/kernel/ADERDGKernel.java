package eu.exahype.kernel;

import java.util.Map;
import java.util.HashMap;
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
  private static final Map<String, String> TYPE_OPTION_IDS = new HashMap<String, String>();
  static {
    TYPE_OPTION_IDS.put("LINEAR_OPTION_ID",     "linear");
    TYPE_OPTION_IDS.put("NONLINEAR_OPTION_ID",  "nonlinear");
    TYPE_OPTION_IDS.put("LEGENDRE_OPTION_ID",   "Legendre");
    TYPE_OPTION_IDS.put("LOBATTO_OPTION_ID",    "Lobatto");
  }
  
  private static final Map<String, String> TERM_OPTION_IDS = new HashMap<String, String>();
  static {
    TERM_OPTION_IDS.put("FLUX_OPTION_ID",              "flux");
    TERM_OPTION_IDS.put("SOURCE_OPTION_ID",            "source");
    TERM_OPTION_IDS.put("NCP_OPTION_ID",               "ncp");
    TERM_OPTION_IDS.put("POINTSOURCES_OPTION_ID",      "pointsources");
    TERM_OPTION_IDS.put("MATERIALPARAMETER_OPTION_ID", "materialparameters");
    TERM_OPTION_IDS.put("PARABOLIC_FLUX_OPTION_ID",     "parabolicflux");
  }
  
  private static final Map<String, String> OPTIMIZATION_OPTION_IDS = new HashMap<String, String>();
  static {
    OPTIMIZATION_OPTION_IDS.put("GENERIC_OPTION_ID",            "generic");
    OPTIMIZATION_OPTION_IDS.put("OPTIMIZED_OPTION_ID",          "optimised");
    OPTIMIZATION_OPTION_IDS.put("NO_TIME_AVG_OPTION_ID",        "notimeavg");
    OPTIMIZATION_OPTION_IDS.put("PATCHWISE_ADJUST_OPTION_ID",   "patchwiseadjust");
    OPTIMIZATION_OPTION_IDS.put("TEMP_VARS_ON_STACK_OPTION_ID", "usestack");
    OPTIMIZATION_OPTION_IDS.put("MAX_PICARD_ITER_OPTION_ID",    "maxpicarditer");
    OPTIMIZATION_OPTION_IDS.put("FUSEDSOURCE_OPTION_ID",        "fusedsource");
    OPTIMIZATION_OPTION_IDS.put("FLUX_VECT_OPTION_ID",          "fluxvect");
    OPTIMIZATION_OPTION_IDS.put("FUSEDSOURCE_VECT_OPTION_ID",   "fusedsourcevect");
    OPTIMIZATION_OPTION_IDS.put("CERK_GUESS_OPTION_ID",         "cerkguess");
    OPTIMIZATION_OPTION_IDS.put("CONVERTER_OPTION_ID",          "converter"); //for debug only, not in guidebook
    OPTIMIZATION_OPTION_IDS.put("FLOPS_OPTION_ID",              "flops");     //for debug only, not in guidebook
  }

  /** 
   * Maps of parsed option
   * associated int is -1 when not present in the spec file 
   */
  private Map<String, Integer> type;
  private Map<String, Integer> terms;
  private Map<String, Integer> optimization;
  
  /**
   * numberOfObservables initialized if using the a LimitingSolver
   */ 
  private int numberOfObservables = -1; // TODO This should not be part of this class!
  
  
  
  public ADERDGKernel(PSolver solver) throws IllegalArgumentException {
    if(solver instanceof AAderdgSolver) {
      type = parseIdsToMap(((AAderdgSolver) solver).getKernelType());
      terms = parseIdsToMap(((AAderdgSolver) solver).getKernelTerms());
      optimization =  parseIdsToMap(((AAderdgSolver) solver).getKernelOpt());
    } else if(solver instanceof ALimitingAderdgSolver) {
      type = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelType());
      terms = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelTerms());
      optimization =  parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelOpt());
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
  
  private static Map<String, Integer> parseIdsToMap(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().collect(Collectors.toMap(
      (e -> (e instanceof AIdentifierId) ? ((AIdentifierId)e).getValue().getText() : ((AIdentifierWithMultId)e).getValue().getText()),
      (e -> (e instanceof AIdentifierId) ? Integer.valueOf(-1) : Integer.valueOf(((AIdentifierWithMultId)e).getMultiplicity().getText()))
    ));
  }
  
  private void validate() throws IllegalArgumentException {
    // check if all parsed arguments are recognized
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
        throw new IllegalArgumentException("Optimization key \""+parsedId+"\" not recognized");
      }
    }
    // must be linear xor nonlinear
    if(!type.containsKey(TYPE_OPTION_IDS.get("LINEAR_OPTION_ID")) ^ type.containsKey(TYPE_OPTION_IDS.get("NONLINEAR_OPTION_ID"))) {//should be only one
      throw new IllegalArgumentException("nonlinear or linear not specified or both specified in the kernel type");
    }
    // pointsource requires an associated int value
    if(usePointSources() && getNumberOfPointSources() < 0) {
      throw new IllegalArgumentException("point sources used but number not specified! In the specification file, use "+TERM_OPTION_IDS.get("POINTSOURCES_OPTION_ID")+":X, with X the number of point sources.");
    }
    // fusedsource requires optimized kernels
    if(useFusedSource() && !(getKernelType() ==  KernelType.OptimisedADERDG)) {
      throw new IllegalArgumentException("The optimization '"+OPTIMIZATION_OPTION_IDS.get("FUSEDSOURCE_OPTION_ID")+"' requires the used of optimized kernels");
    }
    // can't used fusedSource without source
    if(useFusedSource() && !useSource()) {
      throw new IllegalArgumentException("The optimization '"+OPTIMIZATION_OPTION_IDS.get("FUSEDSOURCE_OPTION_ID")+"' requires the PDE term '"+TERM_OPTION_IDS.get("SOURCE_OPTION_ID")+"'");
    }
    // fusedsourcevect include fusedsource
    if(useFusedSource() && !useSource()) {
      throw new IllegalArgumentException("The optimization '"+OPTIMIZATION_OPTION_IDS.get("FUSEDSOURCE_VECT_OPTION_ID")+"' already includes the optimization '"+TERM_OPTION_IDS.get("FUSEDSOURCE_OPTION_ID")+"'");
    }
    // vect PDEs only with otptimized kernels
    if(!(getKernelType() ==  KernelType.OptimisedADERDG) && (useFluxVect())) {//TODO JMG extend with other vect PDE
      throw new IllegalArgumentException("The vectorized PDE optimizations require the optimized kernel term (use '"+OPTIMIZATION_OPTION_IDS.get("OPTIMIZED_OPTION_ID")+"')");
    }
    // can't use opt fluxvect without term flux
    if(useFluxVect() && !useFlux()) {
      throw new IllegalArgumentException("The optimization '"+OPTIMIZATION_OPTION_IDS.get("FLUX_VECT_OPTION_ID")+"' requires the PDE term '"+TERM_OPTION_IDS.get("FLUX_OPTION_ID")+"'");
    }
    // can't use flux and parabolic flux together
    if(useFlux() && useParabolicFlux()) {
<<<<<<< e409deff4e9c56888ffbccaed0916189314ffa50
      throw new IllegalArgumentException("The term '"+TERM_OPTION_IDS.get("FLUX_OPTION_ID")+"' and term '"+TERM_OPTION_IDS.get("PARABOLIC_FLUX_OPTION_ID")+"' can't be used together");
=======
      throw new IllegalArgumentException("The term '"+TERMS_OPTION_ID.get("FLUX_OPTION_ID")+"' and term '"+TERMS_OPTION_ID.get("PARABOLIC_FLUX_OPTION_ID")+"' can't be used together");
>>>>>>> introduce parabolicflux flag in Toolkit
    }
  }

  public enum KernelType {
    GenericADERDG,
    OptimisedADERDG,
    Unknown
  }

  public boolean isLinear() throws IllegalArgumentException {
    return type.containsKey(TYPE_OPTION_IDS.get("LINEAR_OPTION_ID"));
  }

  public KernelType getKernelType() {
    if (optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("OPTIMIZED_OPTION_ID"))) {
      return KernelType.OptimisedADERDG;
    }
    if(optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("GENERIC_OPTION_ID"))) {
      return KernelType.GenericADERDG;
    }
    return  KernelType.Unknown;
  }
  
  public boolean useGaussLobatto() {
    return type.containsKey(TYPE_OPTION_IDS.get("LOBATTO_OPTION_ID"));
  }

  public boolean usesOptimisedKernels() {
    // assert: !optimization.containsKey(GENERIC_OPTION_ID)
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("OPTIMIZED_OPTION_ID"));
  }

  public boolean useFlux() {
    return terms.containsKey(TERM_OPTION_IDS.get("FLUX_OPTION_ID"));
  }

  public boolean useParabolicFlux() {
<<<<<<< e409deff4e9c56888ffbccaed0916189314ffa50
    return terms.containsKey(TERM_OPTION_IDS.get("PARABOLIC_FLUX_OPTION_ID"));
=======
    return terms.containsKey(TERMS_OPTION_ID.get("PARABOLIC_FLUX_OPTION_ID"));
>>>>>>> introduce parabolicflux flag in Toolkit
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
  
  public boolean noTimeAveraging() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("NO_TIME_AVG_OPTION_ID"));
  }
  
  public boolean patchwiseAdjust() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("PATCHWISE_ADJUST_OPTION_ID"));
  }
  
  public boolean tempVarsOnStack() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("TEMP_VARS_ON_STACK_OPTION_ID"));
  }
  
  public int maxPicardIterations() {
    if(useMaxPicardIterations()) {
      return optimization.get(OPTIMIZATION_OPTION_IDS.get("MAX_PICARD_ITER_OPTION_ID"));
    }
    return -1;
  }
  
  public boolean useFusedSource() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("FUSEDSOURCE_OPTION_ID"));
  }
  
  public boolean useFusedSourceVect() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("FUSEDSOURCE_VECT_OPTION_ID"));
  }
  
  public boolean useFluxVect() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("FLUX_VECT_OPTION_ID"));
  }
  
  public boolean useCERKGuess() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("CERK_GUESS_OPTION_ID"));
  }
  
  public boolean useMaxPicardIterations() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("MAX_PICARD_ITER_OPTION_ID")) &&
           optimization.get(OPTIMIZATION_OPTION_IDS.get("MAX_PICARD_ITER_OPTION_ID"))!=-1;
  }
  
  public boolean useConverterDebug() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("CONVERTER_OPTION_ID"));
  }
  
  public boolean useFlopsDebug() {
    return optimization.containsKey(OPTIMIZATION_OPTION_IDS.get("FLOPS_OPTION_ID"));
  }
  
  public int getNumberOfPointSources() {
    if(usePointSources()) {
      return terms.get(TERM_OPTION_IDS.get("POINTSOURCES_OPTION_ID"));
    }
    return -1;
  }
  
  //Used set the numberOfObservable for LimitingSolver
  public void setNumberOfObservables(int obs) { // TODO This should not be part of this class!
    this.numberOfObservables = obs;
  }
  
  public int getNumberOfObservables() { // TODO This should not be part of this class!
    return numberOfObservables; // -1 by default
  }
  
  // useLimiter only if the two limiter parameter have been set
  public boolean useLimiter() { // TODO This should not be part of this class!
    return numberOfObservables > -1;
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
