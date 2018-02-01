package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.IllegalArgumentException;
import java.util.Set;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

// template engine
import minitemp.Context;
import minitemp.TemplateEngine;

import eu.exahype.kernel.ADERDGKernel;
import eu.exahype.io.IOUtils;
import eu.exahype.io.SourceTemplate;

public class GenericADERDG implements Solver {
  private String _solverName;
  private Context context;
  private TemplateEngine templateEngine;

  public GenericADERDG(String projectName, String solverName, int dimensions, int numberOfVariables, int numberOfParameters, Set<String> namingSchemeNames,
      int order, boolean enableProfiler, boolean hasConstants, boolean isFortran, ADERDGKernel kernel) {

    _solverName         = solverName;
    
    final boolean isLinear           = kernel.isLinear();
    final boolean useFlux            = kernel.useFlux();
    final boolean useSource          = kernel.useSource();
    final boolean useNCP             = kernel.useNCP();
    final boolean usePointSources    = kernel.usePointSources();
    final boolean useMaterialParam   = kernel.useMaterialParameterMatrix();
    final boolean noTimeAveraging    = kernel.noTimeAveraging();
    final boolean patchwiseAdjust    = kernel.patchwiseAdjust();
    final int numberOfPointSources   = kernel.getNumberOfPointSources();
    
    templateEngine = new TemplateEngine();
    context = new Context();
    
    //String
    context.put("project"           , projectName);
    context.put("solver"            , solverName);
    context.put("abstractSolver"    , getAbstractSolverName());
    context.put("linearOrNonlinear" , isLinear? "Linear" : "Nonlinear");
    context.put("language"          , isFortran? "fortran" : "c");
    
    //int
    context.put("dimensions"        , dimensions);
    context.put("order"             , order);
    context.put("numberOfVariables" , numberOfVariables);
    context.put("numberOfParameters", numberOfParameters);
    context.put("numberOfPointSources", numberOfPointSources);
    
    //boolean
    context.put("enableProfiler"    , enableProfiler);
    context.put("hasConstants"      , hasConstants);
    context.put("isLinear"          , isLinear);
    context.put("isFortran"         , isFortran);
    context.put("useFlux"           , useFlux);
    context.put("useSource"         , useSource);
    context.put("useNCP"            , useNCP);
    context.put("usePointSources"   , usePointSources);
    context.put("useMaterialParam"  , useMaterialParam);
    context.put("noTimeAveraging"   , noTimeAveraging);
    context.put("patchwiseAdjust"   , patchwiseAdjust);

    
    //boolean as String
    context.put("useFlux_s"         , boolToTemplate(useFlux));
    context.put("useSource_s"       , boolToTemplate(useSource));
    context.put("useNCP_s"          , boolToTemplate(useNCP));
    
    //Set<String>
    context.put("namingSchemes"     , namingSchemeNames.stream().map(s -> s.substring(0, 1).toUpperCase()+s.substring(1)).collect(Collectors.toSet())); //capitalize
    
    //List<Integer> , range used by for loops
    context.put("range_0_nDim"      , IntStream.range(0, dimensions)                          .boxed().collect(Collectors.toList()));
    context.put("range_0_nVar"      , IntStream.range(0, numberOfVariables)                   .boxed().collect(Collectors.toList()));
    context.put("range_0_nVarParam" , IntStream.range(0, numberOfVariables+numberOfParameters).boxed().collect(Collectors.toList()));
  }
  
  @Override
  public String getSolverName() {
    return _solverName;
  }
  
  private String boolToTemplate(boolean b) {
    return b? "true" : "false";
  }

  @Override
  public void writeHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericADERDGSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericADERDGSolverInCUserCode.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/AbstractGenericADERDGSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractImplementation(BufferedWriter writer) throws IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/AbstractGenericADERDGSolverImplementation.template");
    writer.write(templateEngine.render(template, context));
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

  public boolean supportsVariables() {
    return true;
  }
}
