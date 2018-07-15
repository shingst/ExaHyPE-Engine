package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Set;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

// template engine
import minitemp.Context;
import minitemp.TemplateEngine;

import eu.exahype.io.IOUtils;
import eu.exahype.kernel.FiniteVolumesKernel;


public class GenericFiniteVolumesInC implements Solver {
  
  private String _solverName;
  private Context context;
  private TemplateEngine templateEngine;

  public GenericFiniteVolumesInC(
      String type, String projectName, String solverName, int dimensions, int numberOfVariables, int numberOfParameters, Set<String> namingSchemeNames, int patchSize,
      boolean enableProfiler, boolean hasConstants,
      FiniteVolumesKernel kernel) {
        
    _solverName         = solverName;
    
    final boolean useFlux            = kernel.useFlux() || kernel.useParabolicFlux();
    final boolean useSource          = kernel.useSource();
    final boolean useNCP             = kernel.useNCP();
    final boolean usePointSources    = kernel.usePointSources();
    final boolean useParabolicFlux   = kernel.useParabolicFlux();
    final boolean tempVarsOnStack    = kernel.tempVarsOnStack();
    final int     ghostLayerWidth    = kernel.getGhostLayerWidth();
    
    templateEngine = new TemplateEngine();
    context = new Context();

    //String
    context.put("finiteVolumesType" , type);
    context.put("project"           , projectName);
    context.put("solver"            , solverName);
    context.put("abstractSolver"    , getAbstractSolverName());
    
    //int
    context.put("dimensions"        , dimensions);
    context.put("numberOfVariables" , numberOfVariables);
    context.put("numberOfParameters", numberOfParameters);
    context.put("patchSize"         , patchSize);
    context.put("ghostLayerWidth"   , ghostLayerWidth);
    
    //boolean
    context.put("enableProfiler"    , enableProfiler);
    //context.put("hasConstants"      , hasConstants);
    context.put("useFlux"           , useFlux);
    context.put("useSource"         , useSource);
    context.put("useNCP"            , useNCP);
    context.put("usePointSources"   , usePointSources);
    context.put("useParabolicFlux"  , useParabolicFlux);
    context.put("tempVarsOnStack"   , tempVarsOnStack);
    
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

  public void writeHeader(java.io.BufferedWriter writer) throws IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/GenericFiniteVolumesSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  public void writeUserImplementation(java.io.BufferedWriter writer) throws IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/GenericFiniteVolumesSolverInCUserCode.template");
    writer.write(templateEngine.render(template, context));
  }

  @Override
  public void writeAbstractHeader(BufferedWriter writer) throws IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractGenericFiniteVolumesSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractImplementation(BufferedWriter writer) throws IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractGenericFiniteVolumesSolverInCImplementation.template");
    writer.write(templateEngine.render(template, context));
  }

  public void writeUserPDE(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a PDE.f90.\n");
  }


  public void writeTypesDef(java.io.BufferedWriter writer)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a typesDef.f90.\n");
  }
  
  @Override
  public boolean supportsVariables() {
    return true;
  }
}
