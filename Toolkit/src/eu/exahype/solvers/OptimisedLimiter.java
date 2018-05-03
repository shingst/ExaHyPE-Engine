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
import eu.exahype.CodeGeneratorHelper;


public class OptimisedLimiter implements Solver {
  //Internal states
  //--------------- 
  private String         solverName;
  private Context        context;
  private TemplateEngine templateEngine;
  
  public OptimisedLimiter(String projectName, String solverName, Solver ADERDGSolver, Solver FVSolver, boolean countFlops) 
      throws IOException, IllegalArgumentException {    
    
    this.solverName                 = solverName;
    
    final String aderdgSolverName   = ADERDGSolver.getSolverName();
    final String optKernelPath      = CodeGeneratorHelper.getInstance().getIncludePath(projectName, aderdgSolverName);
    final String optNamespace       = CodeGeneratorHelper.getInstance().getNamespace(projectName, aderdgSolverName);
    
    if(optKernelPath == null || optNamespace == null) {
      throw new IllegalArgumentException("Optimized code not found!");
    }
    
    templateEngine = new TemplateEngine();
    context = new Context();
    
    //String
    context.put("project"             , projectName);
    context.put("solver"              , solverName);
    context.put("abstractSolver"      , getAbstractSolverName());
    context.put("ADERDGAbstractSolver", ADERDGSolver.getAbstractSolverName());
    context.put("FVAbstractSolver"    , FVSolver.getAbstractSolverName());
    context.put("optKernelPath"       , optKernelPath);
    context.put("optNamespace"        , optNamespace);
    
    //boolean
    context.put("countFlops"          , countFlops);
  }
    
  @Override
  public String getSolverName() {
    return solverName;
  }
  
  @Override
  public void writeHeader(java.io.BufferedWriter writer) throws IOException, IllegalArgumentException {
    throw new IllegalArgumentException("eu.exahype.solvers.OptimisedLimiter::writeHeader should not be called"); //No user implementation required
  }
  
  @Override
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
      throw new IllegalArgumentException("eu.exahype.solvers.OptimisedLimiter::writeHeader should not be called"); //No user implementation required
  }
  
  @Override
  public void writeAbstractHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {      
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractOptimisedLimiterSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractOptimisedLimiterSolverImplementation.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public boolean supportsVariables() {
    return false;
  }
}
