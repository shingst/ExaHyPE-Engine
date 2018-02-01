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


public class Limiter implements Solver {
  //Internal states
  //--------------- 
  private String         solverName;
  private Context        context;
  private TemplateEngine templateEngine;
  
  public Limiter(String projectName, String solverName, Solver ADERDGSolver, Solver FVSolver) 
      throws IOException, IllegalArgumentException {    
    
    this.solverName                 = solverName;
    
    templateEngine = new TemplateEngine();
    context = new Context();
    
    //String
    context.put("project"             , projectName);
    context.put("solver"              , solverName);
    context.put("abstractSolver"      , getAbstractSolverName());
    context.put("ADERDGAbstractSolver", ADERDGSolver.getAbstractSolverName());
    context.put("FVAbstractSolver"    , FVSolver.getAbstractSolverName());
  }
    
  @Override
  public String getSolverName() {
    return solverName;
  }
  
  @Override
  public void writeHeader(java.io.BufferedWriter writer) throws IOException, IllegalArgumentException {
    throw new IllegalArgumentException("eu.exahype.solvers.Limiter::writeHeader should not be called"); //No user implementation required
  }
  
  @Override
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
      throw new IllegalArgumentException("eu.exahype.solvers.Limiter::writeHeader should not be called"); //No user implementation required
  }
  
  @Override
  public void writeAbstractHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {      
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractLimiterSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractLimiterSolverImplementation.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public boolean supportsVariables() {
    return false;
  }
}
