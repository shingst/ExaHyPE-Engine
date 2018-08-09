package eu.exahype;

import java.io.IOException;
import java.io.BufferedWriter;
import java.util.List;
import java.util.LinkedList;

import eu.exahype.io.FileSearch;
import eu.exahype.io.IOUtils;

// template engine
import minitemp.Context;
import minitemp.TemplateEngine;

public class CreateReadme {

  public boolean valid = true;
  private List<Context> solversContexts;
  
  //Singleton pattern (to be able to access the instance everywhere)
  //-----------------
  private static volatile CreateReadme instance = null;

  private CreateReadme() {
    solversContexts = new LinkedList<Context>();
  }

  public static CreateReadme getInstance() {
      if (instance == null) {
          synchronized(CreateReadme.class) {
              if (instance == null) {
                  instance = new CreateReadme();
              }
          }
      }
      return instance;
  }
  
  //Add a solver's context to the list
  //----------------------------------
  public void addSolverContext(Context context) {
    solversContexts.add(context);
  }
  
  //Generate README
  //---------------
  
  public void writeReadme(String outputPath) {
    try {      
      final TemplateEngine templateEngine = new TemplateEngine();
      final Context context = new Context();
      
      //generate the solvers README part
      final String solverPartTemplate = IOUtils.convertRessourceContentToString("eu/exahype/readmeTemplates/solverPart.template");
      final List<String> renderedSolversPart = new LinkedList<String>();
      for(Context solverContext : solversContexts) {
        fillMissingContextValues(solverContext);
        renderedSolversPart.add(templateEngine.render(solverPartTemplate, solverContext));
      }
      context.put("solverParts", renderedSolversPart);

      //generate the README
      java.io.File readmeFile = FileSearch.relocatableFile(outputPath + "/README_generated.md");
      if (readmeFile.exists()) {
        System.out.println("README file ... does exist already. Is overwritten");
      }
      final BufferedWriter writer = new BufferedWriter(new java.io.FileWriter(readmeFile));
      final String baseTemplate = IOUtils.convertRessourceContentToString("eu/exahype/readmeTemplates/base.template");
      writer.write(templateEngine.render(baseTemplate, context));
      writer.close();
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }
  
  // fill the context with default value if not already set
  private void fillMissingContextValues(Context context) {
    //boolean
    context.putIfNotSet("enableProfiler"        , false);
    context.putIfNotSet("hasConstants"          , false);
    context.putIfNotSet("isLinear"              , false);
    context.putIfNotSet("useFlux"               , false);
    context.putIfNotSet("useFluxVect"           , false);
    context.putIfNotSet("useSource"             , false);
    context.putIfNotSet("useFusedSource"        , false);
    context.putIfNotSet("useFusedSourceVect"    , false);
    context.putIfNotSet("useNCP"                , false);
    context.putIfNotSet("usePointSources"       , false);
    context.putIfNotSet("useMaterialParam"      , false);
    context.putIfNotSet("noTimeAveraging"       , false);
    context.putIfNotSet("patchwiseAdjust"       , false);
    context.putIfNotSet("tempVarsOnStack"       , false);
    context.putIfNotSet("useMaxPicardIterations", false);
    context.putIfNotSet("countFlops"            , false);
  }
}
