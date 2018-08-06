package eu.exahype;

// template engine
import minitemp.Context;
import minitemp.TemplateEngine;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.kernel.ADERDGKernel;
import eu.exahype.node.*;
import eu.exahype.io.IOUtils;

public class SetupBuildEnvironment extends DepthFirstAdapter {
  public Boolean valid = true;

  private DirectoryAndPathChecker directoryAndPathChecker;
  private Context                 context;
  private TemplateEngine          templateEngine;
  
  public SetupBuildEnvironment(DirectoryAndPathChecker directoryAndPathChecker) {
    this.directoryAndPathChecker = directoryAndPathChecker;
    templateEngine               = new TemplateEngine();
    context                      = new Context();
    
    //default context value
    context.put("useSharedMem",      false);
    context.put("useDistributedMem", false);
    context.put("useFortran",        false);
    context.put("useOptKernel",      false);
    context.put("useLikwid",         false);
    context.put("useIpcm",           false);
  }

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    int dimensions = Integer.parseInt( node.getDimension().toString().trim() );
    if (dimensions==2 || dimensions==3) {
      System.out.print(dimensions+"d experiment ... ok\n");
      context.put("dimensions", dimensions);
    }
    else {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
      valid = false;
    }
  }

  @Override
  public void inASharedMemory(ASharedMemory node) {
    context.put("useSharedMem", true);
    System.out.print("shared memory ... TBB (switch to OpenMP manually as indicated below)\n");
    if (!System.getenv().containsKey("TBB_INC")) {
      System.out.print(
          "WARNING: environment variable TBB_INC not set but required if code is built with TBB\n");
    }
    if (!System.getenv().containsKey("TBB_SHLIB")) {
      System.out.print(
          "WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB\n");
    }
  }

  @Override
  public void inADistributedMemory(ADistributedMemory node) {
    context.put("useDistributedMem", true);
    System.out.print("mpi ... switched on \n");
  }

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    handleAderdgSolver(node.getLanguage().getText().trim(), node, node.getName().getText().trim()); 
  }
  
  @Override
  public void inALimitingAderdgSolver(ALimitingAderdgSolver node) {
    handleAderdgSolver(node.getLanguage().getText().trim(), node, node.getName().getText().trim());
  }

  private void handleAderdgSolver(String language, PSolver node, String name) {
    if (language.equals("C")) {
      //do nothing, default behavior
    } else if (language.equals("Fortran")) {
      context.put("useFortran",   true);
    } else {
      System.err.println("ERROR: unknown language for solver " + name
          + ". Supported languages are C and Fortran");
      valid = false;
    }
    
    ADERDGKernel kernel = ADERDGKernel.noExceptionContructor(node);
    if(language.equals("C") && (kernel.usesOptimisedKernels())) {
      context.put("useOptKernel", true);
    }
  }
  
  @Override
  public void inAProfiling(AProfiling node) {
    if (node.getLikwidInc() != null) {
      context.put("useLikwid", true);
      context.put("likwidInc", node.getLikwidInc().toString().trim());
    }
    if(node.getLikwidLib() != null) {
      context.put("likwidLib", node.getLikwidLib().toString().trim());
    }

    if (node.getIpcmInc() != null) {
      context.put("useIpcm", true);
      context.put("ipcmInc", node.getIpcmInc().toString().trim());
    }
    if (node.getIpcmLib() != null) {
      context.put("ipcmLib", node.getIpcmLib().toString().trim());
    }
  }

  @Override
  public void outAProject(AProject node) {
    if(!valid)
      return; // don't try to write if not valid
    try {      
      context.put("project",           node.getName());
      context.put("executableName",   "ExaHyPE-" + node.getName());
      
      //paths
      context.put("peanoToolboxPath", directoryAndPathChecker.peanoKernelPath.getCanonicalPath());
      context.put("exahypePath",      directoryAndPathChecker.exahypePath.getCanonicalPath());
      context.put("outputPath",       directoryAndPathChecker.outputDirectory.getCanonicalPath());
      System.out.print("store pathes and default settings in makefile ... ok");
    
      String architecture = "noarch";
      if (node.getArchitecture()!=null) {
        architecture = node.getArchitecture().toString().trim().toLowerCase();
      }
      context.put("architecture",     architecture);
      
      int alignment = 16;
      if (architecture.equals("snb")) {
        alignment=32;
      } else if (architecture.equals("hsw")) {
        alignment=32;
      } else if (architecture.equals("knc")) {
        alignment=64;
      } else if (architecture.equals("knl")) {
        alignment=64;
      }
      context.put("alignment",        alignment);
    
      // echo toolkit help message
      final String messageTemplate = IOUtils.convertRessourceContentToString("eu/exahype/buildEnvironment/templates/message.template");
      System.out.println("\n\n\n\n");
      System.out.print(templateEngine.render(messageTemplate, context));

      
      // write Makefile
      java.io.BufferedWriter writer = new java.io.BufferedWriter(new java.io.FileWriter(new java.io.File(directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/Makefile")));
      final String template = IOUtils.convertRessourceContentToString("eu/exahype/buildEnvironment/templates/Makefile.template");
      writer.write(templateEngine.render(template, context));
      writer.close();
      
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }
}
