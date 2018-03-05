package eu.exahype;

// Command line parsing
import org.apache.commons.cli.*;

import java.util.List;

import eu.exahype.node.Start;
import eu.exahype.node.Node;

import eu.exahype.fromtosable.FromSableToStructured;

public class Main {
  boolean interactive;
  boolean verbose;

  public static void printHeader() {
    System.out.println("================================");
    System.out.println(" ___          _  _      ___ ___");
    System.out.println("/ __|_ ____ _| || |_  _/ _ \\ __|");
    System.out.println("| _|\\ \\ / _` | __ | || |  _/ _| ");
    System.out.println("\\___/_\\_\\__,_|_||_|\\_, |_| \\___|");
    System.out.println("                   |__/         ");
    System.out.println("================================");
    System.out.println("");
    System.out.println(" www.exahype.eu ");
    System.out.println("");
    System.out.println("================================");
    System.out.println("");
    System.out.println("The project has received funding from the European Union's ");
    System.out.println("Horizon 2020 research and innovation programme under grant ");
    System.out.println("agreement No 671698 (ExaHyPE). It is based upon the PDE ");
    System.out.println("framework Peano (www.peano-framework.org).");
    System.out.println("");
    System.out.println("");
  }
  
  public void printHelp(Options options) {
    String header = "\n"+
       "This is the ExaHyPE toolkit, a glue code generator for setting up and updating "+
       "an ExaHyPE application. Please refer to the ExaHyPE guidebook for usage examples. "+
       "As input, a specification file is needed.\n";
    String footer = "\nFor further information, go to http://www.exahype.eu/";

    HelpFormatter formatter = new HelpFormatter();
    //formatter.printHelp("myapp", header, options, footer, true);
    formatter.printHelp("java -jar ExaHyPE.jar [-flags] <YourSpecfile.exahype>", header, options, footer, false);
    System.exit(-2);
  }
  
  /**
   * Using apache.commons.cli Command line parsing:
   **/
  public Main(String[] args) {
    Options options = new Options();
    CommandLine line;

    options.addOption("c", "clean-opt", false, "Clean optimized kernels (only applicable if optimized kernels are used in Specfile)");
    options.addOption("n", "not-interactive", false, "Run non-interactively in non-stop mode");
    options.addOption("i", "interactive", false, "Run interactively. This is the default. [Deprecated, use --non-interactive instead]");
    options.addOption("q", "quiet", false, "Be quiet, do not say so much");
    options.addOption("h", "help", false, "Show usage information");
    options.addOption("e", "export", false, "Export specification file contents into something else, dump on stdout. Does not involve the toolkit afterwards.");

    CommandLineParser parser = new DefaultParser();

    try {
        line = parser.parse(options, args);
        if(line.hasOption("help")) printHelp(options);
        
        // Non-consumed arguments: Actual file list
        List<String> remaining_args = line.getArgList();
        if(remaining_args.size() != 1) {
          System.err.println("ERROR: Please provide input file as argument.");
          printHelp(options);
        }
      
        // backward compatibility:
        if(line.hasOption("interactive")) interactive = true;
        else if(line.hasOption("not-interactive")) interactive = false;
        else {
          System.out.println(
              "INFO: You might want to add --interactive or --non-interactive as command");
          System.out.println("      line argument to control whether script runs interactively");
          interactive = true; // default value
        }
      
        verbose = !line.hasOption("quiet");

        String inputFileName = remaining_args.get(0);
        Node document = parseFile(inputFileName);
        
        if(line.hasOption("export")) {
          exportSpecfile(document);
          System.exit(0);
        }

        println("Start to interpret script ... ");
        waitForInteraction();
        runToolkit(document, inputFileName);
    } catch(ParseException exp ) {
        System.err.println("Parsing failed.  Reason: " + exp.getMessage());
        System.exit(-1);
    } catch (Exception e) {
      System.err.println("Error: " + e.toString());
      if(verbose) e.printStackTrace();
      System.err.println("ExaHyPE Toolkit failed");
      System.exit(-2);
    }
  } // end of main function
  
  public void print(String msg) {
    if(verbose) System.out.print(msg);
  }
  public void println(String msg) {
    if(verbose) System.out.println(msg);
  }
  
  public void exportSpecfile(Node document) {
    // Begin: Sven testing JSON export.
	FromSableToStructured exporter = new FromSableToStructured();
	exporter.setIncludeMissingOptionals(true).dump(document);
	//exporter.dump();
	System.out.println("Done, Slurped everything.\n");
	
	// Exit here for testing
	return;
  }

  public Node parseFile(String inputFileName) {
    eu.exahype.parser.Parser parser = null;
    eu.exahype.node.Start document = null;
    try {
      print("read input file " + inputFileName + " .");
      parser = new eu.exahype.parser.Parser(new eu.exahype.lexer.Lexer(
          new java.io.PushbackReader(new java.io.FileReader(inputFileName),5000))); //TODO Dominic fix the PushbackReader buffer size
      document = parser.parse();
      println(".. ok");
      println("\n\n\n\n");
    } catch (Exception e) {
      println(".. failed ");
      println("\n\n\n\n");
      System.err.println("Error parsing specification file: " + e.toString());
      System.exit(-3);
    }
    return document;
  }

  public void waitForInteraction() {
    if (interactive) {
      System.out.println("<press Enter>");
      try {
        System.in.read();
      } catch (Exception e) {
      }
      for (int i = 0; i < 50; ++i) System.out.println();
      printHeader();
    }
  }

  /**
   * Run the ExaHyPE toolkit with its glue code generator.
   *
   * This works as following: The specfile is parsed with the parser generated statically
   * by SableCC. Then the AST is visited multiple times by several code generation steps.
   *
   **/
  public void runToolkit(Node document, String inputFileName) throws Exception {
    //
    // Usually, I write the header directly before a new algorithm phase, but
    // not for the first phase
    //    
    if(verbose) printHeader();

    DirectoryAndPathChecker directoryAndPathChecker = null;

    //
    // Check directories and pathes
    //
    directoryAndPathChecker = new DirectoryAndPathChecker();
    document.apply(directoryAndPathChecker);

    println("\n\n\n\n");
    if (!directoryAndPathChecker.valid) {
      System.err.println("ERROR: Some directories did not exist and/or could not be created");
      System.err.println("ExaHyPE script failed ");
      System.exit(-4);
    }
    println("validated and configured pathes ... ok");
    waitForInteraction();
    
    // Create the solvers
    CreateSolverClasses createSolverClasses = new CreateSolverClasses(directoryAndPathChecker);
    document.apply(createSolverClasses);

    println("\n\n\n\n");
    if (!createSolverClasses.valid) {
      System.err.println("ERROR: Could not create application's solver classes");
      System.err.println("ExaHyPE script failed ");
      System.exit(-6);
    }
    println("generate application-specific solver classes ... ok");
    waitForInteraction();

    // Create the plotters
    CreatePlotterClasses createPlotterClasses = new CreatePlotterClasses(directoryAndPathChecker);
    document.apply(createPlotterClasses);

    println("\n\n\n\n");
    if (!createPlotterClasses.valid) {
      System.err.println("ERROR: Could not create application's plotter classes");
      System.err.println("ExaHyPE script failed ");
      System.exit(-8);
    }
    println("generate application-specific plotter classes ... ok");
    waitForInteraction();

    // Create the kernel calls
    GenerateSolverRegistration generateKernelCalls =
        new GenerateSolverRegistration(directoryAndPathChecker, inputFileName);
    document.apply(generateKernelCalls);

    println("\n\n\n\n");
    if (!generateKernelCalls.valid) {
      System.err.println("ERROR: Could not create ExaHyPE's kernel calls");
      System.err.println("ExaHyPE script failed ");
      System.exit(-10);
    }
    println("generate computational kernel calls ... ok");
    waitForInteraction();

    //
    // Setup build environment, i.e. makefiles
    //
    SetupBuildEnvironment setupBuildEnvironment =
        new SetupBuildEnvironment(directoryAndPathChecker);
    document.apply(setupBuildEnvironment);

    println("\n\n\n\n");
    if (!setupBuildEnvironment.valid) {
      System.err.println("ERROR: Could not create ExaHyPE's build environment");
      System.err.println("ExaHyPE script failed ");
      System.exit(-12);
    }
    println("setup build environment ... ok");
    waitForInteraction();
  } // end of runToolkit()
  
  public static void main(String[] args) {
    // switch from static to class context
    Main m = new Main(args);
  }
}
