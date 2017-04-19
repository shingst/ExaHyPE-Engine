package eu.exahype;

import java.util.LinkedList;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;
import eu.exahype.node.PSolver;
import eu.exahype.io.FileSearch;

public class GenerateSolverRegistration extends DepthFirstAdapter {
  public Boolean valid = true;

  private java.io.BufferedWriter _writer;
  private java.io.StringWriter _methodBodyWriter;
  private java.io.StringWriter _versionBodyWriter;

  private DirectoryAndPathChecker _directoryAndPathChecker;

  private int _kernelNumber;
  private int _plotterNumber;
  private int _couplingNumber;
  
  private boolean _enableProfiler = false;
  
  private String _solverName;
  private String _solverType;
  
  private String _projectName;
  
  private boolean _useOptimisedKernels = false; //at least one solver uses optimised kernels

  private boolean _inALimitingADERDGSolver;

  public GenerateSolverRegistration(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
    _kernelNumber            = 0;
    _couplingNumber          = 0;
  }

  /// Write a single line to the toolkit registration information
  /// Example: writeVersionString("foo bar baz");
  public void writeVersionString(String singleline) {
      _versionBodyWriter.write("\n  ostream << \""+singleline+"\\n\";");
  }
  
  /// Write a key:value pair to the toolkit registration information
  /// Example: writeVersionString("foo", false);
  public void writeVersionString(String key, Object value) {
      _versionBodyWriter.write("\n  ostream << \""+key+": "+String.valueOf(value)+"\\n\";");
  }
  
  /// Write any C++ code in the key:value scheme at the value side.
  /// Example: writeVersionCode("foo", "baz::bar();");
  public void writeVersionCode(String key, String code) {
      _versionBodyWriter.write("\n  ostream << \""+key+": \";");
      _versionBodyWriter.write("\n  "+code);
      _versionBodyWriter.write("\n  ostream << \"\\n\";");
  }
  
  @Override
  public void inAProfiling(AProfiling node) {
    _enableProfiler = !node.getProfiler().getText().equals("NoOpProfiler");
    writeVersionString("enableProfiler", _enableProfiler);
  }
  
  @Override
  public void inAProject(AProject node) {
    _projectName = node.getName().getText();
    LinkedList<PSolver> solvers = node.getSolver();
    for(PSolver psolver : solvers) {
      if(psolver instanceof AAderdgSolver) {
        AAderdgSolver asolver = (AAderdgSolver) psolver;
        _useOptimisedKernels = _useOptimisedKernels || (asolver.getLanguage().getText().equals("C") 
                                  && (asolver.getKernel().getText().equals( eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier )
                                      ||  asolver.getKernel().getText().equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )));
      }
    }

    try {
      java.io.File logFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/KernelCalls.cpp");

      _writer = new java.io.BufferedWriter(new java.io.FileWriter(logFile));
      _methodBodyWriter = new java.io.StringWriter();
      _versionBodyWriter = new java.io.StringWriter();
      
      writeVersionString("projectName", _projectName);

      _writer.write("// This file is generated by the ExaHyPE toolkit.\n");
      _writer.write("// Please do not modify - it will be overwritten by the next\n");
      _writer.write("// ExaHyPE toolkit call.\n");
      _writer.write("// \n");
      _writer.write("// ========================\n");
      _writer.write("//   www.exahype.eu\n");
      _writer.write("// ========================\n\n");
      _writer.write("#include <sstream>\n");
      _writer.write("#include <ostream>\n");

      _writer.write("#include \"exahype/plotters/Plotter.h\"\n");
      _writer.write("#include \"exahype/profilers/ProfilerFactory.h\"\n");
      _writer.write("#include \"exahype/solvers/Solver.h\"\n");
      _writer.write("#include \"exahype/solvers/SolverCoupling.h\"\n\n");
      _writer.write("#include \"kernels/KernelCalls.h\"\n\n");

      _writer.write("#include \"kernels/GaussLegendreQuadrature.h\"\n");
      _writer.write("#include \"kernels/GaussLobattoQuadrature.h\"\n");
      _writer.write("#include \"kernels/LimiterProjectionMatrices.h\"\n");
      _writer.write("#include \"kernels/DGMatrices.h\"\n");
      _writer.write("#include \"kernels/DGBasisFunctions.h\"\n");
      _writer.write("#include \"buildinfo.h\"\n\n");
      if(_useOptimisedKernels) {
        _writer.write("#include \"kernels/aderdg/optimised/GaussLegendreQuadrature.h\"\n");
        _writer.write("#include \"kernels/aderdg/optimised/DGMatrices.h\"\n");
        writeVersionString("useOptimisedKernels", "YES");
      } else {
        writeVersionString("useOptimisedKernels", "no");
      }

      _methodBodyWriter.write("void kernels::initSolvers(exahype::Parser& parser, std::vector<std::string>& cmdlineargs) {\n");
      if (node.getSolver().size() == 0) {
        System.out.println("no solvers specified - create empty kernel calls ... ok");
      }
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
  

  private void writeProfilerCreation() {
      _methodBodyWriter.write("  std::string profiler_identifier = parser.getProfilerIdentifier();\n");
      _methodBodyWriter.write("  std::string metrics_identifier_list = parser.getMetricsIdentifierList();\n");
      _methodBodyWriter.write("  std::string profiling_output = parser.getProfilingOutputFilename();\n\n");
/*
      _methodBodyWriter.write(
          "  assertion1(metrics_identifier_list.find_first_of(\"{\") == 0,\n"+
          "           metrics_identifier_list);\n");

      _methodBodyWriter.write(
          "  assertion1(metrics_identifier_list.find_last_of(\"}\") ==\n"+
            "                 metrics_identifier_list.size() - 1,\n"+
          "             metrics_identifier_list);\n\n");
*/
      _methodBodyWriter.write("  // Split \"metric1,metric2...\" into \"metric1\", \"metric2\", ...\n");
      _methodBodyWriter.write("  std::vector<std::string> metrics_vector;\n");
      _methodBodyWriter.write("  std::stringstream sstream;\n");
      _methodBodyWriter.write("  sstream << metrics_identifier_list.substr(0, metrics_identifier_list.size() );\n");
      _methodBodyWriter.write("  std::string metric;\n");
      _methodBodyWriter.write(
          "  while (std::getline(sstream, metric, ',')) {\n"+
          "    metrics_vector.emplace_back(std::move(metric));\n"+
          "  }\n\n");

      _methodBodyWriter.write("  // Create profiler\n");
      _methodBodyWriter.write(
          "  auto profiler = exahype::profilers::ProfilerFactory::getInstance().create(\n"+
          "    profiler_identifier, metrics_vector, profiling_output);\n\n");
  }

  @Override
  public void inAAderdgSolver(eu.exahype.node.AAderdgSolver node) {
    try {
      _inALimitingADERDGSolver = false;
      
      _solverName = node.getName().getText();
      _solverType = "ADERDG";

      _writer.write("#include \"" + _solverName + ".h\"\n");

      _methodBodyWriter.write("  {\n");
      
      if (_enableProfiler) { writeProfilerCreation(); }
      
      _methodBodyWriter.write("  // Create and register solver\n");
      _methodBodyWriter.write("  exahype::solvers::RegisteredSolvers.push_back( new " + _projectName +
                          "::" + _solverName + "(parser.getMaximumMeshSize("+_kernelNumber+"), parser.getMaximumMeshDepth("+_kernelNumber+"), parser.getTimeStepping("+_kernelNumber+"), cmdlineargs"+
                           (_enableProfiler ? ", std::move(profiler)": ""));
      if (node.getConstants()!=null) {
          _methodBodyWriter.write( "  , parser.getParserView(" +  _kernelNumber + ")\n");
      }
      _methodBodyWriter.write( "  ));\n");
      _methodBodyWriter.write("  parser.checkSolverConsistency("+_kernelNumber+");\n\n");
      _methodBodyWriter.write("  \n");
  
      _methodBodyWriter.write("  }\n");
      
      _kernelNumber++;
      _plotterNumber = 0;

      writeVersionString("Kernel["+_kernelNumber+"].registration", "AderdgSolver");
      writeVersionString("Kernel["+_kernelNumber+"].type", _projectName+"::"+_solverName);
      writeVersionCode  ("Kernel["+_kernelNumber+"].parent", _projectName+"::Abstract"+_solverName+"::constantsToString(ostream);");
      writeVersionString("Kernel["+_kernelNumber+"].hasConstants", node.getConstants() != null);

      System.out.println("added creation of solver " + _solverName + " ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  };

  @Override
  public void inAFiniteVolumesSolver(eu.exahype.node.AFiniteVolumesSolver node) {
    try {
      _inALimitingADERDGSolver = false;
      
      _solverName = node.getName().getText();
      _solverType = "FiniteVolumes";
      
      int patchSize     = Integer.parseInt(node.getPatchSize().getText());
      
      _writer.write("#include \"" + _solverName + ".h\"\n");

      _methodBodyWriter.write("  {\n"); 
      
      if (_enableProfiler) { writeProfilerCreation(); }

      _methodBodyWriter.write("  // Create and register solver\n");
      _methodBodyWriter.write("  exahype::solvers::RegisteredSolvers.push_back( new " + _projectName +
                          "::" + _solverName + "(parser.getMaximumMeshSize("+_kernelNumber+"), parser.getMaximumMeshDepth("+_kernelNumber+"), parser.getTimeStepping("+_kernelNumber+")"+
                          (_enableProfiler ? ", std::move(profiler)": ""));
      if (node.getConstants()!=null) {
        _methodBodyWriter.write( "  , parser.getParserView(" +  _kernelNumber + ")\n");
      }
      _methodBodyWriter.write( "  , cmdlineargs));\n");
      _methodBodyWriter.write("  parser.checkSolverConsistency("+_kernelNumber+");\n\n");
      _methodBodyWriter.write("  \n");
      
      _methodBodyWriter.write("  }\n");
      
      _kernelNumber++;
      _plotterNumber = 0;
      
      writeVersionString("Kernel["+_kernelNumber+"].registration", "FiniteVolumesSolver");
      writeVersionString("Kernel["+_kernelNumber+"].type", _projectName+"::"+_solverName);
      writeVersionCode  ("Kernel["+_kernelNumber+"].parent", _projectName+"::Abstract"+_solverName+"::constantsToString(ostream);");
      writeVersionString("Kernel["+_kernelNumber+"].hasConstants", node.getConstants());
      writeVersionString("Kernel["+_kernelNumber+"].patchSize", patchSize);

      System.out.println("added creation of solver " + _solverName + " ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }

  @Override
  public void inALimitingAderdgSolver(ALimitingAderdgSolver node) {
    try {
      _inALimitingADERDGSolver = true;
      
      _solverName  = node.getName().getText();
      _solverType  = "ADERDG"; // TODO(Dominic): We currently use the ADERDG degrees of freedom as the plotter input for the LimitingADERDGSolver

      _writer.write("#include \"exahype/solvers/LimitingADERDGSolver.h\"\n");
      
      // Consistency: We have the same definition at CreateSolverClasses.inALimitingAderdgSolver()
      String solverNameADERDG = _solverName+"_ADERDG";
      String solverNameFV     = _solverName+"_FV";
      
      _writer.write("#include \"" + solverNameADERDG + ".h\"\n");
      _writer.write("#include \"" + solverNameFV + ".h\"\n");

      _methodBodyWriter.write("  {\n\n");
      _methodBodyWriter.write("  // Create and register solver\n");
      
      _methodBodyWriter.write("  exahype::solvers::Solver* solver = nullptr;\n\n");

      writeVersionString("Kernel["+_kernelNumber+"].registration", "LimitingAderdgSolver");
      //writeVersionString("Kernel["+_kernelNumber+"].type", "exahype::solvers::LimitingADERDGSolver");
      writeVersionString("Kernel["+_kernelNumber+"].type[FV]", _projectName + "::" + solverNameFV);
      writeVersionString("Kernel["+_kernelNumber+"].type[ADERDG]", _projectName + "::" + solverNameADERDG);
      writeVersionCode  ("Kernel["+_kernelNumber+"].abstract[FV]", _projectName+"::Abstract"+solverNameFV+"::constantsToString(ostream);");
      writeVersionCode  ("Kernel["+_kernelNumber+"].abstract[ADERDG]", _projectName+"::Abstract"+solverNameADERDG+"::constantsToString(ostream);");
      writeVersionString("Kernel["+_kernelNumber+"].hasConstants", node.getConstants() != null);
      
      // ADER-DG
      _methodBodyWriter.write("  {\n");
      
      if (_enableProfiler) { writeProfilerCreation(); }
      
      _methodBodyWriter.write("  solver = new " + _projectName +
                          "::" + solverNameADERDG+"(parser.getMaximumMeshSize("+_kernelNumber+"), parser.getMaximumMeshDepth("+_kernelNumber+"), parser.getTimeStepping("+_kernelNumber+")"+
                           (_enableProfiler ? ", std::move(profiler)": "")+
                           ", cmdlineargs");
      if (node.getConstants()!=null) {
          _methodBodyWriter.write( "  , parser.getParserView(" +  _kernelNumber + ")\n");
        }
      _methodBodyWriter.write( ");\n");
  
      _methodBodyWriter.write("  }\n");
      
      _methodBodyWriter.write("  std::unique_ptr<exahype::solvers::ADERDGSolver> aderdgSolver(static_cast<exahype::solvers::ADERDGSolver*>(solver));\n");
      
      // FV
      _methodBodyWriter.write("  {\n");
      
      if (_enableProfiler) { writeProfilerCreation(); }

      _methodBodyWriter.write("  solver = new " + _projectName +
                          "::" + solverNameFV+"(parser.getMaximumMeshSize("+_kernelNumber+"), parser.getTimeStepping("+_kernelNumber+")"+
                          (_enableProfiler ? ", std::move(profiler)": "")+
                          ",cmdlineargs");
      if (node.getConstants()!=null) {
        _methodBodyWriter.write( "  , parser.getParserView(" +  _kernelNumber + ")\n");
      }
      _methodBodyWriter.write( ");\n");
      _methodBodyWriter.write("  }\n");
      
      _methodBodyWriter.write("  std::unique_ptr<exahype::solvers::FiniteVolumesSolver> finiteVolumesSolver(static_cast<exahype::solvers::FiniteVolumesSolver*>(solver));\n");
      
      // Limiting ADER-DG
      _methodBodyWriter.write("  \n");
      _methodBodyWriter.write("  exahype::solvers::RegisteredSolvers.push_back(\n"
          + "    new exahype::solvers::LimitingADERDGSolver(\""+_solverName+"\",std::move(aderdgSolver),std::move(finiteVolumesSolver),parser.getDMPRelaxationParameter(0),parser.getDMPDifferenceScaling(0)) );\n");
      
      _methodBodyWriter.write("  parser.checkSolverConsistency("+_kernelNumber+");\n");
      _methodBodyWriter.write("  }\n\n");
      _methodBodyWriter.write("  \n");
      
      _kernelNumber++;
      _plotterNumber = 0;
      
      System.out.println("added creation of solver " + _solverName + " ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
  
  @Override
  public void inAPlotSolution(eu.exahype.node.APlotSolution node) {
    try {
      String plotterName = node.getName().getText();
      String plotterType = node.getPlotterType().getText();

      _writer.write(FileSearch.PPinclude(plotterName, _directoryAndPathChecker.outputDirectory.getAbsolutePath()));
      
      writeVersionString("Kernel["+(_kernelNumber-1)+"].Plotter["+_plotterNumber+"]", _projectName + "::" + plotterName);

      if (plotterType.equals("user::defined")) {
        _methodBodyWriter.write(
            "  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter("
                + (_kernelNumber - 1) + "," + _plotterNumber + ",parser, new " + _projectName + "::" + plotterName + "()) );\n\n");
      } else {
        if (_inALimitingADERDGSolver) {
          _methodBodyWriter.write(
              "  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter("
                  + (_kernelNumber - 1) + "," + _plotterNumber + ",parser,new " + _projectName + "::" + plotterName + "(  *static_cast<exahype::solvers::LimitingADERDGSolver*>(exahype::solvers::RegisteredSolvers[" + (_kernelNumber-1) + "])) ));\n\n");
        } else {
          _methodBodyWriter.write(
              "  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter("
                  + (_kernelNumber - 1) + "," + _plotterNumber + ",parser,new " + _projectName + "::" + plotterName + "(  *static_cast<" + _projectName + "::" + _solverName + "*>(exahype::solvers::RegisteredSolvers[" + (_kernelNumber-1) + "])) ));\n\n");
        }
      }
      _plotterNumber++;
      System.out.println("added plotter ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }

  @Override
  public void outAProject(AProject node) {    
    try {
      _methodBodyWriter.write("\n");
      _methodBodyWriter.write(
          "  std::set<int> orders;\n"+
          "  for (const auto p : exahype::solvers::RegisteredSolvers) {\n"+
          "    orders.insert(p->getNodesPerCoordinateAxis()-1);\n"+
          "  }\n"+
          "  kernels::initGaussLegendreNodesAndWeights(orders);\n"+
          "  kernels::initGaussLobattoNodesAndWeights(orders);\n"+
          "  kernels::initLimiterProjectionMatrices(orders);\n"+
          "  kernels::initDGMatrices(orders);\n" +
          "  kernels::initBasisFunctions(orders);\n");
      if(_useOptimisedKernels) {
        _methodBodyWriter.write("  kernels::aderdg::optimised::initGaussLegendreNodesAndWeights(orders);\n");
        _methodBodyWriter.write("  kernels::aderdg::optimised::initDGMatrices(orders);\n");
      }
      _methodBodyWriter.write("}\n"); // close initSolvers(...)
      _methodBodyWriter.write("\n");
      _methodBodyWriter.write("\n");
      _methodBodyWriter.write("void kernels::finalise() {\n");
      _methodBodyWriter.write(
          "  std::set<int> orders;\n"+
          "  for (const auto p : exahype::solvers::RegisteredSolvers) {\n" +
          "    orders.insert(p->getNodesPerCoordinateAxis()-1);\n" +
          "  }\n" +
          "  kernels::freeGaussLegendreNodesAndWeights(orders);\n"+
          "  kernels::freeGaussLobattoNodesAndWeights(orders);\n"+
          "  kernels::freeLimiterProjectionMatrices(orders);\n"+
          "  kernels::freeDGMatrices(orders);\n"+
          "  kernels::freeBasisFunctions(orders);\n");
      if(_useOptimisedKernels) {
        _methodBodyWriter.write("  kernels::aderdg::optimised::freeGaussLegendreNodesAndWeights(orders);\n");
        _methodBodyWriter.write("  kernels::aderdg::optimised::freeDGMatrices(orders);\n");
      }    
      _methodBodyWriter.write(
          "\n"+
          "  for (auto solver : exahype::solvers::RegisteredSolvers) {\n"+
          "    delete solver;\n"+
          "  }\n"+
          "  exahype::solvers::RegisteredSolvers.clear();\n\n"+
          "  for (auto plotter : exahype::plotters::RegisteredPlotters) {\n"+
          "    delete plotter;\n"+
          "  }\n"+
          "  exahype::plotters::RegisteredPlotters.clear();\n"+
          "  for (auto coupling : exahype::solvers::RegisteredSolverCouplings) {\n"+
          "    delete coupling;\n"+
          "  }\n"+
          "  exahype::solvers::RegisteredSolverCouplings.clear();\n");
      _methodBodyWriter.write("}\n");
      _writer.write("\n");
      _writer.write("\n");
      _writer.write("\n");
      _writer.write(_methodBodyWriter.toString());
      _writer.write("\n");
      _writer.write("\n");
      _writer.write("\n");
      _writer.write("void kernels::toString(std::ostream& ostream) {\n");
      _writer.write("/* Generated SolverRegistration code by the toolkit */\n");
      _writer.write(_versionBodyWriter.toString());
      _writer.write("\n}");
      _writer.write("\n");

      System.out.println("configured all solver solvers ... ok");

      _writer.write("\n\n");
      _writer.close();
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
}
