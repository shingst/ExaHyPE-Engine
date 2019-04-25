/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"
#include "peano/parallel/JoinDataBufferPool.h"


#include "exahype/main.h"
#include "exahype/parser/Parser.h"
#include "exahype/Vertex.h"
#include "exahype/runners/Runner.h"

#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include <vector>
#include <cstdlib> // getenv, exit
#include <iostream>
#include <cstdio>

#include <exahype/main.h>

tarch::logging::Log _log("exahype");

namespace muq{

      bool firstRun=true;
      exahype::parser::Parser parser;
      std::vector<std::string> cmdlineargs;

int init(int argc, char** argv) {
  firstRun=true;
  //
  //   Parse config file
  // =====================
  //
  std::string progname = argv[0];
  
  if (argc < 2) {
    logError("main()", "Usage: " << progname << " --help");
    return -1;
  }

  // cmdlineargs contains all argv expect the progname.
  std::vector<std::string> cmd_copy(argv + 1, argv + argc);
  cmdlineargs = cmd_copy;
  std::string firstarg = cmdlineargs[0];

  bool showHelp    = firstarg == "-h" || firstarg == "--help";
  bool showVersion = firstarg == "-v" || firstarg == "--version";
  bool runTests    = firstarg == "-t" || firstarg == "--tests";
  bool runPingPong = firstarg == "-p" || firstarg == "--pingpong";
  bool showCompiledSpecfile = firstarg == "--show-specfile";
  bool runCompiledSpecfile  = firstarg == "--built-in-specfile";

  //
  //   Early standalone options
  //   ========================
  //

  if(showHelp) {
      exahype::help(progname);
      return EXIT_SUCCESS;
  }

  if(showVersion) {
    std::cout << exahype::version(progname);
    return EXIT_SUCCESS;
  }
  
  if(showCompiledSpecfile) {
    // Unfortunately, we cannot avoid here to get the output dirtied by the
    // tarch::parallel::Node<static>::reserveFreeTag() log outputs.
    // The only alternative to get the clean specfile would be to dump it to
    // a file.
    
    // if this line does not compile for you, rebuild and rerun the toolkit.
    std::cout << std::string(kernels::compiledSpecfile());
    return EXIT_SUCCESS;
  }

  //
  //   Setup environment
  //   =================
  //
  peano::fillLookupTables();
  
  int parallelSetup = peano::initParallelEnvironment(&argc, &argv);
  if (parallelSetup != 0) {
#ifdef Parallel
    // Please do not use the logging if MPI doesn't work properly.
    std::cerr << "mpi initialisation wasn't successful. Application shut down"
              << std::endl;
#else
    _log.error("main()",
               "mpi initialisation wasn't successful. Application shut down");
#endif
    return parallelSetup;
 }

  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
  if (sharedMemorySetup != 0) {
    logError("main()",
             "shared memory initialisation wasn't successful. Application shut "
             "down");
    return sharedMemorySetup;
  }

  if (runPingPong) {
    return exahype::pingPongTest();
  }

  if (runTests) {
    //
    //   Run tests
    // =============
    // Our unit tests do cover the generic ADER-DG kernels. The generic kernels do
    // parallelise. As a consequence, they connect to the autotuning feature.
    // Autotuning however is not set up yet, so this will fail. We therefore
    // disable the unit tests in shared memory mode.
    //

    //#if (defined(Debug) || defined(Asserts)) && !defined(SharedMemoryParallelisation)
    //if(! std::getenv("EXAHYPE_SKIP_TESTS")) { // cf issue #74
    tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
    int testExitCode = tarch::tests::TestCaseRegistry::getInstance()
                           .getTestCaseCollection()
                           .getNumberOfErrors();

    if (testExitCode != 0) {
      logError("main()", "unit tests failed. Quit.");
      return -2;
    }
    else {
      return EXIT_SUCCESS;
    }
  }

  //
  //   Parse specification file
  // =====================================
  //


  std::stringstream specfile;
  std::string specFileName;
  if(runCompiledSpecfile) {
    specFileName = "builtin";
    specfile.str(std::string(kernels::compiledSpecfile()));
  } else {
    specFileName = firstarg;
    specfile.str(kernels::readSpecificationFileToJSON(specFileName));
  }
  parser.readFile(specfile, specFileName);

  if (!parser.isValid()) {
    logError("main()", "invalid config file. Quit");
    return -2;
  }

  //
  //   Init solver registries
  // =====================================
  //
  kernels::registerSolvers(parser);

  //
  //   Configure the logging
  // =========================
  //
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  #if defined(Parallel) || defined(PerformanceAnalysis)
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      true,   // logMachineName
      true,   // logMessageType
      true,   // logTrace
      parser.getLogFileName() );
  #elif defined(Asserts) || defined(Debug)
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      true,   // logTrace
      parser.getLogFileName() );
  #else
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      false,   // logTrace
      parser.getLogFileName() );
  #endif

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::LogFilterFileReader::parsePlainTextFile( "exahype.log-filter" );

//    tarch::logging::CommandLineLogger::getInstance().clearFilterList();
/*
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("info", false));
    #if !defined(Asserts)
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry(
            "info", -1, "peano::grid", true));
    #endif
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", true));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", -1,
                                                             "exahype", false));
*/
  return 0;
  }

  int run_exahype(){
	if(!firstRun){
		tarch::parallel::NodePool::getInstance().shutdown();
		tarch::parallel::NodePool::getInstance().init();
	}
	exahype::runners::Runner runner(parser, cmdlineargs); 
	runner.initHeaps();
	std::cout << "starting first run" << std::endl;
	int programExitCode = runner.run();
	assert(programExitCode == 0);
        firstRun = false;
	std::cout << "done run 1" << std::endl;
	runner.shutdownHeaps();
        peano::parallel::JoinDataBufferPool::getInstance().releaseMessages();
	//std::cout << "done run 2" << std::endl;
	//tarch::parallel::NodePool::getInstance().terminate();
	//std::cout << "done run 3" << std::endl;
	//tarch::parallel::NodePool::getInstance().restart();
	//std::cout << "done run 4" << std::endl;
	return programExitCode;
}

int finalize(){
  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();
  kernels::finalise();
  return 0;//programExitCode;
}

bool setCommunicator(MPI_Comm communicator){
  return tarch::parallel::Node::getInstance().setCommunicator(communicator);
}



}
