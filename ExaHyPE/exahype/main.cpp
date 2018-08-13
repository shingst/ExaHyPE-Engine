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

tarch::logging::Log _log("");

/**
 * The ping pong test has to be triggered by main very very early. There should
 * be no other message in the MPI subsystem.
 */
int exahype::pingPongTest() {
  bool correct = true;
  #if defined(Parallel) && defined(Asserts)
  logInfo( "run()", "start ping pong test .... if test fails, please retranslate with -DnoPackedRecords" );
  exahype::Vertex::initDatatype();
  exahype::Vertex sendVertex[5];

  if (tarch::parallel::Node::getInstance().getNumberOfNodes()>1) {
    if (tarch::parallel::Node::getInstance().getRank()==0) {
      sendVertex[0].setPosition( tarch::la::Vector<DIMENSIONS,double>(2.0), 4);
      sendVertex[0].setAdjacentRank( 0, 10 );
      sendVertex[0].setAdjacentRank( 1, 11 );
      sendVertex[0].setAdjacentRank( 2, 12 );
      sendVertex[0].setAdjacentRank( 3, 13 );
      sendVertex[1].setPosition( tarch::la::Vector<DIMENSIONS,double>(3.0), 5);
      sendVertex[1].setAdjacentRank( 0, 20 );
      sendVertex[1].setAdjacentRank( 1, 21 );
      sendVertex[1].setAdjacentRank( 2, 22 );
      sendVertex[1].setAdjacentRank( 3, 23 );
      sendVertex[2].setPosition( tarch::la::Vector<DIMENSIONS,double>(4.0), 6);
      sendVertex[2].setAdjacentRank( 0, 30 );
      sendVertex[2].setAdjacentRank( 1, 31 );
      sendVertex[2].setAdjacentRank( 2, 32 );
      sendVertex[2].setAdjacentRank( 3, 33 );

      sendVertex[0].send(1,100,false,Vertex::MPIDatatypeContainer::ExchangeMode::Blocking);
      sendVertex[0].send(1,100,true,Vertex::MPIDatatypeContainer::ExchangeMode::Blocking);
      MPI_Send( sendVertex, 3, exahype::Vertex::MPIDatatypeContainer::Datatype, 1, 100, tarch::parallel::Node::getInstance().getCommunicator() );
    }
    if (tarch::parallel::Node::getInstance().getRank()==1) {
      exahype::Vertex receivedVertex;

      receivedVertex.receive(0,100,false,Vertex::MPIDatatypeContainer::ExchangeMode::Blocking);
      assertion1( receivedVertex.getLevel()==4, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(0)==2.0, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(1)==2.0, receivedVertex.toString() );
      correct &=  receivedVertex.getLevel()==4 && receivedVertex.getX()(0)==2.0 && receivedVertex.getX()(1)==2.0;
      #ifdef Dim3
      assertion1( receivedVertex.getX()(2)==2.0, receivedVertex.toString() );
      correct &=  receivedVertex.getX()(2)==2.0;
      #endif


      receivedVertex.receive(0,100,true,Vertex::MPIDatatypeContainer::ExchangeMode::Blocking);
      assertion1( receivedVertex.getLevel()==4, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(0)==2.0, receivedVertex.toString() );
      assertion1( receivedVertex.getX()(1)==2.0, receivedVertex.toString() );
      correct &=  receivedVertex.getLevel()==4 && receivedVertex.getX()(0)==2.0 && receivedVertex.getX()(1)==2.0;
      #ifdef Dim3
      assertion1( receivedVertex.getX()(2)==2.0, receivedVertex.toString() );
      correct &=  receivedVertex.getX()(2)==2.0;
      #endif

      exahype::Vertex receivedVertices[5];
      MPI_Status status;
      MPI_Recv( receivedVertices, 3, exahype::Vertex::MPIDatatypeContainer::Datatype, 0, 100, tarch::parallel::Node::getInstance().getCommunicator(), &status );
      assertion3( receivedVertices[0].getLevel()==4,  receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[0].getX()(0)==2.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[0].getX()(1)==2.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      correct &=  receivedVertices[0].getLevel()==4 && receivedVertices[0].getX()(0)==2.0 && receivedVertices[0].getX()(1)==2.0;
      #ifdef Dim3
      assertion3( receivedVertices[0].getX()(2)==2.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      correct &=  receivedVertices[0].getX()(2)==2.0;
      #endif

      assertion3( receivedVertices[1].getLevel()==5,  receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[1].getX()(0)==3.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[1].getX()(1)==3.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      correct &=  receivedVertices[1].getLevel()==5 && receivedVertices[1].getX()(0)==3.0 && receivedVertices[1].getX()(1)==3.0;
      #ifdef Dim3
      assertion3( receivedVertices[1].getX()(2)==3.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      correct &=  receivedVertices[1].getX()(2)==3.0;
      #endif

      assertion3( receivedVertices[2].getLevel()==6,  receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[2].getX()(0)==4.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      assertion3( receivedVertices[2].getX()(1)==4.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      correct &=  receivedVertices[2].getLevel()==6 && receivedVertices[2].getX()(0)==4.0 && receivedVertices[2].getX()(1)==4.0;
      #ifdef Dim3
      assertion3( receivedVertices[2].getX()(2)==4.0, receivedVertices[0].toString(), receivedVertices[1].toString(), receivedVertices[2].toString() );
      correct &=  receivedVertices[2].getX()(2)==4.0;
      #endif
    }
    MPI_Barrier( tarch::parallel::Node::getInstance().getCommunicator() );

    if (tarch::parallel::Node::getInstance().getRank()==0) {
      //exahype::Vertex* heapVertices = new exahype::Vertex[5];
      exahype::Vertex heapVertices[5];
      MPI_Recv( heapVertices, 3, exahype::Vertex::MPIDatatypeContainer::Datatype, 1, 1, tarch::parallel::Node::getInstance().getCommunicator(), MPI_STATUS_IGNORE );

      assertion3( heapVertices[0].getLevel()==4,  heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      assertion3( heapVertices[0].getX()(0)==2.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      assertion3( heapVertices[0].getX()(1)==2.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      correct &=  heapVertices[0].getLevel()==4 && heapVertices[0].getX()(0)==2.0 && heapVertices[0].getX()(1)==2.0;
      #ifdef Dim3
      assertion3( heapVertices[0].getX()(2)==2.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      correct &=  heapVertices[0].getX()(2)==2.0;
      #endif

      assertion3( heapVertices[1].getLevel()==5,  heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      assertion3( heapVertices[1].getX()(0)==3.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      assertion3( heapVertices[1].getX()(1)==3.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      correct &=  heapVertices[1].getLevel()==5 && heapVertices[1].getX()(0)==3.0 && heapVertices[1].getX()(1)==3.0;
      #ifdef Dim3
      assertion3( heapVertices[1].getX()(2)==3.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      correct &=  heapVertices[1].getX()(2)==3.0;
      #endif

      assertion3( heapVertices[2].getLevel()==6,  heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      assertion3( heapVertices[2].getX()(0)==4.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      assertion3( heapVertices[2].getX()(1)==4.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      correct &=  heapVertices[2].getLevel()==6 && heapVertices[2].getX()(0)==4.0 && heapVertices[2].getX()(1)==4.0;
      #ifdef Dim3
      assertion3( heapVertices[2].getX()(2)==4.0, heapVertices[0].toString(), heapVertices[1].toString(), heapVertices[2].toString() );
      correct &=  heapVertices[2].getX()(2)==4.0;
      #endif
    }
    if (tarch::parallel::Node::getInstance().getRank()==1) {
      exahype::Vertex* heapVertices = new exahype::Vertex[5];
      heapVertices[0].setPosition( tarch::la::Vector<DIMENSIONS,double>(2.0), 4);
      heapVertices[1].setPosition( tarch::la::Vector<DIMENSIONS,double>(3.0), 5);
      heapVertices[2].setPosition( tarch::la::Vector<DIMENSIONS,double>(4.0), 6);
      MPI_Send( heapVertices, 3, exahype::Vertex::MPIDatatypeContainer::Datatype, 0, 1, tarch::parallel::Node::getInstance().getCommunicator() );
      delete[] heapVertices;
    }
    MPI_Barrier( tarch::parallel::Node::getInstance().getCommunicator() );
  }

  logInfo( "run()", " ping pong test ok" );
  #elif defined(Parallel)
  logInfo( "run()", "ping pong tests requires compile with -DAsserts as it uses assert data to validate that right content is exchanged" );
  #endif

  return correct ? EXIT_SUCCESS : -1;
}


int exahype::main(int argc, char** argv) {
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
  std::vector<std::string> cmdlineargs(argv + 1, argv + argc);
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
    help(progname);
    return EXIT_SUCCESS;
  }

  if(showVersion) {
    std::cout << version(progname);
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
    return pingPongTest();
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

  exahype::parser::Parser parser;

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
  if (!tarch::logging::LogFilterFileReader::parsePlainTextFile(
          "exahype.log-filter")) {
    tarch::logging::CommandLineLogger::getInstance().clearFilterList();
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
  }

  exahype::runners::Runner runner(parser, cmdlineargs);
  int programExitCode = runner.run();

  if (programExitCode == 0) {
#ifdef Parallel
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      logInfo("main()", "Peano terminates successfully");
    }
#else
    logInfo("main()", "Peano terminates successfully");
#endif
  } else {
    logInfo("main()", "quit with error code " << programExitCode);
  }

  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();

  kernels::finalise();

  return programExitCode;
}


void exahype::help(const std::string& programname) {
  std::cout << "Usage: " << programname << " [-hvt] <YourApplication.exahype>\n";
  std::cout << "\n";
  std::cout << "   where YourApplication.exahype is an ExaHyPE specification file.\n";
  std::cout << "   Note that you should have compiled ExaHyPE with this file as there\n";
  std::cout << "   are some compile time constants.\n";
  std::cout << "\n";
  std::cout << "   Other possible parameters:\n";
  std::cout << "\n";
  std::cout << "    --help    | -h       Show this help message\n";
  std::cout << "    --version | -v       Show version and other hard coded information\n";
  std::cout << "    --tests   | -t       Run the unit tests\n";
  std::cout << "    --show-specfile      Show the specification file the binary was built with\n";
  std::cout << "    --built-in-specfile  Run with the spec. file the binary was built with\n";
  std::cout << "\n";
}


#ifndef EXAHYPE_LATE_TAKEOVER

/**
 * By default, ExaHyPE provides the main function entrance of the program.
 * If you want to embed ExaHyPE however as an engine in some other program,
 * you can call the exahype::main(argc,argv) at any later time.
 *
 * Thus you can treat ExaHyPE similar to a GUI toolkit or game engine even loop.
 * To do so, just define EXAHYPE_LATE_TAKEOVER. Don't forget to start ExaHyPE
 * finally with calling exahype::main on yourself.
 **/
int main(int argc, char**argv) {
	exahype::main(argc, argv);
}

#endif
