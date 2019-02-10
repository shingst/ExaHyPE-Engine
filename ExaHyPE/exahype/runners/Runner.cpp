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

#include "exahype/runners/Runner.h"

#include <cmath>

#include "../../../Peano/mpibalancing/HotspotBalancing.h"

#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"
#include "exahype/mappings/LoadBalancing.h"

#include "tarch/Assertions.h"

#include "tarch/logging/CommandLineLogger.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#include "tarch/multicore/Core.h"
#include "tarch/multicore/MulticoreDefinitions.h"
#include "tarch/multicore/Jobs.h"

#ifdef SharedTBB
#include "tarch/multicore/tbb/Jobs.h"
#endif

#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"

#include "peano/geometry/Hexahedron.h"

#include "peano/utils/UserInterface.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"

#include "exahype/parser/ParserView.h"

#ifdef Parallel
#include "mpibalancing/GreedyBalancing.h"
#include "mpibalancing/FairNodePoolStrategy.h"
#include "mpibalancing/SFCDiffusionNodePoolStrategy.h"
#endif
#include "exahype/plotters/Plotter.h"

#include "exahype/mappings/MeshRefinement.h"
#include "exahype/mappings/RefinementStatusSpreading.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/performanceanalysis/SpeedupLaws.h"

#include "peano/datatraversal/TaskSet.h"

#ifdef TBBInvade
#include "shminvade/SHMController.h"
#include "shminvade/SHMSharedMemoryBetweenTasks.h"
#include "shminvade/SHMOccupyAllCoresStrategy.h"
#include "shminvade/SHMMultipleRanksPerNodeStrategy.h"
#endif


tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner(exahype::parser::Parser& parser, std::vector<std::string>& cmdlineargs) :
    _parser(parser),
    _cmdlineargs(cmdlineargs),
    _boundingBoxSize(0.0),
    _meshRefinements(0),
    _localRecomputations(0),
    _predictorReruns(0) {
  #ifdef TBBInvade
  _shmInvade = nullptr;
  #endif
}

exahype::runners::Runner::~Runner() {}

void exahype::runners::Runner::initDistributedMemoryConfiguration() {
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    //
    // Configure answering behaviour of global node pool
    // =================================================
    //
    if ( _parser.compareNodePoolStrategy( "FCFS" ) ) {
      tarch::parallel::NodePool::getInstance().setStrategy(
        new tarch::parallel::FCFSNodePoolStrategy()
      );
      logInfo("initDistributedMemoryConfiguration()", "load balancing relies on FCFS answering strategy");
    }
    else if (_parser.compareNodePoolStrategy( "fair" )) {
      int ranksPerNode = _parser.getRanksPerNode();
      if (ranksPerNode<=0) {
        logError( "initDistributedMemoryConfiguration()", "please inform fair balancing how many ranks per node you use through value \"" << _parser.getRanksPerNode() << ":XXX\". Read value " << ranksPerNode << " is invalid" );
        ranksPerNode = 1;
      }
      if ( ranksPerNode>tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
        logWarning( "initDistributedMemoryConfiguration()", "value \"" << _parser.getRanksPerNode() << ":XXX\" exceeds total rank count. Reset to 1" );
        ranksPerNode = 1;
      }
      tarch::parallel::NodePool::getInstance().setStrategy(
        new mpibalancing::FairNodePoolStrategy(ranksPerNode,1,_parser.getNodePoolAnsweringTimeout())
      );
      logInfo("initDistributedMemoryConfiguration()", "load balancing relies on fair answering strategy with " << ranksPerNode << " rank(s) per node") ;
    }
    else if (_parser.compareNodePoolStrategy( "sfc-diffusion" )) {
      int ranksPerNode = _parser.getRanksPerNode();
      if (ranksPerNode<=0) {
        logError( "initDistributedMemoryConfiguration()", "please inform SFC balancing how many ranks per node you use through value \"RanksPerNode:XXX\". Read value " << ranksPerNode << " is invalid" );
        ranksPerNode = 1;
      }
      if ( ranksPerNode>tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
        logWarning( "initDistributedMemoryConfiguration()", "value \"" << _parser.getRanksPerNode() << ":XXX\" exceeds total rank count. Reset to 1" );
        ranksPerNode = 1;
      }
      if (tarch::parallel::Node::getInstance().getNumberOfNodes() % ranksPerNode != 0) {
        logError( "initDistributedMemoryConfiguration()", "Value of \"RanksPerNode:XXX\" does not fit to total number of ranks. ExaHyPE requires homogeneous rank distribution" );
        ranksPerNode = 1;
      }
      int primaryRanksPerNode = _parser.getIntFromPath("/distributed_memory/primary_ranks_per_node");
      if (primaryRanksPerNode<=0) {
        logError( "initDistributedMemoryConfiguration()", "please inform SFC balancing how many primary ranks per node you use through value \"primary-ranks-per-node:XXX\". Read value " << primaryRanksPerNode << " is invalid" );
        primaryRanksPerNode = 1;
      }
      if ( ranksPerNode<primaryRanksPerNode ) {
        logWarning( "initDistributedMemoryConfiguration()", "value " << _parser.getRanksPerNode() << " is smaller than primary-ranks-per-node. Reset to 1" );
        primaryRanksPerNode = 1;
      }
      tarch::parallel::NodePool::getInstance().setStrategy(
        new mpibalancing::SFCDiffusionNodePoolStrategy(ranksPerNode,primaryRanksPerNode,_parser.getNodePoolAnsweringTimeout())
      );
      logInfo("initDistributedMemoryConfiguration()",
        "load balancing relies on an sfc-diffusion answering strategy with " << ranksPerNode <<
        " rank(s) per node while " << primaryRanksPerNode << " rank(s) per node are primary ranks" );
    }
    else {
      logError("initDistributedMemoryConfiguration()", "no valid load balancing answering strategy specified");
      _parser.invalidate();
    }
  }

  // load balancing type
  switch (_parser.getMPILoadBalancingType()) {
    case parser::Parser::MPILoadBalancingType::Static:
      mappings::MeshRefinement::DynamicLoadBalancing = false;
      break;
    case parser::Parser::MPILoadBalancingType::Dynamic:
      mappings::MeshRefinement::DynamicLoadBalancing = true;
      break;
    default:
      logError("initDistributedMemoryConfiguration()", "no valid load balancing type specified");
      _parser.invalidate();
      break;
  }

  // load balancing strategy
  // basically a switch-case
  if ( _parser.compareMPILoadBalancingStrategy( "greedy_naive" )) {
    exahype::mappings::LoadBalancing::setLoadBalancingAnalysis( exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Greedy );
  }
  else if ( _parser.compareMPILoadBalancingStrategy( "greedy_regular" )) {
    exahype::mappings::LoadBalancing::setLoadBalancingAnalysis( exahype::mappings::LoadBalancing::LoadBalancingAnalysis::GreedyWithRegularityAnalysis );
  }
  else if ( _parser.compareMPILoadBalancingStrategy( "hotspot" )) {
    exahype::mappings::LoadBalancing::setLoadBalancingAnalysis( exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Hotspot );
  }
  else {
    logError("initDistributedMemoryConfiguration()", "no valid load balancing analysis type specified");
    _parser.invalidate();
  }

  //
  // Configure answering behaviour of global node pool
  // =================================================
  //
  if ( _parser.getMPILoadBalancingType()==exahype::parser::Parser::MPILoadBalancingType::Static ) {
    switch ( exahype::mappings::LoadBalancing::getLoadBalancingAnalysis() ) {
      case exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Greedy:
        logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins");
        peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
            new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(false)
        );
        break;
      case exahype::mappings::LoadBalancing::LoadBalancingAnalysis::GreedyWithRegularityAnalysis:
        logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins (mpibalancing/GreedyBalancing)");
        peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new mpibalancing::GreedyBalancing(
            getCoarsestGridLevelForLoadBalancing(_boundingBoxSize),
            getFinestUniformGridLevelForLoadBalancing(_boundingBoxSize) /*boundary regularity*/
          )
        );
        break;
      case exahype::mappings::LoadBalancing::LoadBalancingAnalysis::Hotspot:
        logInfo("initDistributedMemoryConfiguration()", "use global hotspot elimination without joins (mpibalancing/StaticBalancing)");
        peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
            new mpibalancing::HotspotBalancing(
                false,getFinestUniformGridLevelForLoadBalancing(_boundingBoxSize), /*boundary regularity*/
                tarch::parallel::Node::getInstance().getNumberOfNodes()/THREE_POWER_D
          )
        );
        break;
    }
  }
  // Dynamic load balancing
  // ----------------------
  // @todo Missing
  else {
    logError("initDistributedMemoryConfiguration()", "only MPI static load balancing supported so far. ");
    _parser.invalidate();
  }

  if ( _parser.isValid() ) {
    tarch::parallel::NodePool::getInstance().restart();

    tarch::parallel::Node::getInstance().setDeadlockTimeOut(_parser.getMPITimeOut());
    tarch::parallel::Node::getInstance().setTimeOutWarning(_parser.getMPITimeOut()/2);
    logInfo("initDistributedMemoryConfiguration()", "use MPI time out of " << _parser.getMPITimeOut() << " (warn after half the timeout span)");

    const int bufferSize = _parser.getMPIBufferSize();
    peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(bufferSize);
    peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(bufferSize);
    logInfo("initDistributedMemoryConfiguration()", "use MPI buffer size of " << bufferSize);

    if ( _parser.getSkipReductionInBatchedTimeSteps() ) {
      logInfo("initDistributedMemoryConfiguration()", "allow ranks to skip reduction and broadcasts within a batch" );
    }
    else {
      logWarning("initDistributedMemoryConfiguration()", "ranks are not allowed to skip any reduction (might harm performance). Use optimisation section to switch feature on" );
    }

    tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();
  }
  #endif
}


void exahype::runners::Runner::shutdownDistributedMemoryConfiguration() {
#ifdef Parallel
  tarch::parallel::NodePool::getInstance().terminate();
  exahype::repositories::RepositoryFactory::getInstance().shutdownAllParallelDatatypes();
#endif
}

void exahype::runners::Runner::initSharedMemoryConfiguration() {
  #ifdef SharedMemoryParallelisation

  #ifdef TBBInvade
  tarch::multicore::Core::getInstance().configure( shminvade::SHMController::getInstance().getMaxAvailableCores() );
  #elif SharedTBB
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads,_parser.getThreadStackSize());
  #elif SharedCPP
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);
  #else
  #error Unknown shared memory variant
  #endif

  if ( _parser.useManualPinning() ) {
    #if defined(TBBInvade)
    logWarning("initSharedMemoryConfiguration()", "TBBInvade always pins threads automatically, i.e. manual pinning is ignored" );
    #elif defined(SharedTBB) || defined(SharedCPP)
    logInfo("initSharedMemoryConfiguration()", "manual pinning switched on" );
    tarch::multicore::Core::getInstance().pinThreads( true );
    #else
    logWarning("initSharedMemoryConfiguration()", "manual pinning only supported for TBB" );
    #endif
  }

  // neighbour merge task sets
  exahype::Vertex::SpawnNeighbourMergeAsThread = _parser.getSpawnNeighbourMergeAsThread();

  // background jobs
  solvers::Solver::MaxNumberOfRunningBackgroundJobConsumerTasksDuringTraversal = _parser.getNumberOfBackgroundJobConsumerTasks();
  tarch::multicore::jobs::Job::setMaxNumberOfRunningBackgroundThreads(solvers::Solver::MaxNumberOfRunningBackgroundJobConsumerTasksDuringTraversal);

  #if defined(SharedTBB)
  tarch::multicore::jobs::setMinMaxNumberOfJobsToConsumeInOneRush(
      _parser.getMinBackgroundJobsInARush(), _parser.getMaxBackgroundJobsInARush() );

  if ( _parser.getProcessHighPriorityBackgroundJobsInAnRush() ) { // high priority behaviour
    if ( _parser.getRunLowPriorityJobsOnlyIfNoHighPriorityJobIsLeft() ) { // low priority behaviour
      tarch::multicore::jobs::setHighPriorityJobBehaviour(
          tarch::multicore::jobs::HighPriorityTaskProcessing::ProcessAllHighPriorityTasksInARushAndRunBackgroundTasksOnlyIfNoHighPriorityTasksAreLeft);
    } else {
      tarch::multicore::jobs::setHighPriorityJobBehaviour(
          tarch::multicore::jobs::HighPriorityTaskProcessing::ProcessAllHighPriorityTasksInARush);
    }
  } else if ( _parser.getSpawnHighPriorityBackgroundJobsAsATask() ) {
    if ( _parser.getRunLowPriorityJobsOnlyIfNoHighPriorityJobIsLeft() ) { // low priority behaviour
      logWarning("initSharedMemoryConfiguration()","There exists no high priority job queue if we spawn high priority jobs directly as TBB tasks. "<<
                  "Fall back to 'run_always' low priority job processing strategy.");
    }
    if ( _parser.getSpawnLowPriorityBackgroundJobsAsATask() ){
      tarch::multicore::jobs::setHighPriorityJobBehaviour(
          tarch::multicore::jobs::HighPriorityTaskProcessing::MapHighPriorityAndBackgroundTasksToRealTBBTasks);
    } else {
      tarch::multicore::jobs::setHighPriorityJobBehaviour(
          tarch::multicore::jobs::HighPriorityTaskProcessing::MapHighPriorityTasksToRealTBBTasks);
    }
  }
  else {
    if ( _parser.getRunLowPriorityJobsOnlyIfNoHighPriorityJobIsLeft() ) {
      tarch::multicore::jobs::setHighPriorityJobBehaviour(
          tarch::multicore::jobs::HighPriorityTaskProcessing::ProcessOneHighPriorityTasksAtATimeAndRunBackgroundTasksOnlyIfNoHighPriorityTasksAreLeft);
    } else {
      tarch::multicore::jobs::setHighPriorityJobBehaviour(
          tarch::multicore::jobs::HighPriorityTaskProcessing::ProcessOneHighPriorityTasksAtATime);
    }
  }
  #endif

  // NOTE: Adjusting the grain size might hurt the intermixing of compute-heavy background jobs, e.g. the PredictionJobs,
  // with bandwith-bound ones, e.g. Riemann solves.

  switch (_parser.getMulticoreOracleType()) {
  case exahype::parser::Parser::MulticoreOracleType::Dummy:
    logInfo("initSharedMemoryConfiguration()",
        "use dummy shared memory oracle");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
      new peano::datatraversal::autotuning::OracleForOnePhaseDummy(
         true,  //   bool useMultithreading                  = true,
         0,     //   int  grainSizeOfUserDefinedRegions      = 0,
         peano::datatraversal::autotuning::OracleForOnePhaseDummy::SplitVertexReadsOnRegularSubtree::Split,
         true, //  bool pipelineDescendProcessing
         true, //   bool pipelineAscendProcessing
         27, //   int  smallestProblemSizeForAscendDescend  = tarch::la::aPowI(DIMENSIONS,3*3*3*3/2),
         3, //   int  grainSizeForAscendDescend          = 3,
         1, //   int  smallestProblemSizeForEnterLeaveCell = tarch::la::aPowI(DIMENSIONS,9/2),
         1, //   int  grainSizeForEnterLeaveCell         = 2,
         1, //   int  smallestProblemSizeForTouchFirstLast = tarch::la::aPowI(DIMENSIONS,3*3*3*3+1),
         1, //   int  grainSizeForTouchFirstLast         = 64,
         1, //   int  smallestProblemSizeForSplitLoadStore = tarch::la::aPowI(DIMENSIONS,3*3*3),
         1  //   int  grainSizeForSplitLoadStore         = 8,
      )
    );
    break;
  case exahype::parser::Parser::MulticoreOracleType::AutotuningWithRestartAndLearning:
    logInfo("initSharedMemoryConfiguration()",
        "use learning autotuning shared memory oracle and allow restarts");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(
          tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1,
          true
        ));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  case exahype::parser::Parser::MulticoreOracleType::AutotuningWithoutLearning:
    logInfo("initSharedMemoryConfiguration()",
        "use autotuning shared memory oracle configuration but disable machine learning algorithm");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(false,false));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  case exahype::parser::Parser::MulticoreOracleType::AutotuningWithLearningButWithoutRestart:
    logInfo("initSharedMemoryConfiguration()",
        "use autotuning shared memory oracle but disable search restarts");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(
          tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1,
          false
        ));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  case exahype::parser::Parser::MulticoreOracleType::GrainSizeSampling:
    logInfo("initSharedMemoryConfiguration()",
        "use shared memory oracle sampling");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling(
            64,
            true    // logarithmicDistribution
        ));
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
    break;
  }
  #endif


  #if  defined(TBBInvade)
  shminvade::SHMSharedMemoryBetweenTasks::getInstance().cleanUp();
  switch ( _parser.getTBBInvadeStrategy() ) {
    case exahype::parser::Parser::TBBInvadeStrategy::Undef:
      logError( "preProcessTimeStepInSharedMemoryEnvironment()", "none or no valid invasion statement found in configuration " << _parser.getSharedMemoryConfiguration() );
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::OccupyAllCores:
      shminvade::SHMStrategy::setStrategy( new shminvade::SHMOccupyAllCoresStrategy() );
	  logInfo( "initSharedMemoryConfiguration()", "selected SHMInvade's OccupyAllCores strategy" );
	  if ( _parser.getRanksPerNode()<=0 or _parser.getRanksPerNode()>=shminvade::SHMController::getInstance().getMaxAvailableCores()) {
		logError( "initSharedMemoryConfiguration()", "no ranks-per-node set. Mandatory for SHMInvade" );
	  }
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::NoInvade:
    case exahype::parser::Parser::TBBInvadeStrategy::NoInvadeButAnalyseDistribution:
    case exahype::parser::Parser::TBBInvadeStrategy::InvadeBetweenTimeSteps:
    case exahype::parser::Parser::TBBInvadeStrategy::InvadeThroughoutComputation:
    case exahype::parser::Parser::TBBInvadeStrategy::InvadeAtTimeStepStartupPlusThroughoutComputation:
      shminvade::SHMStrategy::setStrategy( new shminvade::SHMMultipleRanksPerNodeStrategy() );
	  logInfo( "initSharedMemoryConfiguration()", "selected SHMInvade's MultipleRanksPerNode strategy" );
      break;
  }

  shminvade::SHMController::getInstance().init(_parser.getRanksPerNode(),tarch::parallel::Node::getInstance().getRank());

  // This initialisation with dummies is most likely not required at all
  double localData[3] = { 0.0, 1.0, 1.0 };
  shminvade::SHMSharedMemoryBetweenTasks::getInstance().setSharedUserData(localData,3*sizeof(double));
  logInfo( "initSharedMemoryConfiguration()", "initialised local shared memory region with dummies" );

  #if defined(Asserts)
  for (int k=0; k<_parser.getRanksPerNode(); k++) {
    logInfo( "initSharedMemoryConfiguration()", "getSharedUserData<double>(k,0)=" << (shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,0)) );
    logInfo( "initSharedMemoryConfiguration()", "getSharedUserData<double>(k,1)=" << (shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,1)) );
    logInfo( "initSharedMemoryConfiguration()", "getSharedUserData<double>(k,2)=" << (shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,2)) );
  }
  #endif
  #endif
}


void exahype::runners::Runner::initDataCompression() {
  exahype::solvers::Solver::CompressionAccuracy = _parser.getDoubleCompressionFactor();

  if (exahype::solvers::Solver::CompressionAccuracy==0.0) {
    logInfo( "initDataCompression()", "switched off any data compression");
  }
  else {
    exahype::solvers::Solver::SpawnCompressionAsBackgroundJob = _parser.getSpawnDoubleCompressionAsBackgroundTask();
    logInfo( "initDataCompression()", "store all data with accuracy of " << exahype::solvers::Solver::CompressionAccuracy << ". Use background threads for data conversion=" << exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundJob);
  }
}


void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
  #ifdef SharedMemoryParallelisation
  tarch::multicore::jobs::plotStatistics();

  switch (_parser.getMulticoreOracleType()) {
  case exahype::parser::Parser::MulticoreOracleType::AutotuningWithoutLearning:
    break;
  case exahype::parser::Parser::MulticoreOracleType::Dummy:
  case exahype::parser::Parser::MulticoreOracleType::AutotuningWithRestartAndLearning:
  case exahype::parser::Parser::MulticoreOracleType::AutotuningWithLearningButWithoutRestart:
  case exahype::parser::Parser::MulticoreOracleType::GrainSizeSampling:
    #ifdef Parallel
    if (
      tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1
      &&
      tarch::multicore::Core::getInstance().getNumberOfThreads()>1
    ) {
      logInfo("shutdownSharedMemoryConfiguration()",
          "wrote statistics into file " << _parser.getMulticorePropertiesFile()
          << ". Dump from all other ranks subpressed to avoid file races"
      );
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
    }
    #else
    if ( tarch::multicore::Core::getInstance().getNumberOfThreads()>1 ) {
      logInfo("shutdownSharedMemoryConfiguration()",
        "wrote statistics into file "
        << _parser.getMulticorePropertiesFile());
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
        _parser.getMulticorePropertiesFile());
    }
    #endif
    break;
  }
  #endif
}

int exahype::runners::Runner::getCoarsestGridLevelOfAllSolvers(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  double hMax = exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();

  const int peanoLevel = exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(hMax,boundingBoxSize[0]).second;

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << peanoLevel << " (1 means a single cell)");
  return peanoLevel;
}

int exahype::runners::Runner::getCoarsestGridLevelForLoadBalancing(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  return std::max( 3, getCoarsestGridLevelOfAllSolvers(boundingBoxSize) );
}

int exahype::runners::Runner::getFinestUniformGridLevelOfAllSolvers(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  double hMax = exahype::solvers::Solver::getFinestMaximumMeshSizeOfAllSolvers();

  const int peanoLevel = exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(hMax,boundingBoxSize[0]).second;

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << peanoLevel << " (1 means a single cell)");
  return peanoLevel;
}

int exahype::runners::Runner::getFinestUniformGridLevelForLoadBalancing(
    tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize) const {
  return std::max( 3, getFinestUniformGridLevelOfAllSolvers(boundingBoxSize) );
}

double
exahype::runners::Runner::determineCoarsestMeshSize(
    const tarch::la::Vector<DIMENSIONS, double>& boundingBoxSize) const {
  const double coarsestUserMeshSize = exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();
  const double maxBoundingBoxExtent = tarch::la::max(boundingBoxSize);

  const int coarsestMeshLevel =
      exahype::solvers::Solver::computeCoarsestMeshSizeAndLevel(coarsestUserMeshSize,maxBoundingBoxExtent).second;
  return maxBoundingBoxExtent / threePowI(coarsestMeshLevel-1);
}

tarch::la::Vector<DIMENSIONS, double>
exahype::runners::Runner::determineScaledDomainSize(
    const tarch::la::Vector<DIMENSIONS, double>& domainSize,
    const double meshSize) const {
  tarch::la::Vector<DIMENSIONS, double> scaledDomainSize =
      domainSize / meshSize;
  for(int i=0; i<DIMENSIONS; i++) {
    scaledDomainSize[i] = std::ceil(scaledDomainSize[i]) * meshSize;
  }
  return scaledDomainSize;
}

/**
 * @return Bounding box size. If we have a non-cubical domain,
 *         then the bounding box still is cubical and all of its entries are
 *         the biggest dimension along one coordinate axis.
 */
tarch::la::Vector<DIMENSIONS, double> exahype::runners::Runner::determineBoundingBoxSize(
    const tarch::la::Vector<DIMENSIONS, double>& domainSize) const {
    double longestH = tarch::la::max(domainSize);
    return tarch::la::Vector<DIMENSIONS, double>(longestH);
}

exahype::repositories::Repository* exahype::runners::Runner::createRepository() {
  // Geometry is static as we need it to survive the whole simulation time.
  _domainOffset = _parser.getOffset();
  _domainSize   = _parser.getDomainSize();
  _boundingBoxSize  = determineBoundingBoxSize(_domainSize);

  const int coarsestUserMeshLevel = getCoarsestGridLevelOfAllSolvers(_boundingBoxSize);
  int boundingBoxMeshLevel = coarsestUserMeshLevel;
  tarch::la::Vector<DIMENSIONS,double> boundingBoxOffset = _domainOffset;

  #ifdef Parallel
  if (exahype::State::VirtuallyExpandBoundingBox) {
    const double coarsestUserMeshSpacing =
        exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();
    const double maxDomainExtent = tarch::la::max(_domainSize);

    double boundingBoxScaling         = 0;
    double boundingBoxExtent          = 0;
    double boundingBoxMeshSpacing     = std::numeric_limits<double>::max();

    int level = coarsestUserMeshLevel; // level=1 means a single cell
    while (boundingBoxMeshSpacing > coarsestUserMeshSpacing) {
      const double boundingBoxMeshCells = std::pow(3,level-1);
      boundingBoxScaling                = boundingBoxMeshCells / ( boundingBoxMeshCells - 2 );
      boundingBoxExtent                 = boundingBoxScaling * maxDomainExtent;
      boundingBoxMeshSpacing            = boundingBoxExtent/boundingBoxMeshCells;
      level++;
    }
    level--; // decrement result since boundingBox was computed using level-1

    assertion6(boundingBoxScaling>=1.0,boundingBoxScaling,boundingBoxExtent,boundingBoxMeshSpacing,boundingBoxMeshLevel,coarsestUserMeshSpacing,maxDomainExtent);

    boundingBoxMeshLevel = level;
    _boundingBoxSize    *= boundingBoxScaling;
    boundingBoxOffset   -= boundingBoxMeshSpacing;
  }
  #endif

  const double coarsestUserMeshSize = exahype::solvers::Solver::getCoarsestMaximumMeshSizeOfAllSolvers();
  const double coarsestMeshSize     = determineCoarsestMeshSize(_boundingBoxSize);
  tarch::la::Vector<DIMENSIONS,double> scaledDomainSize =
      determineScaledDomainSize(_domainSize,coarsestMeshSize);

  if ( tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank() ) {
    if (!tarch::la::equals(_domainSize,scaledDomainSize)) {
      logInfo("createRepository(...)",
          "scale domain size artificially to " << scaledDomainSize << " from "
          << _domainSize << " since non-cubic domain was specified");
    }
    logInfo("createRepository(...)",
        "coarsest mesh size was chosen as " << coarsestMeshSize << " based on user's maximum mesh size "<<
        coarsestUserMeshSize << " and length of longest edge of domain " << tarch::la::max(scaledDomainSize));
    if (boundingBoxMeshLevel!=coarsestUserMeshLevel) {
      logInfo("createRepository(...)",
          "We will need to refine the grid " << boundingBoxMeshLevel-coarsestUserMeshLevel << " more time(s) than expected "
          " in order to satisfy user's maximum mesh size criterion while scaling the bounding box");
    }

    logInfo(
        "createRepository(...)",
        "summary: create computational domain at " << _domainOffset <<
        " of width/size " << scaledDomainSize <<
        ". bounding box has offset " << boundingBoxOffset <<
        " and size " << _boundingBoxSize <<
        ". grid regular up to level " << boundingBoxMeshLevel << " (1 means a single cell)");
  }

  _domainSize = scaledDomainSize;

  static peano::geometry::Hexahedron geometry(
      _domainSize,
      _domainOffset);

  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _boundingBoxSize,
      boundingBoxOffset
  );
}

void exahype::runners::Runner::initHeaps() {
  exahype::DataHeap::getInstance().setName("DataHeap");
  exahype::CompressedDataHeap::getInstance().setName("compressed-data");
  exahype::solvers::ADERDGSolver::Heap::getInstance().setName("ADERDGCellDescriptionHeap");
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().setName("FiniteVolumesCellDescriptionHeap");
  if ( tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank() ) {
    logInfo("initHeaps()","initialised DataHeap="<<exahype::DataHeap::getInstance().toString());
    logInfo("initHeaps()","initialised CompressedDataHeap::Heap="<<exahype::CompressedDataHeap::getInstance().toString());
    logInfo("initHeaps()","initialised ADERDGSolver::Heap="<<exahype::solvers::ADERDGSolver::Heap::getInstance().toString());
    logInfo("initHeaps()","initialised FiniteVolumesSolver::Heap="<<exahype::solvers::FiniteVolumesSolver::Heap::getInstance().toString());
  }
  #ifdef Parallel
  exahype::MetadataHeap::getInstance().setName("MetadataHeap");
  if ( tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank() ) {
    logInfo("initHeaps()","initialised MetadataHeap="<<exahype::MetadataHeap::getInstance().toString());
  }
  #endif
}

void exahype::runners::Runner::shutdownHeaps() {
  logInfo("shutdownHeaps()","shutdown all heaps");
  exahype::DataHeap::getInstance().shutdown();
  exahype::solvers::ADERDGSolver::Heap::getInstance().shutdown();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().shutdown();
  #ifdef Parallel
  exahype::MetadataHeap::getInstance().shutdown();
  #endif
}

void exahype::runners::Runner::initHPCEnvironment() {
  peano::performanceanalysis::Analysis::getInstance().enable(false);

  solvers::Solver::SwitchOffNeighbourMergePerformedCheck = _parser.getProfilingTarget()==parser::Parser::ProfilingTarget::Update;
}

void exahype::runners::Runner::initOptimisations() const {
  exahype::solvers::Solver::FuseADERDGPhases         = _parser.getFuseAlgorithmicSteps();
  exahype::solvers::Solver::WeightForPredictionRerun = _parser.getFuseAlgorithmicStepsFactor();

  exahype::solvers::Solver::configurePredictionPhase(
      _parser.getSpawnPredictionAsBackgroundThread(),
      _parser.getSpawnProlongationAsBackgroundThread());

  exahype::solvers::Solver::SpawnUpdateAsBackgroundJob =
      _parser.getSpawnUpdateAsBackgroundThread();

  exahype::solvers::Solver::SpawnAMRBackgroundJobs =
      _parser.getSpawnAMRBackgroundThreads();

  exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps =
      _parser.getDisablePeanoNeighbourExchangeInTimeSteps();
  exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps =
      _parser.getDisableMetadataExchangeInBatchedTimeSteps();

  if ( tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank() ) {
    logInfo("parseOptimisations()","use the following global optimisations:");
      logInfo("parseOptimisations()","\tfuse-algorithmic-steps="                  << (exahype::solvers::Solver::FuseADERDGPhases ? "on" : "off"));
      logInfo("parseOptimisations()","\tfuse-algorithmic-steps-factor="           << exahype::solvers::Solver::WeightForPredictionRerun);
      logInfo("parseOptimisations()","\tspawn-predictor-as-background-thread="    << (exahype::solvers::Solver::SpawnPredictionAsBackgroundJob ? "on" : "off"));
      logInfo("parseOptimisations()","\tspawn-prolongation-as-background-thread=" << (exahype::solvers::Solver::SpawnProlongationAsBackgroundJob ? "on" : "off"));
      logInfo("parseOptimisations()","\tspawn-update-as-background-thread="       << (exahype::solvers::Solver::SpawnUpdateAsBackgroundJob ? "on" : "off"));
      logInfo("parseOptimisations()","\tspawn-amr-background-threads="            << (exahype::solvers::Solver::SpawnAMRBackgroundJobs ? "on" : "off"));
      logInfo("parseOptimisations()","\tdisable-vertex-exchange-in-time-steps="   << (exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps ? "on" : "off"));
      logInfo("parseOptimisations()","\tbatching enabled="<<
          ( _parser.getSkipReductionInBatchedTimeSteps() &&
              exahype::solvers::Solver::allSolversUseTimeSteppingScheme(exahype::solvers::Solver::TimeStepping::GlobalFixed)
              ? "yes" : "no" ) );
      logInfo("parseOptimisations()","\t\ttime-step-batch-factor="<<_parser.getTimestepBatchFactor());
      logInfo("parseOptimisations()","\t\tall solvers use 'globalfixed' time stepping="<<
          ( exahype::solvers::Solver::allSolversUseTimeSteppingScheme(exahype::solvers::Solver::TimeStepping::GlobalFixed)
          ? "yes" : "no" ) );
      logInfo("parseOptimisations()","\tdisable-metadata-exchange-in-batched-time-steps="<<(exahype::solvers::Solver::DisableMetaDataExchangeInBatchedTimeSteps ? "on" : "off"));
  }
}

int exahype::runners::Runner::run() {
  int result = 0;
  if ( _parser.isValid() ) {
    initOptimisations();

    initHeaps();

    #ifdef Parallel
    exahype::State::VirtuallyExpandBoundingBox =_parser.getScaleBoundingBox();
    exahype::State::BroadcastInThisIteration   = true;
    exahype::State::ReduceInThisIteration      = false;
    #endif

    auto* repository = createRepository();
    // must come after repository creation
    initSolvers();

    if (_parser.isValid() && _parser.getMeasureCellProcessingTimes() ) {
      measureCellProcessingTimes();
    }

    if ( _parser.isValid() )
      initDistributedMemoryConfiguration();
    if ( _parser.isValid() )
      initSharedMemoryConfiguration();
    if ( _parser.isValid() )
      initDataCompression();
    if ( _parser.isValid() )
      initHPCEnvironment();

    if ( _parser.isValid() ) {
      if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
        result = runAsMaster(*repository);
      }
      #ifdef Parallel
      else {
        result = runAsWorker(*repository);
      }
      #endif
    }

    if ( _parser.isValid() )
      shutdownSharedMemoryConfiguration();
    if ( _parser.isValid() )
      shutdownDistributedMemoryConfiguration();

    shutdownHeaps();

    delete repository;
  }
  else {
    logError( "run(...)", "do not run code as parser reported errors" );
    result = 1;
  }

  return result;
}

void exahype::runners::Runner::initSolvers() const {
  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    solver->initSolver(
      0.0,_domainOffset,_domainSize,_boundingBoxSize,
      _cmdlineargs,_parser.createParserView(solverNumber)
    );
  }

  // profiling functionality
  if ( _parser.getMeasureCellProcessingTimes() ) {

  }
}

void exahype::runners::Runner::measureCellProcessingTimes() const {
  int solverNumber = 0;
  for ( auto* solver : solvers::RegisteredSolvers ) {
    solvers::Solver::CellProcessingTimes result =
        solver->measureCellProcessingTimes(_parser.getMeasureCellProcessingTimesIterations());

    std::ostringstream stringstr;
    stringstr << "cell processing times for solver "<<solverNumber<<":"<<std::endl;
    result.toString(stringstr,1e6,3,"\u00B5s","\t\t");
    logInfo("measureCellProcessingTimes()",stringstr.str());

    solverNumber++;
  }
}

void exahype::runners::Runner::printMeshSetupInfo(
    exahype::repositories::Repository& repository,
    const int meshSetupIterations) const {
  #if defined(TrackGridStatistics) && defined(Asserts)
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", max-level=" << repository.getState().getMaxLevel() <<
      ", state=" << repository.getState().toString() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() <<
      ", vertical solver communication=" << repository.getState().getVerticalExchangeOfSolverDataRequired() <<
      ", continue to construct grid=" << repository.getState().continueToConstructGrid()
  );
  #elif defined(Asserts)
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", state=" << repository.getState().toString() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() <<
      ", vertical solver communication=" << repository.getState().getVerticalExchangeOfSolverDataRequired() <<
      ", continue to construct grid=" << repository.getState().continueToConstructGrid()
  );
  #elif defined(TrackGridStatistics)
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", max-level=" << repository.getState().getMaxLevel() <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() <<
      ", vertical solver communication=" << repository.getState().getVerticalExchangeOfSolverDataRequired() <<
      ", continue to construct grid=" << repository.getState().continueToConstructGrid()
  );
  #else
  logInfo("createGrid()",
      "grid setup iteration #" << meshSetupIterations <<
      ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() <<
      ", vertical solver communication=" << repository.getState().getVerticalExchangeOfSolverDataRequired() <<
      ", continue to construct grid=" << repository.getState().continueToConstructGrid() );
  #endif

  #if !defined(Parallel)
  logInfo("createGrid(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  #else
  if (tarch::parallel::Node::getInstance().getNumberOfNodes()==1) {
    logInfo("createGrid(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  }
  #endif

  #ifdef Asserts
  if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
    DataHeap::getInstance().plotStatistics();
    peano::heap::PlainCharHeap::getInstance().plotStatistics();
  }
  #endif
}

bool exahype::runners::Runner::createMesh(exahype::repositories::Repository& repository) {
  bool meshUpdate = false;

  int meshSetupIterations = 0;
  repository.switchToMeshRefinement();

  repository.getState().setMeshRefinementHasConverged(false);
  while ( repository.getState().continueToConstructGrid() ) {
    repository.iterate(1,true);
    meshSetupIterations++;

    repository.getState().endedGridConstructionIteration( getFinestUniformGridLevelOfAllSolvers(_boundingBoxSize) );

    printMeshSetupInfo(repository,meshSetupIterations);

    meshUpdate = true;
  }

  logInfo("createGrid(Repository)", "finished grid setup after " << meshSetupIterations << " iterations" );

  if (
    tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0
    &&
    tarch::parallel::Node::getInstance().getNumberOfNodes()>1
  ) {
    logWarning( "createGrid(Repository)", "there are still " << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() << " ranks idle" )
    //  #ifdef Parallel
    //  // Might be too restrictive for later runs. Remove but keep warning from above
    //  assertion( tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()==0 );
    //  #endif
  }

  return meshUpdate;
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  if (!exahype::solvers::RegisteredSolvers.empty()) {
    initialiseMesh(repository);

    logInfo( "runAsMaster(...)", "initialised all data and computed first time step size" );

    bool communicatePeanoVertices =
        !exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps;

    repository.switchToInitialPrediction();
    repository.iterate( exahype::solvers::Solver::PredictionSweeps, communicatePeanoVertices );
    logInfo("runAsMaster(...)","computed first predictor (number of predictor sweeps: "<<exahype::solvers::Solver::PredictionSweeps<<")");

    // configure time stepping loop
    double simulationEndTime   = std::numeric_limits<double>::max();
    int simulationTimeSteps = std::numeric_limits<int>::max();
    if (_parser.foundSimulationEndTime()) {
      simulationEndTime = _parser.getSimulationEndTime();
    } else {
      simulationTimeSteps = _parser.getSimulationTimeSteps();
    }
    parser::Parser::ProfilingTarget profilingTarget = _parser.getProfilingTarget();
    if ( profilingTarget==parser::Parser::ProfilingTarget::WholeCode ) {
      printTimeStepInfo(-1,repository);
      validateInitialSolverTimeStepData(exahype::solvers::Solver::FuseADERDGPhases);
    }
    const bool skipReductionInBatchedTimeSteps  = _parser.getSkipReductionInBatchedTimeSteps();

    // run time stepping loop
    int timeStep = 0;
    while (
        tarch::la::greater(solvers::Solver::getMinTimeStepSizeOfAllSolvers(), 0.0) &&
        solvers::Solver::getMinTimeStampOfAllSolvers() < simulationEndTime         &&
        timeStep < simulationTimeSteps
    ) {
      bool plot = exahype::plotters::checkWhetherPlotterBecomesActive(
          solvers::Solver::getMinTimeStampOfAllSolvers()); // has no side effects

      preProcessTimeStepInSharedMemoryEnvironment();

      // fused time stepping
      int numberOfStepsToRun = 1;
      if ( profilingTarget==parser::Parser::ProfilingTarget::WholeCode && exahype::solvers::Solver::FuseADERDGPhases  ) {
        if (plot) {
          numberOfStepsToRun = 0;
        }
        else if (
            skipReductionInBatchedTimeSteps &&
            solvers::Solver::allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping::GlobalFixed) &&
            repository.getState().getVerticalExchangeOfSolverDataRequired()==false // known after mesh update
        ) {
          numberOfStepsToRun = determineNumberOfBatchedTimeSteps(timeStep);
        }

        runTimeStepsWithFusedAlgorithmicSteps(repository,numberOfStepsToRun);
        printTimeStepInfo(numberOfStepsToRun,repository);
      // straightforward realisation
      } else if ( profilingTarget==parser::Parser::ProfilingTarget::WholeCode ) {
        runOneTimeStepWithThreeSeparateAlgorithmicSteps(repository, plot);
      }
      // profiling of isolated adapters
      else if ( profilingTarget==parser::Parser::ProfilingTarget::Prediction ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun "<<solvers::Solver::PredictionSweeps<<" iteration with PredictionRerun adapter");
        printGridStatistics(repository);
        runPredictionInIsolation(repository);
      } else if ( profilingTarget==parser::Parser::ProfilingTarget::NeigbhourMerge ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun one iteration with MergeNeighours adapter");
        printGridStatistics(repository);
        repository.switchToMergeNeighbours();
        repository.iterate(1,false);
      } else if ( profilingTarget==parser::Parser::ProfilingTarget::Update ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun one iteration with UpdateAndReduce adapter");
        printGridStatistics(repository);
        repository.switchToUpdateAndReduce();
        repository.iterate(1,false);
      }

      postProcessTimeStepInSharedMemoryEnvironment();

      #if !defined(Parallel)
      logInfo("runAsMaster(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
      #else
      if (tarch::parallel::Node::getInstance().getNumberOfNodes()==1) {
        logInfo("runAsMaster(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
      }
      #endif

      timeStep += numberOfStepsToRun==0 ? 1 : numberOfStepsToRun;
    }

    if ( tarch::la::equals(solvers::Solver::getMinTimeStepSizeOfAllSolvers(), 0.0)) {
      logWarning("runAsMaster(...)","Minimum solver time step size is zero (up to machine precision).");
    }

    repository.switchToBroadcastAndDropNeighbourMessages();
    repository.iterate( 1, communicatePeanoVertices );

    printStatistics();
    repository.logIterationStatistics(false);
  }

  repository.terminate();

  return 0;
}

int exahype::runners::Runner::determineNumberOfBatchedTimeSteps(const int& currentTimeStep) {
  int batchSize = 0;

  double simulationEndTime   = std::numeric_limits<double>::max();
  int simulationTimeSteps = std::numeric_limits<int>::max();
  if (_parser.foundSimulationEndTime()) {
    simulationEndTime = _parser.getSimulationEndTime();
  } else {
    simulationTimeSteps = _parser.getSimulationTimeSteps();
  }

  const double minTimeStepSize          = solvers::Solver::getMinTimeStepSizeOfAllSolvers();
  const double maxTimeStamp             = solvers::Solver::getMaxTimeStampOfAllSolvers();
  const double timeIntervalTillNextPlot = exahype::plotters::getTimeOfNextPlot() - maxTimeStamp;
  const bool   haveActivePlotters       = exahype::plotters::getTimeOfNextPlot() < std::numeric_limits<double>::max();

  if (_parser.foundSimulationEndTime()) {
    const double timeIntervalTillEndTime = simulationEndTime - maxTimeStamp;

    batchSize =
      static_cast<int>(
        std::min (
            _parser.getTimestepBatchFactor() * std::min(timeIntervalTillNextPlot, simulationEndTime),
             timeIntervalTillEndTime
        ) /  minTimeStepSize
      );

  } else {
    batchSize =
      std::min(
        static_cast<int>(simulationTimeSteps * _parser.getTimestepBatchFactor()),
        simulationTimeSteps - currentTimeStep
      );
    if (haveActivePlotters) {
      const int stepsTillNextPlot = static_cast<int>(
        _parser.getTimestepBatchFactor() *
        timeIntervalTillNextPlot / minTimeStepSize
      );
      batchSize = std::min(batchSize,stepsTillNextPlot);
    }
  }
  return batchSize<1 ? 1 : batchSize;
}

void exahype::runners::Runner::preProcessTimeStepInSharedMemoryEnvironment() {
  #if defined(TBBInvade)
  const int myIndexWithinSharedUserData  = shminvade::SHMSharedMemoryBetweenTasks::getInstance().getProcessIndexInSharedDataTable();
  const int ranksOnThisNode              = shminvade::SHMSharedMemoryBetweenTasks::getInstance().getNumberOfRegisteredProcesses();

  #ifdef Asserts
  logInfo(
    "preProcessTimeStepInSharedMemoryEnvironment()",
	"ranks on this node=" << ranksOnThisNode <<
	", my index=" << myIndexWithinSharedUserData
  );
  #endif

  std::vector<double>  t1;
  std::vector<double>  f;
  std::vector<double>  s;
  for (int k=0; k<ranksOnThisNode; k++) {
    t1.push_back( shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,0) );
    f.push_back(  shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,1) );
    s.push_back(  shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,2) );

    #ifdef Asserts
    logInfo( "preProcessTimeStepInSharedMemoryEnvironment()", "getSharedUserData<double>(k,0)=" << (shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,0)) );
    logInfo( "preProcessTimeStepInSharedMemoryEnvironment()", "getSharedUserData<double>(k,1)=" << (shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,1)) );
    logInfo( "preProcessTimeStepInSharedMemoryEnvironment()", "getSharedUserData<double>(k,2)=" << (shminvade::SHMSharedMemoryBetweenTasks::getInstance().getSharedUserData<double>(k,2)) );
    #endif
  }

  switch ( _parser.getTBBInvadeStrategy() ) {
    case exahype::parser::Parser::TBBInvadeStrategy::Undef:
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::NoInvade:
      {
    	if ( _shmInvade==nullptr ) {
    	  const int cores = shminvade::SHMController::getInstance().getMaxAvailableCores() / _parser.getRanksPerNode() -1;
          logInfo(
            "preProcessTimeStepInSharedMemoryEnvironment()",
			"try to acquire SHMInvade object for " << cores << " core(s) in runner (max cores=" <<
			shminvade::SHMController::getInstance().getMaxAvailableCores() << ", ranks-per-node=" << _parser.getRanksPerNode() <<
			")"
	      );
          _shmInvade = new shminvade::SHMInvade( cores );
    	}
        shminvade::SHMController::getInstance().switchOff();
      }
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::NoInvadeButAnalyseDistribution:
      {
    	if ( _shmInvade==nullptr ) {
    	  const int cores = shminvade::SHMController::getInstance().getMaxAvailableCores() / _parser.getRanksPerNode() -1;
          logInfo(
            "preProcessTimeStepInSharedMemoryEnvironment()",
			"try to acquire SHMInvade object for " << cores << " core(s) in runner (max cores=" <<
			shminvade::SHMController::getInstance().getMaxAvailableCores() << ", ranks-per-node=" << _parser.getRanksPerNode() <<
			")"
	      );
          _shmInvade = new shminvade::SHMInvade( cores );
    	}
        shminvade::SHMController::getInstance().switchOff();

        // ask for an optimal number of cores for local rank
        int optimalNumberOfThreads = std::max(
          peano::performanceanalysis::SpeedupLaws::getOptimalNumberOfThreads(
            myIndexWithinSharedUserData,
            t1,f,s,
            shminvade::SHMController::getInstance().getMaxAvailableCores(),
            tarch::parallel::Node::getInstance().getRank() % _parser.getRanksPerNode()
          ),
          2
        );

        logInfo(
          "preProcessTimeStepInSharedMemoryEnvironment()",
			"optimal number of cores would be " << optimalNumberOfThreads
	      );
      }
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::OccupyAllCores:
      {
      	if ( _shmInvade==nullptr ) {
      	  const int cores = shminvade::SHMInvade::MaxCores;
          logInfo( "preProcessTimeStepInSharedMemoryEnvironment()", "try to acquire SHMInvade object for " << cores << " (all cores) in runner" );
          _shmInvade = new shminvade::SHMInvade( cores );
      	}
        shminvade::SHMController::getInstance().switchOff();
      }
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::InvadeBetweenTimeSteps:
      {
        // ask for an optimal number of cores for local rank
        int optimalNumberOfThreads = std::max(
          peano::performanceanalysis::SpeedupLaws::getOptimalNumberOfThreads(
            myIndexWithinSharedUserData,
            t1,f,s,
            shminvade::SHMController::getInstance().getMaxAvailableCores(),
            tarch::parallel::Node::getInstance().getRank() % _parser.getRanksPerNode()
          ),
          2
        );
        assertion(_shmInvade==nullptr);
        _shmInvade = new shminvade::SHMInvade( optimalNumberOfThreads );
        shminvade::SHMController::getInstance().switchOff();
      }
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::InvadeThroughoutComputation:
        assertion(_shmInvade==nullptr);
      _shmInvade = new shminvade::SHMInvade( 2 );
      shminvade::SHMController::getInstance().switchOn();
      break;
    case exahype::parser::Parser::TBBInvadeStrategy::InvadeAtTimeStepStartupPlusThroughoutComputation:
      {
        // ask for an optimal number of cores for local rank
        int optimalNumberOfThreads = std::max(
          peano::performanceanalysis::SpeedupLaws::getOptimalNumberOfThreads(
            myIndexWithinSharedUserData,
            t1,f,s,
            shminvade::SHMController::getInstance().getMaxAvailableCores(),
            tarch::parallel::Node::getInstance().getRank() % _parser.getRanksPerNode()
          ),
          2
        );
        assertion(_shmInvade==nullptr);
        _shmInvade = new shminvade::SHMInvade( optimalNumberOfThreads );
        shminvade::SHMController::getInstance().switchOn();
      }
      break;
  }

  #endif
}


void exahype::runners::Runner::postProcessTimeStepInSharedMemoryEnvironment() {
  #if  defined(SharedMemoryParallelisation) && defined(PerformanceAnalysis) && !defined(Parallel)
  if (sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::hasLearnedSinceLastQuery()) {
    static int dumpCounter = -1;
    dumpCounter++;
    peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics( _parser.getMulticorePropertiesFile() + "-dump-" + std::to_string(dumpCounter) );
  }
  #endif

  #if  defined(TBBInvade)
  static tarch::timing::Watch                     invasionWatch("exahype::Runner", "postProcessTimeStepInSharedMemoryEnvironment()", false);
  static peano::performanceanalysis::SpeedupLaws  amdahlsLaw;

  switch ( _parser.getTBBInvadeStrategy() ) {
     case exahype::parser::Parser::TBBInvadeStrategy::Undef:
     case exahype::parser::Parser::TBBInvadeStrategy::NoInvade:
     case exahype::parser::Parser::TBBInvadeStrategy::OccupyAllCores:
       break;
     case exahype::parser::Parser::TBBInvadeStrategy::InvadeBetweenTimeSteps:
     case exahype::parser::Parser::TBBInvadeStrategy::NoInvadeButAnalyseDistribution:
     case exahype::parser::Parser::TBBInvadeStrategy::InvadeAtTimeStepStartupPlusThroughoutComputation:
       {
      	 logInfo( "postProcessTimeStepInSharedMemoryEnvironment()", "release SHMInvade object" );
         assertion( _shmInvade != nullptr );
         delete _shmInvade;
         _shmInvade = nullptr;

    	 // adopt my local performance model
    	 invasionWatch.stopTimer();
    	 amdahlsLaw.addMeasurement(
           shminvade::SHMController::getInstance().getBookedCores(),
    	   invasionWatch.getCalendarTime()
    	 );
    	 amdahlsLaw.relaxAmdahlsLawWithThreadStartupCost();

    	 //
    	 //
    	 //
    	 double localData[3] = { amdahlsLaw.getSerialTime(), amdahlsLaw.getSerialCodeFraction(), amdahlsLaw.getStartupCostPerThread() };
    	 shminvade::SHMSharedMemoryBetweenTasks::getInstance().setSharedUserData(localData,3);

    	 logDebug( "postProcessTimeStepInSharedMemoryEnvironment()", "localData[0]=" << localData[0] );
    	 logDebug( "postProcessTimeStepInSharedMemoryEnvironment()", "localData[1]=" << localData[1] );
    	 logDebug( "postProcessTimeStepInSharedMemoryEnvironment()", "localData[2]=" << localData[2] );

    	 assertion2( amdahlsLaw.getSerialTime()>=0.0,           amdahlsLaw.getSerialTime(),            amdahlsLaw.toString() );
    	 assertion2( amdahlsLaw.getSerialCodeFraction()>=0.0,   amdahlsLaw.getSerialCodeFraction(),    amdahlsLaw.toString() );
    	 assertion2( amdahlsLaw.getSerialCodeFraction()<=1.0,   amdahlsLaw.getSerialCodeFraction(),    amdahlsLaw.toString() );
         assertion2( amdahlsLaw.getStartupCostPerThread()>=0.0, amdahlsLaw.getStartupCostPerThread(),  amdahlsLaw.toString() );
       }
       break;
     case exahype::parser::Parser::TBBInvadeStrategy::InvadeThroughoutComputation:
       {
      	 logInfo( "postProcessTimeStepInSharedMemoryEnvironment()", "release SHMInvade object" );
         assertion( _shmInvade != nullptr );
         delete _shmInvade;
         _shmInvade = nullptr;
       }
       break;
  }

  invasionWatch.startTimer();
  shminvade::SHMController::getInstance().switchOn();
  #endif
}


void exahype::runners::Runner::updateStatistics() {
  _localRecomputations  +=  (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) ? 1 : 0;
  _meshRefinements      +=  (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedGlobalRecomputation()) ? 1 : 0;
  _predictorReruns      +=  (exahype::solvers::Solver::oneSolverViolatedStabilityCondition()) ? 1 : 0;
}

void exahype::runners::Runner::printStatistics() {
  logInfo("printStatistics(...)","number of mesh refinements      = "<<_meshRefinements);
  logInfo("printStatistics(...)","number of local recomputations  = "<<_localRecomputations);
  logInfo("printStatistics(...)","number of predictor reruns      = "<<_predictorReruns);
}

void exahype::runners::Runner::validateInitialSolverTimeStepData(const bool fuseADERDGPhases) const {
  #ifdef Asserts
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    assertionEquals(solver->getMinTimeStamp(),0.0);
    assertion1(std::isfinite(solver->getMinTimeStepSize()),solver->getMinTimeStepSize());
    assertion1(solver->getMinTimeStepSize()>0,solver->getMinTimeStepSize());

    switch(solver->getTimeStepping()) {
      case exahype::solvers::Solver::TimeStepping::Global:
        assertionEquals(solver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
        break;
      case exahype::solvers::Solver::TimeStepping::GlobalFixed:
        break;
    }
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG: {
        auto* aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        if (!exahype::solvers::Solver::FuseADERDGPhases) {
          assertionEquals(aderdgSolver->getPreviousMinCorrectorTimeStepSize(),0.0); // TOOD(Dominic): Revision
        }
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertionEquals(aderdgSolver->getMinCorrectorTimeStamp(),0.0);
        if (fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinPredictorTimeStepSize());
        } else {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),0.0);
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      }
      break;
      case exahype::solvers::Solver::Type::LimitingADERDG: {
        // ADER-DG
        auto* aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        if (!exahype::solvers::Solver::FuseADERDGPhases) {
          assertionEquals(aderdgSolver->getPreviousMinCorrectorTimeStepSize(),0.0); // TODO(Dominic): Revision
        }
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertionEquals(aderdgSolver->getMinCorrectorTimeStamp(),0.0);
        if (fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinPredictorTimeStepSize());
        } else {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),0.0);
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      } break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        auto* finiteVolumesSolver = static_cast<exahype::solvers::FiniteVolumesSolver*>(solver);
        if (!exahype::solvers::Solver::FuseADERDGPhases) {
          assertionEquals(finiteVolumesSolver->getPreviousMinTimeStepSize(),0.0);
        }
        break;
    }
  }
  #endif
}

void exahype::runners::Runner::initialiseMesh(exahype::repositories::Repository& repository) {
  // We refine here using the previous solution (which is valid)
  logInfo("initialiseMesh(...)","create initial grid");
  createMesh(repository);

  logInfo("initialiseMesh(...)","finalise mesh refinement and compute first time step size");
  repository.switchToFinaliseMeshRefinement();
  repository.iterate();
}

void exahype::runners::Runner::updateMeshOrLimiterDomain(
    exahype::repositories::Repository& repository, const bool fusedTimeStepping) {
  // 1. All solvers drop their MPI messages and broadcast time step data
  repository.switchToBroadcastAndDropNeighbourMessages();
  repository.iterate(1,false);

  // 2. Only the solvers with irregular limiter domain change do the limiter status spreading.
  if ( exahype::solvers::Solver::oneSolverRequestedRefinementStatusSpreading() ) {
    logInfo("updateMeshAndSubdomains(...)","pre-spreading of limiter status");
    repository.switchToRefinementStatusSpreading();
    repository.iterate(
        exahype::solvers::Solver::getMaxRefinementStatus()+1,false);
  }

  // 3. Perform a grid update for those solvers that requested refinement
  if ( exahype::solvers::Solver::oneSolverRequestedMeshRefinement() ) {
    logInfo("updateMeshAndSubdomains(...)","perform mesh refinement");
    createMesh(repository);
  }

  // 4. Drop the MPI metadata for all solvers that requested grid refinement
  // Further reinitialse solvers that reported an irregular limiter domain change
  logInfo("updateMeshAndSubdomains(...)","finalise mesh refinement (if applicable)");
  logInfo("updateMeshAndSubdomains(...)","reinitialise cells and send data to neigbours (if applicable)");
  repository.switchToFinaliseMeshRefinementOrLocalRollback();
  repository.iterate(1,false);

  // 5. Compute Prediction for solvers performing mesh update or
  // perform local recomputation for other solvers
  // Do not advance the time stamp if global recomputation/mesh refinement
  // Advance time stamp if a solver performs a local recomputation.
  // Lastly, send out data face data if we perform fused time stepping.
  // Otherwise, just sent time step data.
  // (Note that we compute the predictor locally for the solvers performing
  // a local prediction as well. They are ready to sent here as well.)
  if (fusedTimeStepping ||
      exahype::solvers::Solver::oneSolverRequestedLocalRecomputation()) {
    logInfo("updateMeshAndSubdomains(...)","recompute solution locally (if applicable) and compute new time step size");
    repository.switchToPredictionOrLocalRecomputation(); // do not roll forward here if global recomp.; we want to stay at the old time step
    const int sweeps = (exahype::solvers::Solver::FuseADERDGPhases) ? exahype::solvers::Solver::PredictionSweeps : 1;
    repository.iterate( sweeps ,false ); // local recomputation: has now recomputed predictor in interface cells
  }
}

void exahype::runners::Runner::printTimeStepInfo(int numberOfStepsRanSinceLastCall, const exahype::repositories::Repository& repository) {
  double currentMinTimeStamp    = std::numeric_limits<double>::max();
  double currentMinTimeStepSize = std::numeric_limits<double>::max();
  double nextMinTimeStepSize    = std::numeric_limits<double>::max();

  static int n = 0;
  if (numberOfStepsRanSinceLastCall==0) {
    n++;
  }
  else if (numberOfStepsRanSinceLastCall>0) {
    n+=numberOfStepsRanSinceLastCall;
  }

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinTimeStamp());
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinTimeStepSize());
    nextMinTimeStepSize =
        std::min(nextMinTimeStepSize, p->getMinNextTimeStepSize());

    #if defined(Debug) || defined(Asserts)
    switch(p->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: t_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: dt_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: t_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: dt_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStepSize());
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        logInfo("startNewTimeStep(...)",
                 "\tADER-DG prev correction*:  t_min   =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getPreviousMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prev correction*:  dt_min  =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getPreviousMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: t_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: dt_min        =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: t_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinPredictorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: dt_min        =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinPredictorTimeStepSize());
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        break;
    }
    #endif
  }

  logInfo("startNewTimeStep(...)",
      "step " << n << "\tt_min          =" << currentMinTimeStamp);

  logInfo("startNewTimeStep(...)",
      "\tdt_min         =" << currentMinTimeStepSize);

  #if !defined(Parallel)
  // memory consumption on rank 0 would not make any sense
  logInfo("startNewTimeStep(...)",
      "\tmemoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");

  if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
    DataHeap::getInstance().plotStatistics();
    peano::heap::PlainCharHeap::getInstance().plotStatistics();

    logInfo(
      "startNewTimeStep(...)",
      "\tpiped-uncompressed-bytes=" << exahype::solvers::ADERDGSolver::PipedUncompressedBytes
      << "\tpiped-compressed-bytes=" << exahype::solvers::ADERDGSolver::PipedCompressedBytes
      << "\tcompression-rate=" << (exahype::solvers::ADERDGSolver::PipedCompressedBytes/exahype::solvers::ADERDGSolver::PipedUncompressedBytes)
    );
  }
  #endif

  #if defined(TrackGridStatistics)
  if (repository.getState().getNumberOfInnerCells()>0 and repository.getState().getMaxLevel()>0) {
    logInfo(
      "startNewTimeStep(...)",
      "\tinner cells/inner unrefined cells=" << repository.getState().getNumberOfInnerCells()
      << "/" << repository.getState().getNumberOfInnerLeafCells() );
    logInfo(
      "startNewTimeStep(...)",
      "\tinner max/min mesh width=" << repository.getState().getMaximumMeshWidth()
      << "/" << repository.getState().getMinimumMeshWidth()
      );
    logInfo(
      "startNewTimeStep(...)",
      "\tmax level=" << repository.getState().getMaxLevel()
      );
  }
  #endif

  if (solvers::Solver::getMinTimeStampOfAllSolvers()>std::numeric_limits<double>::max()/100.0) {
    logError("runAsMaster(...)","quit simulation as solver seems to explode" );
    exit(-1);
  }

  #if defined(Debug) || defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().closeOutputStreamAndReopenNewOne();
  #endif
}

void exahype::runners::Runner::runTimeStepsWithFusedAlgorithmicSteps(
    exahype::repositories::Repository& repository, int numberOfStepsToRun) {

  if (numberOfStepsToRun==0) {
    logInfo("runTimeStepsWithFusedAlgorithmicSteps(...)","plot");
  }
  else if (numberOfStepsToRun>1) {
    logInfo("runTimeStepsWithFusedAlgorithmicSteps(...)","run "<<numberOfStepsToRun<< " time steps within one batch");
  }

  bool communicatePeanoVertices =
      !exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps;
  repository.switchToFusedTimeStep();
  if (numberOfStepsToRun==0) {
    repository.iterate( exahype::solvers::Solver::PredictionSweeps,communicatePeanoVertices );
  } else {
    repository.iterate( exahype::solvers::Solver::PredictionSweeps*numberOfStepsToRun,false/*Always disable during batching*/ );
  }

  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) {
    logInfo("runTimeStepsWithFusedAlgorithmicSteps(...)","local recomputation requested by at least one solver");
  }
  if (exahype::solvers::Solver::oneSolverRequestedMeshRefinement()) {
    logInfo("runTimeStepsWithFusedAlgorithmicSteps(...)","mesh update requested by at least one solver");
  }

  if (exahype::solvers::Solver::oneSolverRequestedMeshRefinement() ||
      exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) {
    updateMeshOrLimiterDomain(repository,true);
  }

  if (exahype::solvers::Solver::oneSolverViolatedStabilityCondition()) {
    logInfo("runTimeStepsWithFusedAlgorithmicSteps(...)", "\t\t recompute space-time predictor");
    repository.switchToPredictionRerun();
    repository.iterate( exahype::solvers::Solver::PredictionSweeps, communicatePeanoVertices );
  }

  updateStatistics();
  // ---- reduction/broadcast barrier ----
}

void exahype::runners::Runner::runOneTimeStepWithThreeSeparateAlgorithmicSteps(
    exahype::repositories::Repository& repository, bool plot) {
  // Only one time step (predictor vs. corrector) is used in this case.
  bool communicatePeanoVertices =
        !exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps;
  repository.switchToMergeNeighbours();  // Riemann -> face2face
  repository.iterate( 1, communicatePeanoVertices );
  repository.switchToUpdateAndReduce();  // Face to cell + Inside cell
  repository.iterate( 1, communicatePeanoVertices );

  if (exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) {
    logInfo("runOneTimeStepWithThreeSeparateAlgorithmicSteps(...)","local recomputation requested by at least one solver");
  }
  if (exahype::solvers::Solver::oneSolverRequestedMeshRefinement()) {
    logInfo("runOneTimeStepWithThreeSeparateAlgorithmicSteps(...)","mesh update requested by at least one solver");
  }

  if (exahype::solvers::Solver::oneSolverRequestedMeshRefinement() ||
      exahype::solvers::LimitingADERDGSolver::oneSolverRequestedLocalRecomputation()) {
    updateMeshOrLimiterDomain(repository,false);
  }

  printTimeStepInfo(1,repository);
  repository.switchToPrediction(); // Cell onto faces
  repository.iterate( exahype::solvers::Solver::PredictionSweeps, communicatePeanoVertices );

  updateStatistics();
}

void exahype::runners::Runner::printGridStatistics(repositories::Repository& repository) {
  #if defined(TrackGridStatistics)
  if (repository.getState().getNumberOfInnerCells()>0 and repository.getState().getMaxLevel()>0) {
    logInfo(
      "printGridStatistics(...)",
      "\tinner cells/inner unrefined cells=" << repository.getState().getNumberOfInnerCells()
      << "/" << repository.getState().getNumberOfInnerLeafCells() );
    logInfo(
      "printGridStatistics(...)",
      "\tinner max/min mesh width=" << repository.getState().getMaximumMeshWidth()
      << "/" << repository.getState().getMinimumMeshWidth()
      );
    logInfo(
      "printGridStatistics(...)",
      "\tmax level=" << repository.getState().getMaxLevel()
      );
  }
  #endif
}

void exahype::runners::Runner::runPredictionInIsolation(repositories::Repository& repository) {
  for (auto* solver : solvers::RegisteredSolvers) {
    switch (solver->getType()) {
      case solvers::Solver::Type::ADERDG:
        static_cast<solvers::ADERDGSolver*>(solver)->setStabilityConditionWasViolated(true);
        static_cast<solvers::ADERDGSolver*>(solver)->setNextMeshUpdateEvent(); // reset
        break;
      case solvers::Solver::Type::LimitingADERDG:
        static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver().get()->setStabilityConditionWasViolated(true);
        static_cast<solvers::LimitingADERDGSolver*>(solver)->getSolver().get()->setNextMeshUpdateEvent(); // reset
        break;
      case solvers::Solver::Type::FiniteVolumes:
        // do nothing
        break;
    }
  }
  repository.switchToPredictionRerun(); // This one waits for background job termination
  repository.iterate(exahype::solvers::Solver::PredictionSweeps,false);
}

void exahype::runners::Runner::validateSolverTimeStepDataForThreeAlgorithmicPhases(const bool fuseADERDGPhases) const {
  #ifdef Asserts
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    assertionEquals(solver->getMinTimeStamp(),0.0);
    assertion1(std::isfinite(solver->getMinTimeStepSize()),solver->getMinTimeStepSize());
    assertion1(solver->getMinTimeStepSize()>0,solver->getMinTimeStepSize());

    switch(solver->getTimeStepping()) {
      case exahype::solvers::Solver::TimeStepping::Global:
        assertionEquals(solver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
        break;
      case exahype::solvers::Solver::TimeStepping::GlobalFixed:
        break;
    }
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG: {
        auto* aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(aderdgSolver->getMinCorrectorTimeStamp() > 0.0,aderdgSolver->getMinCorrectorTimeStamp());
        if (!fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

      }  break;
      case exahype::solvers::Solver::Type::LimitingADERDG: {
        // ADER-DG
        auto* aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(aderdgSolver->getMinCorrectorTimeStamp() > 0.0,aderdgSolver->getMinCorrectorTimeStamp());
        if (!fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

        // Finite Volumes
        auto* finiteVolumesSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter().get();
        assertionEquals(finiteVolumesSolver->getMinTimeStamp(),0.0);
        assertion1(std::isfinite(finiteVolumesSolver->getMinTimeStepSize()),finiteVolumesSolver->getMinTimeStepSize());
        assertion1(finiteVolumesSolver->getMinTimeStepSize()>0,finiteVolumesSolver->getMinTimeStepSize());
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(finiteVolumesSolver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      } break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        break;
    }
  }
  #endif
}
