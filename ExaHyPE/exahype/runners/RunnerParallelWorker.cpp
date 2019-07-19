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

#ifdef Parallel
#include "exahype/repositories/Repository.h"
#include "exahype/mappings/FinaliseMeshRefinement.h"
#include "peano/parallel/messages/ForkMessage.h"
#include "peano/utils/Globals.h"
#include "peano/utils/UserInterface.h"
#include "peano/utils/PeanoOptimisations.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/parallel/NodePool.h"

#include "peano/parallel/SendReceiveBufferPool.h"

#include "exahype/Vertex.h"

#if defined(DistributedOffloading)
#include "exahype/offloading/StaticDistributor.h"
#include "exahype/offloading/AggressiveDistributor.h"
#include "exahype/offloading/AggressiveCCPDistributor.h"
#include "exahype/offloading/AggressiveHybridDistributor.h"
#include "exahype/offloading/OffloadingManager.h"
#endif

#if defined(SharedMemoryParallelisation) && defined(PerformanceAnalysis)
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/datatraversal/autotuning/Oracle.h"
#endif

int exahype::runners::Runner::runAsWorker(
    exahype::repositories::Repository& repository) {
  int newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob();
  while (newMasterNode !=
         tarch::parallel::NodePool::JobRequestMessageAnswerValues::Terminate) {
    if (newMasterNode >=
        tarch::parallel::NodePool::JobRequestMessageAnswerValues::NewMaster) {
      peano::parallel::SendReceiveBufferPool::getInstance().createBufferManually<exahype::Vertex>(
          newMasterNode,peano::parallel::SendReceiveBufferPool::BufferAccessType::LIFO); // TODO(Dominic): LIFO or FIFO?

      peano::parallel::messages::ForkMessage forkMessage;
      forkMessage.receive(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          tarch::parallel::NodePool::getInstance().getTagForForkMessages(),
          true, peano::parallel::messages::ForkMessage::ExchangeMode::NonblockingWithPollingLoopOverTests);

      repository.restart(
          forkMessage.getH(), forkMessage.getDomainOffset(),
          forkMessage.getLevel(),
          forkMessage.getPositionOfFineGridCellRelativeToCoarseGridCell());

      bool continueToIterate = true;
      while (continueToIterate) {
        switch (repository.continueToIterate()) {
          case exahype::repositories::Repository::Continue:
            {
             preProcessTimeStepInSharedMemoryEnvironment();

             repository.iterate();

             static int lastMemoryUsageValue = 0;
             int memoryUsageDelta = 0;
             if (lastMemoryUsageValue > 0) {
               memoryUsageDelta = peano::utils::UserInterface::getMemoryUsageMB() - lastMemoryUsageValue;
             }

             lastMemoryUsageValue = peano::utils::UserInterface::getMemoryUsageMB();
              logInfo("runAsWorker(...)",
                "\tmemoryUsage    =" << lastMemoryUsageValue << " MB\t" << "(delta: " << memoryUsageDelta << ")" );

//              logInfo("runAsWorker(...)",
//                "\tmemoryDelta    =" << memoryUsageDelta << " MB");


              #if  defined(SharedMemoryParallelisation) && defined(PerformanceAnalysis)
              if (sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::hasLearnedSinceLastQuery()) {
                static int dumpCounter = -1;
                dumpCounter++;
                peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics( _parser.getMulticorePropertiesFile() + "-dump-" + std::to_string(dumpCounter) );
              }
              #endif

              postProcessTimeStepInSharedMemoryEnvironment();
            }
            break;
          case exahype::repositories::Repository::Terminate:
            continueToIterate = false;
            break;
          case exahype::repositories::Repository::RunGlobalStep:
            runGlobalStep();
            break;
        }
      }

      // insert your postprocessing here
      // -------------------------------

      // -------------------------------

      //repository.logIterationStatistics(false);
      repository.terminate();
    } else if (newMasterNode ==
               tarch::parallel::NodePool::JobRequestMessageAnswerValues::
                   RunAllNodes) {
      runGlobalStep();
    }
    newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob();
  }
  return 0;
}

void exahype::runners::Runner::runGlobalStep() {
  // You might want to remove this assertion, but please consult the
  // documentation before.
  assertion(!peano::parallel::loadbalancing::Oracle::getInstance()
                 .isLoadBalancingActivated());

#if defined(DistributedOffloading)
  // For the static offloading strategy, we need to allgather load information
  // from all ranks once and compute a new target load distribution.
#if defined(OffloadingStrategyStatic)
  logInfo("runner(...)",
          "running global step "<<exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells<<
		  ", "<<exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells );
  exahype::offloading::StaticDistributor::getInstance().computeNewLoadDistribution(
      exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells,
	  exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells);
#elif defined(OffloadingStrategyAggressive)
  logInfo("runner(...)",
          "running global step "<<exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells<<
		  ", "<<exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells );
  exahype::offloading::AggressiveDistributor::getInstance().computeIdealLoadDistribution(
      exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells,
	  exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells);
#elif defined(OffloadingStrategyAggressiveCCP)
  logInfo("runner(...)",
          "running global step "<<exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells<<
      ", "<<exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells );
  exahype::offloading::AggressiveCCPDistributor::getInstance().computeIdealLoadDistribution(
      exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells,
    exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells);
#elif defined(OffloadingStrategyAggressiveHybrid)
  logInfo("runner(...)",
          "running global step "<<exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells<<
      ", "<<exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells );
  exahype::offloading::AggressiveHybridDistributor::getInstance().computeIdealLoadDistribution(
      exahype::mappings::FinaliseMeshRefinement::NumberOfEnclaveCells,
    exahype::mappings::FinaliseMeshRefinement::NumberOfSkeletonCells);
#endif 


#if defined(OffloadingStrategyAggressive)
  exahype::offloading::AggressiveDistributor::getInstance().enable();
#elif defined(OffloadingStrategyAggressiveCCP)
  exahype::offloading::AggressiveCCPDistributor::getInstance().enable();
#elif defined(OffloadingStrategyAggressiveHybrid)
  exahype::offloading::AggressiveHybridDistributor::getInstance().enable();
#endif

#endif

}
#endif
