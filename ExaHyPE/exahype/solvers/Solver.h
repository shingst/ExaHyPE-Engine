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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl
 **/
 
#ifndef _EXAHYPE_SOLVERS_SOLVER_H_
#define _EXAHYPE_SOLVERS_SOLVER_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <deque>

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "peano/utils/PeanoOptimisations.h"
#include "tarch/multicore/MulticoreDefinitions.h"
#include "tarch/la/Vector.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include "peano/utils/Globals.h"
#include "peano/grid/VertexEnumerator.h"
#include "peano/heap/DoubleHeap.h"
#include "peano/heap/CharHeap.h"
#include "peano/heap/HeapAllocator.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/State.h"

#include "exahype/profilers/Profiler.h"
#include "exahype/profilers/simple/NoOpProfiler.h"

#if defined(CompilerICC) && defined(SharedTBB)
// See: https://www.threadingbuildingblocks.org/tutorial-intel-tbb-scalable-memory-allocator
#include <tbb/cache_aligned_allocator.h> // prevents false sharing
#endif

// Some helpers
constexpr int power(int basis, int exp) {
  return (exp == 0) ? 1 : basis * power(basis, exp - 1);
}

#ifdef ALIGNMENT
constexpr int addPadding(const int originalSize) {
  return ALIGNMENT/8 * static_cast<int>((originalSize+(ALIGNMENT/8-1))/(ALIGNMENT/8));
}
#else
constexpr int addPadding(const int originalSize) {
  return originalSize;
}
#endif

namespace exahype {
  // Forward declarations
  class Cell;
  class Vertex;

  namespace parser {
  class ParserView;
  }  // namespace parser
  /**
   * We store the degrees of freedom associated with the ADERDGCellDescription and FiniteVolumesCellDescription
   * instances on this heap.
   * We further use this heap to send and receive face data from one MPI rank to the other.
   *
   * !!! CreateCopiesOfSentData
   *
   * All solvers must store the face data they send to neighbours at persistent addresses.
   */
  #ifdef ALIGNMENT
  #if defined(CompilerICC) && defined(SharedTBB)
  typedef tbb::cache_aligned_allocator<double> AlignedAllocator;
  #else
  typedef peano::heap::HeapAllocator<double, ALIGNMENT > AlignedAllocator;
  #endif
  typedef peano::heap::AlignedDoubleSendReceiveTask<ALIGNMENT> AlignedDoubleSendReceiveTask;
  typedef peano::heap::AlignedCharSendReceiveTask<ALIGNMENT>   AlignedCharSendReceiveTask;
  #endif

  #if defined(ALIGNMENT) and defined(UsePeanosSymmetricBoundaryExchanger)
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::SymmetricBoundaryDataExchanger< double, false, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    std::vector< double, AlignedAllocator >
  >     DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::SymmetricBoundaryDataExchanger< char, false, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    std::vector< char, AlignedAllocator >
  >     CompressedDataHeap;
  #elif defined(ALIGNMENT) and !defined(UsePeanosSymmetricBoundaryExchanger)
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::RLEBoundaryDataExchanger< double, false, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    std::vector< double, AlignedAllocator >
  >     DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::RLEBoundaryDataExchanger< char, false, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    std::vector< char, AlignedAllocator >
  >     CompressedDataHeap;
  #elif !defined(ALIGNMENT) and defined(UsePeanosSymmetricBoundaryExchanger)
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true,  peano::heap::SendReceiveTask<double> >,
    peano::heap::SynchronousDataExchanger< double, true,  peano::heap::SendReceiveTask<double> >,
    peano::heap::SymmetricBoundaryDataExchanger< double, false, peano::heap::SendReceiveTask<double> >
  >     DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SymmetricBoundaryDataExchanger< char, false, peano::heap::SendReceiveTask<char> >
  >     CompressedDataHeap;

  #elif defined(ALIGNMENT) and defined(UsePeanosAggregationBoundaryExchanger)
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::AggregationBoundaryDataExchanger< double, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    std::vector< double, AlignedAllocator >
  >     DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::AggregationBoundaryDataExchanger< char, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    std::vector< char, AlignedAllocator >
  >     CompressedDataHeap;
  #elif defined(ALIGNMENT) and !defined(UsePeanosAggregationBoundaryExchanger)
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< double, true, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
    peano::heap::AggregationBoundaryDataExchanger< double, AlignedDoubleSendReceiveTask, std::vector< double, AlignedAllocator > >,
   std::vector< double, AlignedAllocator >
  >      DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::SynchronousDataExchanger< char, true, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    peano::heap::AggregationBoundaryDataExchanger< char, AlignedCharSendReceiveTask, std::vector< char, AlignedAllocator > >,
    std::vector< char, AlignedAllocator >
  >     CompressedDataHeap;
  #elif !defined(ALIGNMENT) and defined(UsePeanosAggregationBoundaryExchanger)
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true,  peano::heap::SendReceiveTask<double> >,
    peano::heap::SynchronousDataExchanger< double, true,  peano::heap::SendReceiveTask<double> >,
    peano::heap::AggregationBoundaryDataExchanger< double, peano::heap::SendReceiveTask<double>, std::vector<double> >
  >     DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::AggregationBoundaryDataExchanger< char, peano::heap::SendReceiveTask<char>, std::vector<char> >
  >     CompressedDataHeap;
  #else
  typedef peano::heap::DoubleHeap<
    peano::heap::SynchronousDataExchanger< double, true,  peano::heap::SendReceiveTask<double> >,
    peano::heap::SynchronousDataExchanger< double, true,  peano::heap::SendReceiveTask<double> >,
    peano::heap::RLEBoundaryDataExchanger< double, false, peano::heap::SendReceiveTask<double> >
  >     DataHeap;
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::RLEBoundaryDataExchanger< char, false, peano::heap::SendReceiveTask<char> >
  >     CompressedDataHeap;
  #endif

  /**
   * @see waitUntilAllBackgroundTasksHaveTerminated()
   */
  extern tarch::multicore::BooleanSemaphore BackgroundJobSemaphore;

  /**
   * A semaphore for serialising heap access.
   */
  extern tarch::multicore::BooleanSemaphore HeapSemaphore;

#ifdef Parallel
  /**
   * An empty DataHeap message.
   *
   * !!! CreateCopiesOfSentData
   *
   * If we have set CreateCopiesOfSentData to
   * false for the DataHeap, all messages need to
   * have a fixed address as long as the send
   * process takes.
   *
   * Has to be declared extern in C++ standard as
   * it is instantiated in the corresponding cpp file.
   */
  extern DataHeap::HeapEntries EmptyDataHeapMessage;

  /**
   * We abuse this heap to send and receive metadata from one MPI rank to the other.
   * We never actually store data on this heap.
   *
   * !!! CreateCopiesOfSentData
   *
   * It is assumed by the metadata send routines of the solvers that
   * all data exchangers of the MetadataHeap create copies of the data to send.
   *
   * <h2> Implementation </h2>
   *
   * These meta data are not symmetric, i..e we can use the RLE heap but we
   * may not use any symmetric heap.
   */
  #if defined(UsePeanosAggregationBoundaryExchangerForMetaData)
  typedef peano::heap::CharHeap<
      peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
      peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
      peano::heap::AggregationBoundaryDataExchanger< char, peano::heap::SendReceiveTask<char>, std::vector<char> >
  >     MetadataHeap;
  #elif defined(UsePeanosSymmetricBoundaryExchangerForMetaData)
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SymmetricBoundaryDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >
  >     MetadataHeap;
  #else
  typedef peano::heap::CharHeap<
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::SynchronousDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >,
    peano::heap::RLEBoundaryDataExchanger< char, true,  peano::heap::SendReceiveTask<char> >
  >     MetadataHeap;
  #endif

  /**
   * Defines an invalid metadata entry.
   */
  static constexpr int InvalidMetadataEntry = -1;

  /**
   * Defines the length of the metadata
   * we send out per solver.
   *
   * First entry is the cell (description) type.
   * Second entry is the augmentation status,
   * third the helper status and the fourth the
   * limiter status.
   */
  static constexpr int NeighbourCommunicationMetadataPerSolver    = 4;

  static constexpr int NeighbourCommunicationMetadataCellType           = 0;
  static constexpr int NeighbourCommunicationMetadataAugmentationStatus = 1;
  static constexpr int NeighbourCommunicationMetadataCommunicationStatus       = 2;
  static constexpr int NeighbourCommunicationMetadataLimiterStatus      = 3;

  static constexpr int MasterWorkerCommunicationMetadataPerSolver       = 5;

  static constexpr int MasterWorkerCommunicationMetadataCellType           = 0;
  static constexpr int MasterWorkerCommunicationMetadataAugmentationStatus = 1;
  static constexpr int MasterWorkerCommunicationMetadataCommunicationStatus       = 2;
  static constexpr int MasterWorkerCommunicationMetadataLimiterStatus      = 3;
  static constexpr int MasterWorkerCommunicationMetadataSendReceiveData    = 4;

  /**
   * TODO(Dominic): Docu is outdated
   *
   * Encodes the metadata as integer sequence.
   *
   * The first element refers to the number of
   * ADERDGCellDescriptions associated with this cell (nADERG).
   * The next 2*nADERG elements store a pair of
   * solver number, and cell description type (encoded as int)
   * for each ADERDGCellDescription associated with this cell (description).
   *
   * The element 1+2*nADERDG refers to the number of
   * FiniteVolumesCellDescriptions associated with this cell (nFV).
   * The remaining 2*nFV elements store a pair of
   * solver number, and cell description type (encoded as int)
   * for each FiniteVolumesCellDescription associated with this cell
   * (description).
   */
  MetadataHeap::HeapEntries gatherNeighbourCommunicationMetadata(
      const int cellDescriptionsIndex,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest);

  /**
   * TODO(Dominic): Add docu.
   */
  MetadataHeap::HeapEntries gatherMasterWorkerCommunicationMetadata(const int cellDescriptionsIndex);

  /**
   * Send metadata to rank \p toRank.
   */
  void sendNeighbourCommunicationMetadata(
      const int                                   toRank,
      const int                                   cellDescriptionsIndex,
      const tarch::la::Vector<DIMENSIONS,int>&    src,
      const tarch::la::Vector<DIMENSIONS,int>&    dest,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);

  /**
   * Receive metadata to rank \p toRank.
   *
   * \return The index of the received metadata message
   * on the exahype::MetadataHeap.
   */
  int receiveNeighbourCommunicationMetadata(
      const int                                   fromRank,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);

  /**
   * Send metadata to rank \p toRank.
   */
  void sendMasterWorkerCommunicationMetadata(
      const int                                   toRank,
      const int                                   cellDescriptionsIndex,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);

  /**
   * Receive metadata to rank \p toRank.
   *
   * \return The index of the received metadata message
   * on the exahype::MetadataHeap.
   */
  int receiveMasterWorkerCommunicationMetadata(
      const int                                   fromRank,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);

  /**
   * Send a metadata sequence filled with InvalidMetadataEntry
   * to rank \p toRank.
   */
  void sendNeighbourCommunicationMetadataSequenceWithInvalidEntries(
      const int                                   toRank,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);

  /**
   * Send a metadata sequence filled with InvalidMetadataEntry
   * to rank \p toRank.
   */
  void sendMasterWorkerCommunicationMetadataSequenceWithInvalidEntries(
      const int                                   toRank,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);

  /**
   * Drop metadata sent by rank \p fromRank.
   */
  void dropMetadata(
      const int                                   fromRank,
      const peano::heap::MessageType&             messageType,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level);
  #endif

  namespace solvers {
    class Solver;

    typedef std::vector<Solver*> RegisteredSolversEntries;
    /**
     * All the registered solvers. Has to be declared extern in C++ standard as
     * it is instantiated in the corresponding cpp file.
     */
    extern std::vector<Solver*> RegisteredSolvers;

    /**
     * The limiter domain change that was detected after a solution
     * update or during the limiter status spreading.
     */
    enum class LimiterDomainChange {
      /**
       * A regular change of the limiter domain
       * has occurred. This might be either no change at
       * all or a situation where a cell directly next to a
       * troubled cell has been newly marked as troubled.
       */
      Regular,

      /**
       * A cell which is not directly next to a troubled cell
       * has newly been marked as troubled.
       * The cell and all its neighbours with Limiter-Status other than Ok,
       * are on the finest mesh level.
       *
       * <h2>Consequences</h2>
       * Usually, only a limiter status spreading, reinitialisation
       * and recomputation is necessary.
       *
       * However, if we further detect a cell of type Descendant/EmptyDescendant
       * marked with LimiterStatus other than Ok on the finest mesh level,
       * this status will change to IrregularChangeOnCoarserMeshLevel.
       */
      Irregular,

      /**
       * Scenario 1:
       * A cell which is not directly next to a troubled cell
       * has newly been marked as troubled.
       * The cell is not on the finest mesh level.
       *
       * Scenario 2:
       * A cell of type Descendant/EmptyDescendant
       * was marked with LimiterStatus other than Ok.
       * The cell is on the finest mesh level.
       *
       * <h2>Consequences</h2>
       * The runner must then refine the mesh accordingly, and perform a
       * rollback in all cells to the previous solution. It computes
       * a new time step size in all cells. Next, it recomputes the predictor in all
       * cells, troubled or not. Finally, it reruns the whole ADERDG time step in
       * all cells, troubled or not.
       *
       * This can potentially be relaxed for anarchic time stepping where
       * each cell has its own time step size and stamp.
       */
      IrregularRequiringMeshUpdate,
    };

    /**
     * Converts LimiterDomainChange to its double equivalent.
     */
    double convertToDouble(const LimiterDomainChange& limiterDomainChange);

    /**
     * Converts a double to a LimiterDomainChange.
     */
    LimiterDomainChange convertToLimiterDomainChange(const double value);
  }
}

/**
 * Describes one solver.
 */
class exahype::solvers::Solver {
 private:
  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  /**
   * A job that calls Solver::adjustSolutionDuringMeshRefinementBody(...).
   */
  class AdjustSolutionDuringMeshRefinementJob {
  private:
    Solver&       _solver;
    const int     _cellDescriptionsIndex;
    const int     _element;
    const bool    _isInitialMeshRefinement;
  public:
    AdjustSolutionDuringMeshRefinementJob(
        Solver&       solver,
        const int     cellDescriptionsIndex,
        const int     element,
        const bool    isInitialMeshRefinement);

    bool operator()();
  };

 protected:

  void tearApart(int numberOfEntries, int normalHeapIndex, int compressedHeapIndex, int bytesForMantissa) const;
  void glueTogether(int numberOfEntries, int normalHeapIndex, int compressedHeapIndex, int bytesForMantissa) const;

  /**
   * \see body adjustSolutionDuringMeshRefinement(...).
   * Must be implemented by the subclasses.
   */
  virtual void adjustSolutionDuringMeshRefinementBody(
      const int cellDescriptionsIndex,
      const int element,
      const bool isInitialMeshRefinement) = 0;

 public:

  /**
   * TrackGridStatistics is a flag from Peano that I "misuse" here as these
   * data also are grid statistics.
   */
  #ifdef TrackGridStatistics
  static double PipedUncompressedBytes;
  static double PipedCompressedBytes;
  #endif


  /**
   * A flag indicating we fuse the algorithmic
   * phases of all ADERDGSolver and
   * LimitingADERDGSolver instances.
   *
   * TODO(Dominic): Make private and hide in init function
   */
  static bool FuseADERDGPhases;

  /**
   * The weight which is used to scale
   * the stable time step size the fused
   * ADERDG time stepping scheme is
   * reset to after a rerun has become necessary.
   *
   * TODO(Dominic): Further consider to introduce
   * a second weight for the averaging:
   *
   * t_est = 0.5 (t_est_old + beta t_stable), beta<1.
   *
   * fuse-algorithmic-steps-reset-factor
   * fuse-algorithmic-steps-averaging-factor
   *
   * TODO(Dominic): Make private and hide in init function
   */
  static double WeightForPredictionRerun;

  /**
   * If this is set, we can skip sending metadata around during
   * batching iterations.
   */
  static bool DisableMetaDataExchangeInBatchedTimeSteps;

  /**
   * If this is set, we can skip Peano vertex neighbour exchange during batching iterations.
   */
  static bool DisablePeanoNeighbourExchangeInTimeSteps;

  /**
   * Set to 0 if no floating point compression is used. Is usually done in the
   * runner once at startup and from hereon is a read-only variable. The
   * subsequent field SpawnCompressionAsBackgroundThread has no semantics if
   * the present value is set to 0.
   */
  static double CompressionAccuracy;

  static bool SpawnCompressionAsBackgroundJob;

  /**
   * Set to true if the prediction and/or the fused time step
   * should be launched as background job whenever possible.
   *
   * TODO(Dominic): Rename as we start other background jobs as well
   */
  static bool SpawnPredictionAsBackgroundJob;

  /**
   * Set to true if the mesh refinement iterations
   * should run background jobs whenever possible.
   */
  static bool SpawnAMRBackgroundJobs;

  /**
   * The type of a solver.
   */
  enum class Type { ADERDG, FiniteVolumes, LimitingADERDG };

  /**
   * The time stepping mode.
   */
  enum class TimeStepping {
    /**
     * In the global time stepping mode, every cells works with the same time step.
     */
    Global,
    /**
     * In the fixed time stepping mode, we assume that each cell advanced in
     * time with the prescribed time step size. No CFL condition is checked.
     */
    GlobalFixed
    // Local, Anarchic
  };

  /**
   * The refinement control states
   * returned by the user functions.
   */
  enum class RefinementControl { Keep = 0, Refine = 1, Erase = 2 };

  /**
   * TODO(Dominic): Add docu.
   */
  typedef struct UpdateResult {
    double _timeStepSize                     = std::numeric_limits<double>::max();
    LimiterDomainChange _limiterDomainChange = LimiterDomainChange::Regular;
    bool _refinementRequested                = false;

    UpdateResult() {}
  } UpdateResult;

  /**
   * This struct is used in the AMR context
   * to lookup a parent cell description and
   * for computing the subcell position of the child
   * with respect to this parent.
   *
   * TODO(Dominic): Move to more appropriate place.
   */
  typedef struct SubcellPosition {
    int parentCellDescriptionsIndex;
    int parentElement;
    tarch::la::Vector<DIMENSIONS,int> subcellIndex;
    int levelDifference;

    SubcellPosition() :
      parentCellDescriptionsIndex(),
      parentElement(NotFound),
      subcellIndex(-1),
      levelDifference(-1) {}

    void invalidate() {
      parentCellDescriptionsIndex = multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex;
      parentElement = NotFound;
      for (int i=0; i<DIMENSIONS; i++) {
        subcellIndex[i] = -1;
      }
      levelDifference = -1;
    }

    ~SubcellPosition() {}
  } SubcellPosition;

  /**
   * The augmentation control states.
   */
  enum class AugmentationControl {
    /**
     * Indicates that a spacetree cell is next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Cell.
     */
    NextToCell = 0,
    /**
     * Indicates that a spacetree cell is next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor.
     */
    NextToAncestor = 1,
    /**
     * Indicates that a spacetree cell is both, next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor, and
     * a spacetree cell of type
     * exahype::records::ADERDGCellDescription::Cell.
     */
    NextToCellAndAncestor = 2,
    /**
     * Indicates that a spacetree cell is neither, next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor, nor
     * next to a spacetree cell of type exahype::records::ADERDGCellDescription::Cell.
     *
     * A cell of type exahype::records::ADERDGCellDescription::Descendant can then request erasing.
     * A cell of type exahype::records::ADERDGCellDescription::Cell does then not need
     * to request augmenting.
     */
     Default = 3
  };

  /**
   * Default return value of function getElement(...)
   * If we do not find the element in a vector
   * stored at a heap address.
   */
  static const int NotFound;

  /**
   * Moves a DataHeap array, i.e. copies the found
   * data at "fromIndex" to the array at "toIndex" and
   * deletes the "fromIndex" array afterwards.
   */
  static void moveDataHeapArray(const int fromIndex,const int toIndex,bool recycleFromArray);

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  static double getMinTimeStampOfAllSolvers();

  /**
   * Run over all solvers and identify the minimal sum of minimal time stamp
   * plus the minimal time step size.
   *
   * The result is a lower bound of the minimum time stamp
   * that will be obtained in the following time step.
   */
  static double estimateMinNextSolverTimeStampOfAllSolvers();

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  static double getMinTimeStepSizeOfAllSolvers();

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  static double getMaxSolverTimeStepSizeOfAllSolvers();

  /**
   * Run over all solvers and identify the maximal time stamp.
   *
   * On the individual patches, we do use the min time stamp so
   * far, so the routine returns the maximum over all min solver
   * time stamps.
   */
  static double getMaxTimeStampOfAllSolvers();

  static bool allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping scheme);

  static double getCoarsestMaximumMeshSizeOfAllSolvers();
  static double getFinestMaximumMeshSizeOfAllSolvers();

  /**
   * Returns the coarsest level which holds patches of
   * a solver.
   *
   * \note It is very important that initSolvers
   * has been called on all solvers before this
   * method is used.
   */
  static int getCoarsestMeshLevelOfAllSolvers();

  /**
   * Returns the finest level which holds a uniform
   * base grid of a solver.
   *
   * \note It is very important that initSolvers
   * has been called on all solvers before this
   * method is used.
   */
  static int getFinestUniformMeshLevelOfAllSolvers();

  /**
   * Run over all solvers and identify the maximum depth of adaptive
   * refinement employed.
   *
   * This number might correlate with the number
   * of grid iterations to run for performing an erasing operation.
   *
   * \note It is very important that initSolvers
   * has been called on all solvers before this
   * method is used.
   */
  static int getMaxAdaptiveRefinementDepthOfAllSolvers();

  /**
   * Loop over the solver registry and check if no solver
   * performs adaptive mesh refinement.
   */
  static bool allSolversPerformOnlyUniformRefinement();

  /**
   * Loop over the solver registry and check if one
   * of the solvers has requested a mesh update.
   */
  static bool oneSolverRequestedMeshUpdate();

  /**
   * Loop over the solver registry and check if one
   * of the solver's mesh refinement has not
   * attained a stable state yet.
   *
   * TODO(Dominic): Make this a state attribute since
   * we do not need to know which particular solver
   * did not attain a stable state.
   */
  static bool oneSolverHasNotAttainedStableState();

  /**
   * Returns true if one of the solvers used a time step size
   * that violated the CFL condition.
   *
   * TODO(Dominic): Rename. Name can be confused with
   * oneSolverHasNotAttainedStableState.
   */
  static bool oneSolverViolatedStabilityCondition();

  /**
   * Weights the min next predictor time step size
   * by the user's safety factor for the fused time stepping
   * algorithm.
   */
  static void weighMinNextPredictorTimeStepSize(exahype::solvers::Solver* solver);

  /**
   * Reinitialises the corrector and predictor time step sizes of an ADER-DG solver
   * with stable time step sizes if we detect a-posteriori that the CFL condition was
   * harmed by the estimated predictor time step size used in the last iteration.
   */
  static void reinitialiseTimeStepDataIfLastPredictorTimeStepSizeWasInstable(exahype::solvers::Solver* solver);

  /**
   * Starts a new time step on all solvers.
   *
   * \param[in] meshUpdateRequests   flags for each solver indicating if a mesh update is necessary.
   * \param[in] limiterDomainChanges flags for each solver indicating if the limiter domain has changed.
   * \param[in] minTimeStepSizes     the minimum CFL-stable time step size for all solvers.
   * \param[in] minCellSizes         the minimum cell size found in the grid for each solver.
   * \param[in] maxCellSizes         the maximum cell size found in the grid for each solver.
   * \param[in] isFirstIterationOfBatchOrNoBatch we run the first iteration of a batch or no batch at all
   * \param[in] isLastIterationOfBatchOrNoBatch we run the last iteration of a batch or no batch at all
   * \param[in] fusedTimeStepping fused time stepping is used or not
   */
  static void startNewTimeStepForAllSolvers(
      const std::vector<double>& minTimeStepSizes,
      const std::vector<int>& maxLevels,
      const std::vector<bool>& meshUpdateRequests,
      const std::vector<exahype::solvers::LimiterDomainChange>& limiterDomainChanges,
      const bool isFirstIterationOfBatchOrNoBatch,
      const bool isLastIterationOfBatchOrNoBatch,
      const bool fusedTimeStepping);

 /**
  * Ensure that all background jobs (such as prediction or compression jobs) have terminated before progressing
  * further. We have to wait until all tasks have terminated if we want to modify the heap,
  * i.e. insert new data or remove data.
  * Therefore, the wait (as well as the underlying semaphore) belong
  * into this abstract superclass.
  */
 static void ensureAllBackgroundJobsHaveTerminated();

 protected:
  /**
   * @see waitUntilAllBackgroundTasksHaveTerminated()
   */
  static int                                _NumberOfBackgroundJobs;

  /**
   * Each solver has an identifier/name. It is used for debug purposes only.
   */
  const std::string _identifier;

  const Type _type;

  /**
   * The number of state variables of the conservation or balance law.
   */
  const int _numberOfVariables;

  /**
   * The number of parameters, e.g, material parameters.
   */
  const int _numberOfParameters;

  /**
   * The number of nodal basis functions that are employed in each
   * coordinate direction.
   */
  const int _nodesPerCoordinateAxis;

  /**
   * The offset of the computational domain.
   *
   * Is initialised by the initSolver method.
   */
  tarch::la::Vector<DIMENSIONS,double> _domainOffset;

  /**
   * The size of the computational domain.
   * * Is initialised by the initSolver method.
   */
  tarch::la::Vector<DIMENSIONS,double> _domainSize;

  /**
   * The maximum extent a cell is allowed to have in each coordinate direction.
   *
   * \note This is an upper bound specified in the specification file.
   * This is not the actual maximum extent of a cell.
   */
  const double _maximumMeshSize;

  /**
   * The coarsest level of the adaptive mesh that is
   * occupied by this solver.
   */
  int _coarsestMeshLevel;

  /**
   * The maximum depth the adaptive mesh is allowed
   * to occupy (set by the user).
   * Summing this value with _coarsestMeshdLevel results in
   * the finest mesh level the solver might occupy during the
   * simulation.
   */
  const int _maximumAdaptiveMeshDepth;

  /**
   * The maximum tree level which is occupied by
   * cells of this solver.
   *
   * This value needs to be updated every time the grid has been changed.
   */
  int _maxLevel;
  int _nextMaxLevel;

  /**
   * The time stepping mode of this solver.
   */
  const TimeStepping _timeStepping;

  /**
   * A profiler for this solver.
   */
  std::unique_ptr<profilers::Profiler> _profiler;

  /**
   * Flag indicating if a mesh update was
   * requested by this solver.
   *
   * This is the state after the
   * time step size computation.
   *
   * <h2>MPI</h2>
   * This is the state after this rank's
   * solver has merged its state
   * with its workers' worker.
   */
  bool _meshUpdateRequest;

  /**
   * Flag indicating if a mesh update was
   * requested by this solver.
   *
   * This is the state before the
   * time step size computation.
   *
   * <h2>MPI</h2>
   * This is the state before this rank's
   * solver has merged its state
   * with its workers' solver.
   */
  bool _nextMeshUpdateRequest;

  /**
   * Flag indicating if the mesh refinement
   * performed by this solver attained a stable state.
   *
   * <h2>MPI</h2>
   * This is the state after this rank's
   * solver has merged its state
   * with its workers' worker.
   */
  bool _attainedStableState;

  /**
   * Flag indicating if the mesh refinement
   * performed by this solver attained a stable state.
   *
   * This is the state before the
   * time step size computation.
   *
   * <h2>MPI</h2>
   * This is the state before this rank's
   * solver has merged its state
   * with its workers' solver.
   */
  bool _nextAttainedStableState;

 public:
  Solver(const std::string& identifier, exahype::solvers::Solver::Type type,
         int numberOfVariables, int numberOfParameters,
         int nodesPerCoordinateAxis,
         double maximumMeshSize,
         int maximumAdaptiveMeshDepth,
         exahype::solvers::Solver::TimeStepping timeStepping,
         std::unique_ptr<profilers::Profiler> profiler =
             std::unique_ptr<profilers::Profiler>(
                 new profilers::simple::NoOpProfiler("")));

  virtual ~Solver() { _profiler->writeToConfiguredOutput(); }

  // Disallow copy and assignment
  Solver(const Solver& other) = delete;
  Solver& operator=(const Solver& other) = delete;

  /**
   * Return a string representation for the type \p param.
   */
  static std::string toString(const exahype::solvers::Solver::Type& param);

  /**
   * Return a string representation for the time stepping mode \p param.
   */
  static std::string toString(const exahype::solvers::Solver::TimeStepping& param);

  /**
   * Return the mesh level corresponding to the given mesh size with
   * respect to the given domainSize.
   *
   * \note That the domain root cell is actually at Peano mesh level 1
   * since the domain itself is embedded in a 3^d mesh in Peano.
   *
   * \note Load balancing makes only sense for a Peano mesh with
   * at least 3 (Peano) levels.
   * This is not ensured or checked in this routine.
   */
  static int computeMeshLevel(double meshSize, double domainSize);

  /**
   * Returns the maximum extent a mesh cell is allowed to have
   * in all coordinate directions.
   * This maximum mesh size is used both as a
   * constraint on the AMR as well as to set up the initial
   * grid. If you return the extent of the computational domain in
   * each coordinate direction or larger values,
   * you indicate that this solver is not active in the domain.
   */
  double getMaximumMeshSize() const;

  /**
   * The coarsest level of the adaptive mesh that is
   * occupied by this solver.
   */
  int getCoarsestMeshLevel() const;

  /**
   * The maximum depth the adaptive mesh is allowed to
   * occupy (set by the user).
   * Summing this value with _coarsestMeshdLevel results in
   * the finest mesh level the solver might occupy during the
   * simulation.
   */
  int getMaximumAdaptiveMeshDepth() const;

  /**
   * The finest level of the adaptive mesh that might be
   * occupied by this solver.
   */
  int getMaximumAdaptiveMeshLevel() const;

  /**
   * \note methods are virtual in order to enable
   * overriding by LimitingADERDGSolver
   */
  virtual void updateNextMaxLevel(int maxLevel);
  virtual int getNextMaxLevel() const;
  virtual int getMaxLevel() const;

  /**
   * Returns the identifier of this solver.
   */
  std::string getIdentifier() const;

  /**
   * Returns the type of this solver.
   */
  Type getType() const;

  /**
   * Returns the time stepping algorithm this solver is using.
   */
  TimeStepping getTimeStepping() const;

  /**
   * Returns the number of state variables.
   */
  int getNumberOfVariables() const;

  /**
   * Returns the number of parameters, e.g.,material constants etc.
   */
  int getNumberOfParameters() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

  virtual std::string toString() const;

  virtual void toString(std::ostream& out) const;

  /**
   * Reset the mesh update flags.
   *
   * \deprecated
   */
  virtual void resetMeshUpdateRequestFlags();

  /**
   * Update if a mesh update was requested by this solver.
   *
   * <h2>MPI</h2>
   * This is the state before we have send data to the master rank
   * and have merged the state with this rank's workers.
   */
  virtual void updateNextMeshUpdateRequest(const bool& meshUpdateRequest);

  /**
   * Indicates if a mesh update was requested
   * by this solver.
   *
   * <h2>MPI</h2>
   * This is the state before we have send data to the master rank
   * and have merged the state with this rank's workers.
   */
  virtual bool getNextMeshUpdateRequest() const;

  /**
   * Indicates if a mesh update was requested
   * by this solver.
   *
   * <h2>MPI</h2>
   *This is the state before we have send data to the master rank
   * and have merged the state with this rank's workers.
   */
  virtual bool getMeshUpdateRequest() const;

  /**
   * Overwrite the _MeshUpdateRequest flag
   * by the _nextMeshUpdateRequest flag.
   * Reset the _nextMeshUpdateRequest flag
   * to false;
   */
  virtual void setNextMeshUpdateRequest();

  /**
   * Update if the mesh refinement of this solver attained
   * a stable state.
   *
   * <h2>MPI</h2>
   * This is the state before we have send data to the master rank
   * and have merged the state with this rank's workers.
   */
  virtual void updateNextAttainedStableState(const bool& attainedStableState);

  /**
   * Indicates if the mesh refinement of this solver
   * attained a stable state.
   *
   * <h2>MPI</h2>
   * This is the state before we have send data to the master rank
   * and have merged the state with this rank's workers.
   */
  virtual bool getNextAttainedStableState() const;

  /**
   * Indicates if the mesh refinement of this solver
   * attained a stable state.
   *
   * <h2>MPI</h2>
   *This is the state before we have send data to the master rank
   * and have merged the state with this rank's workers.
   */
  virtual bool getAttainedStableState() const;

  /**
   * Overwrite the _attainedStableState flag
   * by the _nextAttainedStableState flag.
   * Reset the _nextAttainedStableState flag
   * to false;
   */
  virtual void setNextAttainedStableState();

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  virtual double getMinTimeStamp() const = 0;

  /**
   * The minimum time step size
   * of all cell descriptions.
   */
  virtual double getMinTimeStepSize() const = 0;

  virtual void updateMinNextTimeStepSize(double value) = 0;

  /**
   * Initialise the solver's time stamps and time step sizes.
   * Further use the bounding box and the already known
   * maximum mesh size to compute the coarsest grid level
   * this solver is placed on.
   *
   * \note It is very important that the domainSize
   * is chosen as an multiple of the coarsest mesh size
   * of all solvers within the grid.
   *
   * The maximum adaptive refinement level is defined
   * with respect to this level.
   */

  virtual void initSolver(
      const double timeStamp,
      const tarch::la::Vector<DIMENSIONS,double>& domainOffset,
      const tarch::la::Vector<DIMENSIONS,double>& domainSize,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
      const std::vector<std::string>& cmdlineargs,
      const exahype::parser::ParserView& parserView) = 0;

  /**
   * \return true if the solver is computing in the current algorithmic section.
   * This depends usually on internal flags of the solver such as ones indicating
   * a mesh update request or a limiter domain change during a previous time stepping
   * iteration.
   *
   * All mappings introduced for a specific job, e.g. limiting, mesh refinement etc.,
   * do not rely on this method. It is used only for mappings which are shared by different
   * algorithm sections. These are
   * BroadcastAndMergeTimeStepData, Merging, Prediction, and TimeStepSizeComputation,
   * MeshRefinement (+FinaliseMeshRefinement)
   *
   * E.g. a time step size computation via mapping TimeStepSizeComputation is required in
   * algorithm section "LocalRecomputationAllSend" for all solvers which
   * have finished a global recomputation or a mesh refinement.
   * It is not required for limiting ADER-DG solvers which are currently performing
   * a local recomputation.
   *
   */
  virtual bool isPerformingPrediction(const exahype::State::AlgorithmSection& section) const = 0;

  /**
   * \return true if this solver needs to merge metadata only(!) in the current algorithm section.
   */
  virtual bool isMergingMetadata(const exahype::State::AlgorithmSection& section) const = 0;

  /**
   * Copies the time stepping data from the global solver onto the patch's time
   * meta data.
   *
   * \param[in] element Index of the cell description in
   *                    the array at address \p cellDescriptionsIndex.
   */
  virtual void synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) const = 0;

  /**
   * This routine is called if we perform
   * a time step update
   */
  virtual void startNewTimeStep() = 0;

  /**
   * Similar as Solver::startNewTimeStep but
   * for the fused time stepping scheme.
   *
   * \param[in] isFirstIterationOfBatch indicates that we are in the first iteration
   *                                    of a batch or not. Note that this must be also set to true
   *                                    in case we run a batch of size 1, i.e. no batch at all.
   *
   * \param[in] isLastIterationOfBatch indicates that we are in the last iteration
   *                                    of a batch or not. Note that this must be also set to true
   *                                    in case we run a batch of size 1, i.e. no batch at all.
   */
  virtual void startNewTimeStepFused(
      const bool isFirstIterationOfBatch,
      const bool isLastIterationOfBatch) = 0;

  /**
   * Similar as Solver::updateTimeStepSizes but
   * for the fused time stepping scheme.
   */
  virtual void updateTimeStepSizesFused() = 0;

  /**
   * In contrast to startNewTimeStep(), this
   * method does not shift the time stamp.
   *
   * It simply updates the time step sizes
   *
   * This method is used after a mesh refinement.
   */
  virtual void updateTimeStepSizes() = 0;

  /**
   * Zeroes all the time step sizes.
   * This method is used by the adaptive mesh refinement mapping.
   * After the mesh refinement, we need to recompute
   * the time step sizes.
   *
   * <h1>ADER-DG</h1>
   * Further resets the predictor time stamp to take
   * the value of the corrector time stamp.
   * The fused must be initialised again after
   * each mesh refinement.
   *
   *   // TODO(Dominic): Still neccessary?
   */
  virtual void zeroTimeStepSizes() = 0;

  /**
   * Rolls the solver time step data back to the
   * previous time step.
   * Note that the newest time step
   * data is lost in this process.
   */
  virtual void rollbackToPreviousTimeStep() = 0;

  /**
   * Same as Solver::rollbackToPreviousTimeStep
   * for fused time stepping.
   */
  virtual void rollbackToPreviousTimeStepFused() = 0;

  virtual double getMinNextTimeStepSize() const=0;

  /**
   * Returns true if the index \p cellDescriptionsIndex
   * is a valid heap index.
   */
  virtual bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const = 0;

  /**
   * If an entry for this solver exists,
   * return the element index of the cell description
   * in the array at address \p cellDescriptionsIndex.
   * Otherwise and if \p cellDescriptionsIndex is an
   * invalid index, return Solver::NotFound.
   */
  virtual int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const = 0;

  /**
   * Modify a cell description in enter cell event.
   * This event should be used for single cell operations
   * like marking for refinement, erasing, augmenting,
   * or deaugmenting.
   *
   * \return a struct of type bool.
   *
   * \note We use this at the moment only
   * for refinement events. We can consider later
   * on to merge the time stepping functionality
   * (solution update, predictor comp.) into
   * this hook.
   */
  virtual bool progressMeshRefinementInEnterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const bool initialGrid,
      const int solverNumber) = 0;

  /**
   * Refinement routine that should be used for
   * collective children-parent operations.
   *
   * \return If a new compute cell was introduced
   * as part of a refinement operation.
   */
  virtual bool progressMeshRefinementInLeaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) = 0;

  /**
   * \return if the vertices around a cell should be erased, kept,
   * or refined.
   */
  virtual exahype::solvers::Solver::RefinementControl eraseOrRefineAdjacentVertices(
      const int& cellDescriptionsIndex,
      const int& solverNumber,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) const = 0;

  /**
   * Returns true if the solver has attained
   * a stable state on the cell description
   *
   * \param fineGridCell a fine grid cell
   * \param fineGridVertices vertices surrounding the fine grid cell
   * \param fineGridVerticesEnumerator a enumerator for the fine grid vertices
   * \param solverNumber a solver number
   */
  virtual bool attainedStableState(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const int solverNumber) const = 0;

  /**
   * This method is called after the
   * mesh refinement iterations where this
   * solver performs states updates in
   * the enterCell() and leaveCell().
   *
   * This method is used to finalise some state
   * updates.
   *
   * TODO(Dominic): Docu
   */
  virtual void finaliseStateUpdates(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) = 0;

  /////////////////////////////////////
  // CELL-LOCAL
  /////////////////////////////////////
  /**
   * Evaluate the refinement criterion after
   * a solution update has performed.
   *
   * We currently only return RefinementType::APrioriRefinement (or
   * RefinementType::APosterrioriRefinement)
   * if a cell requested refinement.
   * ExaHyPE might then stop the
   * time stepping and update the mesh
   * before continuing. Erasing is here not considered.
   *
   * \return True if mesh refinement is requested.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual bool evaluateRefinementCriterionAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Compute and return a new admissible time step
   * size for the cell description
   * \p element in the array at address
   * \p cellDescriptionsIndex.
   *
   * Then, update the time stamps and time step
   * sizes on the cell description accordingly.
   *
   * \return The new admissible time step size if the
   * cell description is of type Cell or
   * std::numeric_limits<double>::max().
   *
   * \note The update of the time stamps
   * and the time step sizes of the cell description
   * performed in this method can be revoked by the
   * solver's time stepping synchronisation mode.
   * If a time stepping mode like global or globalfixed
   * is chosen, then the changes made here are simply
   * overwritten. You might want to take a look
   * into method synchroniseTimeStepping().
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   * not copied once for each thread.
   *
   * \see startNewTimeStep(),
   *      synchroniseTimeStepping(int,int),
   *      synchroniseTimeStepping(CellDescription&)
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Same as \p startNewTimeStep for the fused time stepping scheme.
   *
   * \param[in] isFirstIterationOfBatch indicates that we are in the first iteration
   *                                    of a batch or not. Note that this must be also set to true
   *                                    in case we run a batch of size 1, i.e. no batch at all.
   *
   * \param[in] isLastIterationOfBatch indicates that we are in the last iteration
   *                                    of a batch or not. Note that this must be also set to true
   *                                    in case we run a batch of size 1, i.e. no batch at all.
   */
  virtual double startNewTimeStepFused(
        const int cellDescriptionsIndex,
        const int element,
        const bool isFirstIterationOfBatch,
        const bool isLastIterationOfBatch) = 0;

  /**
   * Computes a new time step size and overwrites
   * a cell description's time stamps and time step sizes
   * with it.
   *
   * In contrast to startNewTimeStep(int,int), this
   * method does not shift time stamps.
   *
   * This method is usually called after mesh refinement
   * was performed.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
    virtual double updateTimeStepSizes(
          const int cellDescriptionsIndex,
          const int element) = 0;

  /**
   * Same as ::updateTimeStepSizes for the fused
   * time stepping.
   */
  virtual double updateTimeStepSizesFused(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Zeroes all the time step sizes.
   * This method is used by the adaptive mesh refinement mapping.
   * After the mesh refinement, we need to recompute
   * the time step sizes.
   *
   * <h1>ADER-DG</h1>
   * Further resets the predictor time stamp to take
   * the value of the corrector time stamp.
   * The fused must be initialised again after
   * each mesh refinement.
   *
   * \note We do not overwrite _minNextTimeStepSize or an
   * equivalent value since this would erase the time
   * step size of the fixed time stepping schemes ("globalfixed" etc.)
   */
   // TODO(Dominic): Still neccessary?
  virtual void zeroTimeStepSizes(const int cellDescriptionsIndex, const int element) const = 0;

  /**
    * Rollback to the previous time step, i.e,
    * overwrite the time step size and time stamp
    * fields of the cell description
    * by previous values.
    */
   virtual void rollbackToPreviousTimeStep(
       const int cellDescriptionsIndex,
       const int solverElement) const = 0;

   /*
    * Same as LimitingADERDGSolver::rollbackToPreviousTimeStep
    * but for the fused time stepping scheme.
    */
   virtual void rollbackToPreviousTimeStepFused(
       const int cellDescriptionsIndex,
       const int solverElement) const = 0;

  /**
   * Impose initial conditions.
   *
   * \note Make sure to reset neighbour merge
   * helper variables in this method call.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  void adjustSolutionDuringMeshRefinement(
      const int cellDescriptionsIndex,
      const int element,
      const bool isInitialMeshRefinement);

  /**
   * Fuse algorithmic phases of the solvers.
   *
   * <h2>FiniteVolumesSolver</h2>
   *
   * This call degenerates to an updateSolution
   * call for the FiniteVolumesSolver.
   *
   * <h2>ADERDGSolver</h2>
   *
   * Runs the triad of updateSolution,performPredictionAndVolumeIntegral
   * plus startNewTimeStep.
   *
   * <h2>LimitingADERDGSolver</h2>
   *
   * Either runs the ADERDGSolver triad or
   * performs an FV update. Performs some additional
   * tasks.
   *
   * \param[in] isFirstIterationOfBatch Indicates that we currently run no batch or
   *                                    we are in the first iteration of a batch.
   * \param[in] isLastIterationOfBatch  Indicates that we currently run no batch or
   *                                    we are in the last iteration of a batch.
   *                                    (If no batch is run, both flags
   *                                    \p isFirstIterationOfBatch and
   *                                    \p isLastIterationOfBatch are true).
   * \param[in] isAtRemoteBoundary Flag indicating that the cell hosting the
   *                                    cell description is adjacent to a remote rank.
   */
  virtual UpdateResult fusedTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      const bool isFirstIterationOfBatch,
      const bool isLastIterationOfBatch,
      const bool isAtRemoteBoundary) = 0;

  /**
   * The nonfused update routine.
   *
   * <h2>FiniteVolumesSolver</h2>
   *
   * This call degenerates to an updateSolution
   * call for the FiniteVolumesSolver.
   *
   * <h2>ADERDGSolver</h2>
   *
   * Update the solution and evaluate the refinement criterion.
   *
   * <h2>LimitingADERDGSolver</h2>
   *
   * Update the ADER-DG and or FV solution and
   * evaluate the limiter and
   * the refinement criteria.
   *
   * \note Make sure to reset neighbour merge
   * helper variables in this method call.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual UpdateResult update(
          const int cellDescriptionsIndex,
          const int element,
          const bool isAtRemoteBoundary) = 0;

  /**
   * Explicitly ask the solver to compress
   * a cell description.
   *
   * \param[in] isAtRemoteBoundary Flag indicating that the cell hosting the
   *                               cell description is adjacent to a remote rank.
   */
  virtual void compress(
      const int cellDescriptionsIndex,
      const int element,
      const bool isAtRemoteBoundary) const = 0;

  /**
    * Prolongates face data from a parent cell description to
    * the cell description at address (cellDescriptionsIndex,element)
    * in case the fine grid cell associated with the cell description is adjacent to
    * the hull of the coarse grid cell associated with the parent cell description.
    *
    * Further zero out the face data of ancestors.
    *
    * \note This function assumes a top-down traversal of the grid and must thus
    * be called from the enterCell(...) mapping method.
    *
    * \note It is assumed that this operation is applied only to helper cell descriptions
    * of type Descendant and Ancestor. No cell description of type Cell
    * must be touched by this operation. Otherwise, we cannot spawn
    * the prediction and/or the compression as background task.
    *
    * \note Has no const modifier since kernels are not const functions yet.
    */
  virtual void prolongateAndPrepareRestriction(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Restrict data to a parent on
   * a coarser level.
   *
   * \note Thread-safety must be ensured by the implementation itself.
   *
   * \note This function assumes a bottom-up traversal of the grid and must thus
   * be called from the leaveCell(...) mapping method.
   *
   * \note Has no const modifier yet since kernels are not
   * const yet.
   */
  virtual void restriction(
      const int cellDescriptionsIndex,
      const int element) = 0;

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////

  /**
   * Merge the metadata of two cell descriptions.
   *
   * \param[in] element Index of the cell description
   *            holding the data to send out in
   *            the array at address \p cellDescriptionsIndex.
   *            This is not the solver number.
   *
   * \see tryGetElement
   */
  virtual void mergeNeighboursMetadata(
        const int                                 cellDescriptionsIndex1,
        const int                                 element1,
        const int                                 cellDescriptionsIndex2,
        const int                                 element2,
        const tarch::la::Vector<DIMENSIONS, int>& pos1,
        const tarch::la::Vector<DIMENSIONS, int>& pos2) const= 0;

  /**
   * Receive solver data from neighbour rank and write
   * it on the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   *  not copied once for each thread.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual void mergeNeighbours(
        const int                                 cellDescriptionsIndex1,
        const int                                 element1,
        const int                                 cellDescriptionsIndex2,
        const int                                 element2,
        const tarch::la::Vector<DIMENSIONS, int>& pos1,
        const tarch::la::Vector<DIMENSIONS, int>& pos2) = 0;

  /**
   * Take the cell descriptions \p element
   * from array at address \p cellDescriptionsIndex
   * and merge it with boundary data.
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   *  not copied once for each thread.
   *
   * \param[in] element Index of the cell description
   *                    at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual void mergeWithBoundaryData(
        const int                                 cellDescriptionsIndex,
        const int                                 element,
        const tarch::la::Vector<DIMENSIONS, int>& posCell,
        const tarch::la::Vector<DIMENSIONS, int>& posBoundary) = 0;

  #ifdef Parallel
  /**
   * If a cell description was allocated at heap address \p cellDescriptionsIndex
   * for solver \p solverNumber, encode metadata of the cell description
   * and push it to the back of the metadata vector \p metadata.
   *
   * Otherwise, push exahype::NeighbourCommunicationMetadataPerSolver
   * times exahype::InvalidMetadataEntry to the back of the vector.
   */
  virtual void appendNeighbourCommunicationMetadata(
      MetadataHeap::HeapEntries& metadata,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int cellDescriptionsIndex,
      const int solverNumber) const = 0;

  /**
   * Merge cell description \p element in
   * the cell descriptions array stored at \p
   * cellDescriptionsIndex with metadata.
   *
   * Currently, the neighbour metadata is only the neighbour
   * type as int \p neighbourTypeAsInt and
   * the neighbour's limiter status as int.
   *
   * <h2> Number of merges </h2>
   * Usually metadata is only merged once between neighbouring cells but
   * at a MPI boundary, it might be merged twice.
   * It is thus important that the inputs do not change after a merge.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithNeighbourMetadata(
      const MetadataHeap::HeapEntries&  metadata,
      const tarch::la::Vector<DIMENSIONS, int>& src,
      const tarch::la::Vector<DIMENSIONS, int>& dest,
      const int cellDescriptionsIndex,
      const int element)  const = 0;

  /**
   * Send solver data to neighbour rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual void sendDataToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to neighbour rank.
   */
  virtual void sendEmptyDataToNeighbour(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Receive solver data from neighbour rank and write
   * it on the cell description \p element in
   * the cell descriptions array stored at \p
   * cellDescriptionsIndex.
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   *  not copied once for each thread.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from neighbour rank.
   */
  virtual void dropNeighbourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  ///////////////////////////////////
  // WORKER<=>MASTER
  ///////////////////////////////////
  /**
   * Finishes outstanding refinement operations
   * and sends solution data down to the worker
   * if required.
   */
  virtual void progressMeshRefinementInPrepareSendToWorker(
      const int workerRank,
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const bool initialGrid,
      const int solverNumber) = 0;

  /**
   * Just receive data or not from the master
   * depending on the refinement event.
   */
  virtual void progressMeshRefinementInReceiveDataFromMaster(
      const int masterRank,
      const int receivedCellDescriptionsIndex,
      const int receivedElement,
      const peano::grid::VertexEnumerator& receivedVerticesEnumerator) const = 0;

  /**
   * Finish prolongation operations started on the master.
   */
  virtual void progressMeshRefinementInMergeWithWorker(
      const int localCellDescriptionsIndex,    const int localElement,
      const int receivedCellDescriptionsIndex, const int receivedElement,
      const bool initialGrid) = 0;

  /**
   * Finish erasing operations on the worker side and
   * send data up to the master if necessary.
   * This data is then picked up to finish restriction
   * operations.
   */
  virtual void progressMeshRefinementInPrepareSendToMaster(
        const int masterRank,
        const int cellDescriptionsIndex, const int element,
        const tarch::la::Vector<DIMENSIONS,double>& x,
        const int level) const override;

  /**
   * Finish erasing operations started on the master which
   * require data from the worker.
   *
   * Veto erasing requests from the coarse grid cell as well.
   */
  virtual bool progressMeshRefinementInMergeWithMaster(
      const int worker,
      const int localCellDescriptionsIndex,    const int localElement,
      const int receivedCellDescriptionsIndex, const int receivedElement,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * If a cell description was allocated at heap address \p cellDescriptionsIndex
   * for solver \p solverNumber, encode metadata of the cell description
   * and push it to the back of the metadata vector \p metadata.
   *
   * Otherwise, push exahype::MasterWorkerCommunicationMetadataPerSolver
   * times exahype::InvalidMetadataEntry to the back of the vector.
   *
   */
  virtual void appendMasterWorkerCommunicationMetadata(
      MetadataHeap::HeapEntries& metadata,
      const int cellDescriptionsIndex,
      const int solverNumber) const = 0;

  /**
   * Send solver data to master or worker rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Send empty solver data to master or worker rank
   * due to fork or join.
   */
  virtual void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Merge with solver data from master or worker rank
   * that was sent out due to a fork or join. Wrote the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Drop solver data from master or worker rank
   * that was sent out due to a fork or join.
   */
  virtual void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////

  /**
   * Determine if the solver has to reduce data for the cell description
   * \p element at heap address \p cellDescriptionsIndex.
   * This is the case if we have adaptive mesh refinement
   * activated and a cell description located at the
   * master worker boundary is of type Ancestor (not EmptyAncestor).
   * In this case the worker will restrict face data up to the
   * master in this iteration.
   *
   * Note that this function must be called in
   * prepareSendToWorker(...) if you plug it into the Peano engine.
   *
   * Note that this function is not in control of determining when to reduce
   * time step data. This should be done outside of this function.
   */
  virtual bool hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) const = 0;

  /**
   * Send data to the master that is not
   * depending on a particular cell description.
   *
   * This operation might be used for the reduction of a global
   * minimum time step size over all MPI ranks.
   *
   * \note We always assume that
   * startNewTimeStep() has been already called on the
   * local solver instance. You thus
   * have to return the updated local time step size.
   *
   * \see startNewTimeStep(), mergeWithWorkerData()
   */
  virtual void sendDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Merge with solver data from worker rank.
   */
  virtual void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Compile a message containing mesh update flags
   * for the master.
   *
   * The initial capacity defaults to 2 but can be modified
   * to attend more data to the message.
   *
   * \see exahype::solvers::Solver::sendMeshUpdateFlagsToMaster,
   * exahype::solvers::LimitingADERDGSolver::sendMeshUpdateFlagsToMaster
   */
  exahype::DataHeap::HeapEntries
  compileMeshUpdateFlagsForMaster(const int capacity=2) const;

  /*
   * Send the rank-local mesh update request and
   * limiter domain change to the master.
   *
   * At the time of sending data to the master,
   * we have already set the next
   * mesh update request locally.
   * We thus need to communicate the
   * current mesh update request to the master.
   */
  void sendMeshUpdateFlagsToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Merge with the workers mesh update flags.
   *
   * \see exahype::solvers::Solver::mergeWithWorkerMeshUpdateFlags,
   * exahype::solvers::LimitingADERDGSolver::mergeWithWorkerMeshUpdateFlags
   */
  void mergeWithWorkerMeshUpdateFlags(const DataHeap::HeapEntries& message);

  /**
   * Merge with the workers mesh update flags.
   *
   * The master has not yet performed swapped
   * the current values with the "next" values.
   * This will happen after the merge.
   */
  void mergeWithWorkerMeshUpdateFlags(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Send solver data to master rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToMaster(
      const int                                    masterRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Send empty solver data to master rank.
   */
  virtual void sendEmptyDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Merge with solver data from worker rank.
   * Write the data to the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual void mergeWithWorkerData(
      const int                                    workerRank,
      const MetadataHeap::HeapEntries&             workerMetadata,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from worker rank.
   */
  virtual void dropWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////
  /**
   * Send data to the worker that is not
   * depending on a particular cell description.
   *
   * This operation might be used for the synchronisation
   * of a global minimum time step size over all MPI ranks.
   *
   * \see startNewTimeStep(), mergeWithMasterData()
   */
  virtual void sendDataToWorker(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Merge with solver data from master rank.
   */
  virtual void mergeWithMasterData(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send solver data to worker rank. Read the data from
   * the cell description \p element in the cell descriptions
   * vector stored at \p cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \note Has no const modifier since kernels are not const functions yet.
   */
  virtual void sendDataToWorker(
      const int                                    workerRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to worker rank.
   */
  virtual void sendEmptyDataToWorker(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const = 0;

  /**
   * Receive solver data from the master rank.
   *
   * \param[inout] heapIndices a queue where the solver
   *                    needs to add DataHeap indices
   *                    of received data.
   */
  virtual void receiveDataFromMaster(
        const int                                    masterRank,
        std::deque<int>&                             receivedDataHeapIndices,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const = 0;

  /**
   * Pop the heap indices from the double ended
   * queue and merge the data - or not.
   * Delete the corresponding heap arrays.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithMasterData(
      const MetadataHeap::HeapEntries&             masterMetadata,
      std::deque<int>&                             receivedDataHeapIndices,
      const int                                    cellDescriptionsIndex,
      const int                                    element) const = 0;

  /**
   * Drop solver data from master rank.
   *
   * Just pop the heap indices from the double ended
   * queue and delete the corresponding heap arrays.
   */
  virtual void dropMasterData(
      std::deque<int>& heapIndices) const = 0;
  #endif
};

#endif
