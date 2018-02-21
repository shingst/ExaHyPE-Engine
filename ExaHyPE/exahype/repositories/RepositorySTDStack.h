// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 


#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellSTDStack.h"
#include "peano/stacks/VertexSTDStack.h"


 #include "exahype/adapters/MeshRefinement.h" 
 #include "exahype/adapters/MeshRefinementAndPlotTree.h" 
 #include "exahype/adapters/FinaliseMeshRefinement.h" 
 #include "exahype/adapters/FinaliseMeshRefinementOrLocalRollback.h" 
 #include "exahype/adapters/FusedTimeStep.h" 
 #include "exahype/adapters/PredictionRerun.h" 
 #include "exahype/adapters/BroadcastAndDropNeighbourMessages.h" 
 #include "exahype/adapters/LimiterStatusSpreading.h" 
 #include "exahype/adapters/PredictionOrLocalRecomputation.h" 
 #include "exahype/adapters/GlobalRollback.h" 
 #include "exahype/adapters/BroadcastAndMergeNeighbours.h" 
 #include "exahype/adapters/UpdateAndReduce.h" 
 #include "exahype/adapters/Prediction.h" 



namespace exahype {
      namespace repositories {
        class RepositorySTDStack;  
      }
}


class exahype::repositories::RepositorySTDStack: public exahype::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellSTDStack<exahype::Cell>       CellStack;
    typedef peano::stacks::VertexSTDStack<exahype::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    exahype::State          _solverState;
    peano::grid::RegularGridContainer<exahype::Vertex,exahype::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MeshRefinement> _gridWithMeshRefinement;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MeshRefinementAndPlotTree> _gridWithMeshRefinementAndPlotTree;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FinaliseMeshRefinement> _gridWithFinaliseMeshRefinement;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FinaliseMeshRefinementOrLocalRollback> _gridWithFinaliseMeshRefinementOrLocalRollback;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FusedTimeStep> _gridWithFusedTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionRerun> _gridWithPredictionRerun;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::BroadcastAndDropNeighbourMessages> _gridWithBroadcastAndDropNeighbourMessages;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LimiterStatusSpreading> _gridWithLimiterStatusSpreading;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionOrLocalRecomputation> _gridWithPredictionOrLocalRecomputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GlobalRollback> _gridWithGlobalRollback;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::BroadcastAndMergeNeighbours> _gridWithBroadcastAndMergeNeighbours;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::UpdateAndReduce> _gridWithUpdateAndReduce;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Prediction> _gridWithPrediction;

     
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureMeshRefinementCPUTime;
    tarch::timing::Measurement _measureMeshRefinementAndPlotTreeCPUTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementCPUTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementOrLocalRollbackCPUTime;
    tarch::timing::Measurement _measureFusedTimeStepCPUTime;
    tarch::timing::Measurement _measurePredictionRerunCPUTime;
    tarch::timing::Measurement _measureBroadcastAndDropNeighbourMessagesCPUTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCPUTime;
    tarch::timing::Measurement _measurePredictionOrLocalRecomputationCPUTime;
    tarch::timing::Measurement _measureGlobalRollbackCPUTime;
    tarch::timing::Measurement _measureBroadcastAndMergeNeighboursCPUTime;
    tarch::timing::Measurement _measureUpdateAndReduceCPUTime;
    tarch::timing::Measurement _measurePredictionCPUTime;

    tarch::timing::Measurement _measureMeshRefinementCalendarTime;
    tarch::timing::Measurement _measureMeshRefinementAndPlotTreeCalendarTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementCalendarTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime;
    tarch::timing::Measurement _measureFusedTimeStepCalendarTime;
    tarch::timing::Measurement _measurePredictionRerunCalendarTime;
    tarch::timing::Measurement _measureBroadcastAndDropNeighbourMessagesCalendarTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCalendarTime;
    tarch::timing::Measurement _measurePredictionOrLocalRecomputationCalendarTime;
    tarch::timing::Measurement _measureGlobalRollbackCalendarTime;
    tarch::timing::Measurement _measureBroadcastAndMergeNeighboursCalendarTime;
    tarch::timing::Measurement _measureUpdateAndReduceCalendarTime;
    tarch::timing::Measurement _measurePredictionCalendarTime;

   
  public:
    RepositorySTDStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositorySTDStack(
      peano::geometry::Geometry&                   geometry
    );
    
    virtual ~RepositorySTDStack();

    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel,
      const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
    );
         
    virtual void terminate();
        
    virtual exahype::State& getState();
    virtual const exahype::State& getState() const;
	
    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);

    virtual void writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* createEmptyCheckpoint(); 

    virtual void switchToMeshRefinement();    
    virtual void switchToMeshRefinementAndPlotTree();    
    virtual void switchToFinaliseMeshRefinement();    
    virtual void switchToFinaliseMeshRefinementOrLocalRollback();    
    virtual void switchToFusedTimeStep();    
    virtual void switchToPredictionRerun();    
    virtual void switchToBroadcastAndDropNeighbourMessages();    
    virtual void switchToLimiterStatusSpreading();    
    virtual void switchToPredictionOrLocalRecomputation();    
    virtual void switchToGlobalRollback();    
    virtual void switchToBroadcastAndMergeNeighbours();    
    virtual void switchToUpdateAndReduce();    
    virtual void switchToPrediction();    

    virtual bool isActiveAdapterMeshRefinement() const;
    virtual bool isActiveAdapterMeshRefinementAndPlotTree() const;
    virtual bool isActiveAdapterFinaliseMeshRefinement() const;
    virtual bool isActiveAdapterFinaliseMeshRefinementOrLocalRollback() const;
    virtual bool isActiveAdapterFusedTimeStep() const;
    virtual bool isActiveAdapterPredictionRerun() const;
    virtual bool isActiveAdapterBroadcastAndDropNeighbourMessages() const;
    virtual bool isActiveAdapterLimiterStatusSpreading() const;
    virtual bool isActiveAdapterPredictionOrLocalRecomputation() const;
    virtual bool isActiveAdapterGlobalRollback() const;
    virtual bool isActiveAdapterBroadcastAndMergeNeighbours() const;
    virtual bool isActiveAdapterUpdateAndReduce() const;
    virtual bool isActiveAdapterPrediction() const;

   
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics(bool logAllAdapters) const;
    virtual void clearIterationStatistics();
};


#endif
