#include "exahype/adapters/MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.h"

#include <sstream>

#include "peano/utils/Loop.h"
#include "peano/grid/CellFlags.h"


#include "multiscalelinkedcell/HangingVertexBookkeeper.h"
#include "exahype/VertexOperations.h"


peano::CommunicationSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::communicationSpecification() const {
  return peano::CommunicationSpecification::getMinimalSpecification();
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::touchVertexFirstTimeSpecification(int level) const { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::enterCellSpecification(int level) const {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::ascendSpecification(int level) const {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::descendSpecification(int level) const {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6() {
}


exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::~MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6(const MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6&  masterThread) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithWorkerThread(const MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6& workerThread) {
}
#endif


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::createHangingVertex(
  exahype::Vertex&     fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
  exahype::Vertex * const   coarseGridVertices,
  const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
  exahype::Cell&       coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  const int level = coarseGridVerticesEnumerator.getLevel()+1;
  
  VertexOperations::writeCellDescriptionsIndex(
    fineGridVertex,
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createHangingVertex(
      fineGridX,level,
      fineGridPositionOfVertex,
      VertexOperations::readCellDescriptionsIndex(coarseGridVerticesEnumerator,coarseGridVertices))
    );
}



void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::destroyHangingVertex(
  const exahype::Vertex&   fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::createInnerVertex(
  exahype::Vertex&               fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  VertexOperations::writeCellDescriptionsIndex(
      fineGridVertex,
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForNewVertex());
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::createBoundaryVertex(
  exahype::Vertex&               fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  VertexOperations::writeCellDescriptionsIndex(
      fineGridVertex,
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForBoundaryVertex());
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::createCell(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::destroyCell(
  const exahype::Cell&           fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().destroyCell(fineGridCell.getCellDescriptionsIndex());
}


#ifdef Parallel
void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  VertexOperations::writeCellDescriptionsIndex(
    vertex,
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(
      vertex.getAdjacentRanks(),
      VertexOperations::readCellDescriptionsIndex(vertex)
    )
  );
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


bool exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithMaster(
  const exahype::Cell&           workerGridCell,
  exahype::Vertex * const        workerGridVertices,
  const peano::grid::VertexEnumerator& workerEnumerator,
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
  const exahype::State&          workerState,
  exahype::State&                masterState
) {
  dfor2(k)
    VertexOperations::writeCellDescriptionsIndex(
      fineGridVertices[ fineGridVerticesEnumerator(k) ],
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(
        fineGridVertices[ fineGridVerticesEnumerator(k) ].getAdjacentRanks(),
        VertexOperations::readCellDescriptionsIndex(fineGridVertices[ fineGridVerticesEnumerator(k) ])
      )
    );
  enddforx
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::receiveDataFromMaster(
      exahype::Cell&                        receivedCell, 
      exahype::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      exahype::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      exahype::Cell&                        receivedCoarseGridCell,
      exahype::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      exahype::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
  VertexOperations::writeCellDescriptionsIndex(
    localVertex,
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(
      localVertex.getAdjacentRanks(),
      VertexOperations::readCellDescriptionsIndex(localVertex)
    )
  );
}
#endif


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::enterCell(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  dfor2(k)
    if (fineGridVertices[fineGridVerticesEnumerator(k)].isHangingNode()) {
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().getAdjacencyEntriesOfVertex( 
        fineGridVerticesEnumerator.getVertexPosition(k),
        fineGridVerticesEnumerator.getLevel()
      )(TWO_POWER_D-kScalar-1) = fineGridCell.getCellDescriptionsIndex();
    }
    else {
      VertexOperations::writeCellDescriptionsIndex(
          fineGridVertices[fineGridVerticesEnumerator(k)], TWO_POWER_D-kScalar-1, fineGridCell.getCellDescriptionsIndex());
    }
  enddforx
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::beginIteration(
  exahype::State&  solverState
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().beginIteration();
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::endIteration(
  exahype::State&  solverState
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().endIteration();
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
}


void exahype::adapters::MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
}
