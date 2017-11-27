#include "exahype/adapters/MeshRefinementAndPlotGrid.h"


peano::CommunicationSpecification   exahype::adapters::MeshRefinementAndPlotGrid::communicationSpecification() const {
  return peano::CommunicationSpecification::getMinimalSpecification()
    &  _map2PreProcessing.communicationSpecification()
    &  _map2MeshRefinement.communicationSpecification()
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.communicationSpecification()
    &  _map2AugmentedAMRTreePlot2d.communicationSpecification()
    &  _map2LoadBalancing.communicationSpecification()
    &  _map2PostProcessing.communicationSpecification()
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.communicationSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2PreProcessing.touchVertexLastTimeSpecification(level)
    &  _map2MeshRefinement.touchVertexLastTimeSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.touchVertexLastTimeSpecification(level)
    &  _map2AugmentedAMRTreePlot2d.touchVertexLastTimeSpecification(level)
    &  _map2LoadBalancing.touchVertexLastTimeSpecification(level)
    &  _map2PostProcessing.touchVertexLastTimeSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.touchVertexLastTimeSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid::touchVertexFirstTimeSpecification(int level) const { 
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2PreProcessing.touchVertexFirstTimeSpecification(level)
    &  _map2MeshRefinement.touchVertexFirstTimeSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.touchVertexFirstTimeSpecification(level)
    &  _map2AugmentedAMRTreePlot2d.touchVertexFirstTimeSpecification(level)
    &  _map2LoadBalancing.touchVertexFirstTimeSpecification(level)
    &  _map2PostProcessing.touchVertexFirstTimeSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.touchVertexFirstTimeSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid::enterCellSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2PreProcessing.enterCellSpecification(level)
    &  _map2MeshRefinement.enterCellSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.enterCellSpecification(level)
    &  _map2AugmentedAMRTreePlot2d.enterCellSpecification(level)
    &  _map2LoadBalancing.enterCellSpecification(level)
    &  _map2PostProcessing.enterCellSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.enterCellSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid::leaveCellSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2PreProcessing.leaveCellSpecification(level)
    &  _map2MeshRefinement.leaveCellSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.leaveCellSpecification(level)
    &  _map2AugmentedAMRTreePlot2d.leaveCellSpecification(level)
    &  _map2LoadBalancing.leaveCellSpecification(level)
    &  _map2PostProcessing.leaveCellSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.leaveCellSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid::ascendSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2PreProcessing.ascendSpecification(level)
    &  _map2MeshRefinement.ascendSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.ascendSpecification(level)
    &  _map2AugmentedAMRTreePlot2d.ascendSpecification(level)
    &  _map2LoadBalancing.ascendSpecification(level)
    &  _map2PostProcessing.ascendSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.ascendSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::MeshRefinementAndPlotGrid::descendSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2PreProcessing.descendSpecification(level)
    &  _map2MeshRefinement.descendSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.descendSpecification(level)
    &  _map2AugmentedAMRTreePlot2d.descendSpecification(level)
    &  _map2LoadBalancing.descendSpecification(level)
    &  _map2PostProcessing.descendSpecification(level)
    &  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.descendSpecification(level)

  ;
}


exahype::adapters::MeshRefinementAndPlotGrid::MeshRefinementAndPlotGrid() {
}


exahype::adapters::MeshRefinementAndPlotGrid::~MeshRefinementAndPlotGrid() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::MeshRefinementAndPlotGrid::MeshRefinementAndPlotGrid(const MeshRefinementAndPlotGrid&  masterThread):
  _map2PreProcessing(masterThread._map2PreProcessing) , 
  _map2MeshRefinement(masterThread._map2MeshRefinement) , 
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2(masterThread._map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2) , 
  _map2AugmentedAMRTreePlot2d(masterThread._map2AugmentedAMRTreePlot2d) , 
  _map2LoadBalancing(masterThread._map2LoadBalancing) , 
  _map2PostProcessing(masterThread._map2PostProcessing) , 
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6(masterThread._map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6) 

{
}


void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithWorkerThread(const MeshRefinementAndPlotGrid& workerThread) {
  _map2PreProcessing.mergeWithWorkerThread(workerThread._map2PreProcessing);
  _map2MeshRefinement.mergeWithWorkerThread(workerThread._map2MeshRefinement);
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithWorkerThread(workerThread._map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2);
  _map2AugmentedAMRTreePlot2d.mergeWithWorkerThread(workerThread._map2AugmentedAMRTreePlot2d);
  _map2LoadBalancing.mergeWithWorkerThread(workerThread._map2LoadBalancing);
  _map2PostProcessing.mergeWithWorkerThread(workerThread._map2PostProcessing);
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithWorkerThread(workerThread._map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6);

}
#endif


void exahype::adapters::MeshRefinementAndPlotGrid::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2PreProcessing.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void exahype::adapters::MeshRefinementAndPlotGrid::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2PreProcessing.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::MeshRefinementAndPlotGrid::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2PreProcessing.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::MeshRefinementAndPlotGrid::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2PreProcessing.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::MeshRefinementAndPlotGrid::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2PreProcessing.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::MeshRefinementAndPlotGrid::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2PreProcessing.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinement.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2AugmentedAMRTreePlot2d.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2LoadBalancing.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::MeshRefinementAndPlotGrid::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2PreProcessing.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinement.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2AugmentedAMRTreePlot2d.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2LoadBalancing.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


#ifdef Parallel
void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2PreProcessing.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2MeshRefinement.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2AugmentedAMRTreePlot2d.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2LoadBalancing.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2PostProcessing.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}


void exahype::adapters::MeshRefinementAndPlotGrid::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2PreProcessing.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2MeshRefinement.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2AugmentedAMRTreePlot2d.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2LoadBalancing.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2PostProcessing.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.prepareSendToNeighbour( vertex, toRank, x, h, level );

}


void exahype::adapters::MeshRefinementAndPlotGrid::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2PreProcessing.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2MeshRefinement.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2AugmentedAMRTreePlot2d.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2LoadBalancing.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2PostProcessing.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}


void exahype::adapters::MeshRefinementAndPlotGrid::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2PreProcessing.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2MeshRefinement.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2AugmentedAMRTreePlot2d.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2LoadBalancing.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2PostProcessing.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}


void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2PreProcessing.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2MeshRefinement.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2AugmentedAMRTreePlot2d.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2LoadBalancing.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2PostProcessing.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}


void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2PreProcessing.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2MeshRefinement.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2AugmentedAMRTreePlot2d.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2LoadBalancing.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2PostProcessing.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}


bool exahype::adapters::MeshRefinementAndPlotGrid::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  bool result = false;
   result |= _map2PreProcessing.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2MeshRefinement.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2AugmentedAMRTreePlot2d.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2LoadBalancing.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2PostProcessing.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}


void exahype::adapters::MeshRefinementAndPlotGrid::prepareSendToMaster(
  exahype::Cell&                       localCell,
  exahype::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
   _map2PreProcessing.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2MeshRefinement.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2AugmentedAMRTreePlot2d.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2LoadBalancing.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2PostProcessing.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithMaster(
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
   _map2PreProcessing.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2MeshRefinement.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2AugmentedAMRTreePlot2d.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2LoadBalancing.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2PostProcessing.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}


void exahype::adapters::MeshRefinementAndPlotGrid::receiveDataFromMaster(
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
   _map2PreProcessing.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2MeshRefinement.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2AugmentedAMRTreePlot2d.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2LoadBalancing.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2PostProcessing.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2PreProcessing.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2MeshRefinement.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2AugmentedAMRTreePlot2d.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2LoadBalancing.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2PostProcessing.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}


void exahype::adapters::MeshRefinementAndPlotGrid::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2PreProcessing.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2MeshRefinement.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2AugmentedAMRTreePlot2d.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2LoadBalancing.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2PostProcessing.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif


void exahype::adapters::MeshRefinementAndPlotGrid::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2PreProcessing.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::MeshRefinementAndPlotGrid::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2PreProcessing.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinement.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2AugmentedAMRTreePlot2d.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2LoadBalancing.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::MeshRefinementAndPlotGrid::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2PreProcessing.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinement.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2AugmentedAMRTreePlot2d.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2LoadBalancing.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::MeshRefinementAndPlotGrid::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  _map2PreProcessing.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinement.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2AugmentedAMRTreePlot2d.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2LoadBalancing.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::MeshRefinementAndPlotGrid::beginIteration(
  exahype::State&  solverState
) {
  _map2PreProcessing.beginIteration( solverState );
  _map2MeshRefinement.beginIteration( solverState );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.beginIteration( solverState );
  _map2AugmentedAMRTreePlot2d.beginIteration( solverState );
  _map2LoadBalancing.beginIteration( solverState );
  _map2PostProcessing.beginIteration( solverState );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.beginIteration( solverState );

}


void exahype::adapters::MeshRefinementAndPlotGrid::endIteration(
  exahype::State&  solverState
) {
  _map2PreProcessing.endIteration( solverState );
  _map2MeshRefinement.endIteration( solverState );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.endIteration( solverState );
  _map2AugmentedAMRTreePlot2d.endIteration( solverState );
  _map2LoadBalancing.endIteration( solverState );
  _map2PostProcessing.endIteration( solverState );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.endIteration( solverState );

}




void exahype::adapters::MeshRefinementAndPlotGrid::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
  _map2PreProcessing.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2MeshRefinement.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2AugmentedAMRTreePlot2d.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2LoadBalancing.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PostProcessing.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void exahype::adapters::MeshRefinementAndPlotGrid::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
  _map2PreProcessing.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2MeshRefinement.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2MeshRefinementAndPlotGrid2VTKGridVisualiser_2.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2AugmentedAMRTreePlot2d.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2LoadBalancing.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PostProcessing.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2MeshRefinementAndPlotGrid2MultiscaleLinkedCell_6.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
