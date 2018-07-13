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

#ifndef ADAPTIVEMESHREFINEMENT_H_
#define ADAPTIVEMESHREFINEMENT_H_

#include "peano/utils/Globals.h"

#include "exahype/Cell.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

namespace exahype {

namespace amr {
  /**
   * Per coordinate direction xi, count the number of shifts
   * of step size \p childSize(xi) necessary to
   * reach \p childOffset from \p parentOffset.
   *
   * \param[in] childOffset  Offset of a child cell.
   * \param[in] childSize    Size of the child cell.
   * \param[in] parentOffset Offset of the parent cell.
   *
   * \see getSubfaceIndex
   */
  tarch::la::Vector<DIMENSIONS,int> computeSubcellIndex(
        const tarch::la::Vector<DIMENSIONS,double>& childOffset,
        const tarch::la::Vector<DIMENSIONS,double>& childSize,
        const tarch::la::Vector<DIMENSIONS,double>& parentOffset);

  /**
   * Collect all the element with index!=d
   * from \p subcellIndex.
   *
   * \see getSubcellIndex
   */
  tarch::la::Vector<DIMENSIONS-1,int> getSubfaceIndex(
        const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
        const int d);

  /**
   * Per coordinate direction xi, check if
   * subcellIndex[xi] is either 0 or 3^levelDelta - 1.
   * If this is the case for at least one subcellIndex[xi]
   * return true. Otherwise return false.
   */
  bool onBoundaryOfParent(
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
      const int levelDelta);

  /**
   * This method deemed to be necessary since I could
   * not trust the Boundary markers.
   *
   * Per coordinate direction xi, check if the cell
   * is at the boundary of the parent cell, i.e. if
   * subcellIndex[xi] is either 0 or 3^levelDelta - 1,
   * and if the parent face is a boundary face.
   * If so, the face of the child is on the boundary as well.
   */
  std::bitset<DIMENSIONS_TIMES_TWO> determineInsideAndOutsideFacesOfChild(
      const tarch::la::Vector<DIMENSIONS,double>& childOffset,
      const tarch::la::Vector<DIMENSIONS,double>& childSize,
      const tarch::la::Vector<DIMENSIONS,double>& parentOffset,
      const int                                   levelDelta,
      const std::bitset<DIMENSIONS_TIMES_TWO>&    parentIsInside);

  /**
   * Determine the position of a Descendant with respect
   * to a  Cell or Descendant that contains data, i.e.,
   * has at least one neighbour that is a real cell.
   *
   * \note This function only makes sense if the
   * Descendant has a valid parentIndex attribute.
   *
   * This method is required for the face data prolongation, the
   * volume data prolongation (!), and the FV volume data prolongation (!).
   *
   * \param topMost Set to true if you want to lookup the
   * rank-local top-most parent of the descendant which
   * might either be of type Cell or a Descendant at the
   * master-worker boundaryl
   */
  template <class CellDescription,class CellDescriptionHeap, bool topMost>
  exahype::solvers::Solver::SubcellPosition
  computeSubcellPositionOfDescendant(const CellDescription& pChild);
}  // namespace amr
}  // namespace exahype


#include "AdaptiveMeshRefinement.cpph"


#endif /* ADAPTIVEMESHREFINEMENT_H_ */
