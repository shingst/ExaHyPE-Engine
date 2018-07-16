#include "exahype/records/ADERDGCellDescription.h"

#if defined(Parallel)
   exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _solverNumber(solverNumber),
   _neighbourMergePerformed(neighbourMergePerformed),
   _adjacentToRemoteRank(adjacentToRemoteRank),
   _hasToHoldDataForMasterWorkerCommunication(hasToHoldDataForMasterWorkerCommunication),
   _faceDataExchangeCounter(faceDataExchangeCounter),
   _parentIndex(parentIndex),
   _parentCellLevel(parentCellLevel),
   _parentOffset(parentOffset),
   _hasVirtualChildren(hasVirtualChildren),
   _type(type),
   _refinementEvent(refinementEvent),
   _level(level),
   _offset(offset),
   _size(size),
   _previousCorrectorTimeStamp(previousCorrectorTimeStamp),
   _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
   _correctorTimeStepSize(correctorTimeStepSize),
   _correctorTimeStamp(correctorTimeStamp),
   _predictorTimeStepSize(predictorTimeStepSize),
   _predictorTimeStamp(predictorTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _update(update),
   _updateAverages(updateAverages),
   _updateCompressed(updateCompressed),
   _extrapolatedPredictor(extrapolatedPredictor),
   _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
   _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
   _fluctuation(fluctuation),
   _fluctuationAverages(fluctuationAverages),
   _fluctuationCompressed(fluctuationCompressed),
   _solutionMin(solutionMin),
   _solutionMax(solutionMax),
   _facewiseAugmentationStatus(facewiseAugmentationStatus),
   _augmentationStatus(augmentationStatus),
   _previousAugmentationStatus(previousAugmentationStatus),
   _facewiseCommunicationStatus(facewiseCommunicationStatus),
   _communicationStatus(communicationStatus),
   _facewiseRefinementStatus(facewiseRefinementStatus),
   _refinementFlag(refinementFlag),
   _refinementStatus(refinementStatus),
   _previousRefinementStatus(previousRefinementStatus),
   _iterationsToCureTroubledCell(iterationsToCureTroubledCell),
   _compressionState(compressionState),
   _bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
   _bytesPerDoFInSolution(bytesPerDoFInSolution),
   _bytesPerDoFInUpdate(bytesPerDoFInUpdate),
   _bytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor),
   _bytesPerDoFInFluctuation(bytesPerDoFInFluctuation) {
      
   }
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription() {
      
   }
   
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords._adjacentToRemoteRank, persistentRecords._hasToHoldDataForMasterWorkerCommunication, persistentRecords._faceDataExchangeCounter, persistentRecords._parentIndex, persistentRecords._parentCellLevel, persistentRecords._parentOffset, persistentRecords._hasVirtualChildren, persistentRecords._type, persistentRecords._refinementEvent, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._previousCorrectorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._facewiseAugmentationStatus, persistentRecords._augmentationStatus, persistentRecords._previousAugmentationStatus, persistentRecords._facewiseCommunicationStatus, persistentRecords._communicationStatus, persistentRecords._facewiseRefinementStatus, persistentRecords._refinementFlag, persistentRecords._refinementStatus, persistentRecords._previousRefinementStatus, persistentRecords._iterationsToCureTroubledCell, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInUpdate, persistentRecords._bytesPerDoFInExtrapolatedPredictor, persistentRecords._bytesPerDoFInFluctuation) {
      
   }
   
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _persistentRecords(solverNumber, neighbourMergePerformed, adjacentToRemoteRank, hasToHoldDataForMasterWorkerCommunication, faceDataExchangeCounter, parentIndex, parentCellLevel, parentOffset, hasVirtualChildren, type, refinementEvent, level, offset, size, previousCorrectorTimeStamp, previousCorrectorTimeStepSize, correctorTimeStepSize, correctorTimeStamp, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, facewiseAugmentationStatus, augmentationStatus, previousAugmentationStatus, facewiseCommunicationStatus, communicationStatus, facewiseRefinementStatus, refinementFlag, refinementStatus, previousRefinementStatus, iterationsToCureTroubledCell, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
      
   }
   
   
   exahype::records::ADERDGCellDescription::~ADERDGCellDescription() { }
   
   std::string exahype::records::ADERDGCellDescription::toString(const CompressionState& param) {
      switch (param) {
         case Uncompressed: return "Uncompressed";
         case CurrentlyProcessed: return "CurrentlyProcessed";
         case Compressed: return "Compressed";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getCompressionStateMapping() {
      return "CompressionState(Uncompressed=0,CurrentlyProcessed=1,Compressed=2)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const RefinementEvent& param) {
      switch (param) {
         case None: return "None";
         case ErasingChildrenRequested: return "ErasingChildrenRequested";
         case ErasingChildren: return "ErasingChildren";
         case ChangeChildrenToVirtualChildrenRequested: return "ChangeChildrenToVirtualChildrenRequested";
         case ChangeChildrenToVirtualChildren: return "ChangeChildrenToVirtualChildren";
         case RefiningRequested: return "RefiningRequested";
         case Refining: return "Refining";
         case Prolongating: return "Prolongating";
         case ErasingVirtualChildrenRequested: return "ErasingVirtualChildrenRequested";
         case ErasingVirtualChildren: return "ErasingVirtualChildren";
         case VirtualRefiningRequested: return "VirtualRefiningRequested";
         case VirtualRefining: return "VirtualRefining";
         case Erasing: return "Erasing";
         case ChangeToVirtualCell: return "ChangeToVirtualCell";
         case ErasingVirtualCell: return "ErasingVirtualCell";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getRefinementEventMapping() {
      return "RefinementEvent(None=0,ErasingChildrenRequested=1,ErasingChildren=2,ChangeChildrenToVirtualChildrenRequested=3,ChangeChildrenToVirtualChildren=4,RefiningRequested=5,Refining=6,Prolongating=7,ErasingVirtualChildrenRequested=8,ErasingVirtualChildren=9,VirtualRefiningRequested=10,VirtualRefining=11,Erasing=12,ChangeToVirtualCell=13,ErasingVirtualCell=14)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const Type& param) {
      switch (param) {
         case Erased: return "Erased";
         case Ancestor: return "Ancestor";
         case Cell: return "Cell";
         case Descendant: return "Descendant";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getTypeMapping() {
      return "Type(Erased=0,Ancestor=1,Cell=2,Descendant=3)";
   }
   
   
   std::string exahype::records::ADERDGCellDescription::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::ADERDGCellDescription::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "adjacentToRemoteRank:" << getAdjacentToRemoteRank();
      out << ",";
      out << "hasToHoldDataForMasterWorkerCommunication:" << getHasToHoldDataForMasterWorkerCommunication();
      out << ",";
      out << "faceDataExchangeCounter:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFaceDataExchangeCounter(i) << ",";
   }
   out << getFaceDataExchangeCounter(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "parentCellLevel:" << getParentCellLevel();
      out << ",";
      out << "parentOffset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getParentOffset(i) << ",";
   }
   out << getParentOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "hasVirtualChildren:" << getHasVirtualChildren();
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out << ",";
      out << "level:" << getLevel();
      out << ",";
      out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
      out << ",";
      out << "previousCorrectorTimeStamp:" << getPreviousCorrectorTimeStamp();
      out << ",";
      out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStamp:" << getCorrectorTimeStamp();
      out << ",";
      out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
      out << ",";
      out << "predictorTimeStamp:" << getPredictorTimeStamp();
      out << ",";
      out << "solution:" << getSolution();
      out << ",";
      out << "solutionAverages:" << getSolutionAverages();
      out << ",";
      out << "solutionCompressed:" << getSolutionCompressed();
      out << ",";
      out << "previousSolution:" << getPreviousSolution();
      out << ",";
      out << "previousSolutionAverages:" << getPreviousSolutionAverages();
      out << ",";
      out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
      out << ",";
      out << "update:" << getUpdate();
      out << ",";
      out << "updateAverages:" << getUpdateAverages();
      out << ",";
      out << "updateCompressed:" << getUpdateCompressed();
      out << ",";
      out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
      out << ",";
      out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
      out << ",";
      out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
      out << ",";
      out << "fluctuation:" << getFluctuation();
      out << ",";
      out << "fluctuationAverages:" << getFluctuationAverages();
      out << ",";
      out << "fluctuationCompressed:" << getFluctuationCompressed();
      out << ",";
      out << "solutionMin:" << getSolutionMin();
      out << ",";
      out << "solutionMax:" << getSolutionMax();
      out << ",";
      out << "facewiseAugmentationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseAugmentationStatus(i) << ",";
   }
   out << getFacewiseAugmentationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "augmentationStatus:" << getAugmentationStatus();
      out << ",";
      out << "previousAugmentationStatus:" << getPreviousAugmentationStatus();
      out << ",";
      out << "facewiseCommunicationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseCommunicationStatus(i) << ",";
   }
   out << getFacewiseCommunicationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "communicationStatus:" << getCommunicationStatus();
      out << ",";
      out << "facewiseRefinementStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseRefinementStatus(i) << ",";
   }
   out << getFacewiseRefinementStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "refinementFlag:" << getRefinementFlag();
      out << ",";
      out << "refinementStatus:" << getRefinementStatus();
      out << ",";
      out << "previousRefinementStatus:" << getPreviousRefinementStatus();
      out << ",";
      out << "iterationsToCureTroubledCell:" << getIterationsToCureTroubledCell();
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
      out << ",";
      out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
      out << ",";
      out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
      out <<  ")";
   }
   
   
   exahype::records::ADERDGCellDescription::PersistentRecords exahype::records::ADERDGCellDescription::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::ADERDGCellDescriptionPacked exahype::records::ADERDGCellDescription::convert() const{
      return ADERDGCellDescriptionPacked(
         getSolverNumber(),
         getNeighbourMergePerformed(),
         getAdjacentToRemoteRank(),
         getHasToHoldDataForMasterWorkerCommunication(),
         getFaceDataExchangeCounter(),
         getParentIndex(),
         getParentCellLevel(),
         getParentOffset(),
         getHasVirtualChildren(),
         getType(),
         getRefinementEvent(),
         getLevel(),
         getOffset(),
         getSize(),
         getPreviousCorrectorTimeStamp(),
         getPreviousCorrectorTimeStepSize(),
         getCorrectorTimeStepSize(),
         getCorrectorTimeStamp(),
         getPredictorTimeStepSize(),
         getPredictorTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getUpdate(),
         getUpdateAverages(),
         getUpdateCompressed(),
         getExtrapolatedPredictor(),
         getExtrapolatedPredictorAverages(),
         getExtrapolatedPredictorCompressed(),
         getFluctuation(),
         getFluctuationAverages(),
         getFluctuationCompressed(),
         getSolutionMin(),
         getSolutionMax(),
         getFacewiseAugmentationStatus(),
         getAugmentationStatus(),
         getPreviousAugmentationStatus(),
         getFacewiseCommunicationStatus(),
         getCommunicationStatus(),
         getFacewiseRefinementStatus(),
         getRefinementFlag(),
         getRefinementStatus(),
         getPreviousRefinementStatus(),
         getIterationsToCureTroubledCell(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInUpdate(),
         getBytesPerDoFInExtrapolatedPredictor(),
         getBytesPerDoFInFluctuation()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::ADERDGCellDescription::_log( "exahype::records::ADERDGCellDescription" );
      
      MPI_Datatype exahype::records::ADERDGCellDescription::Datatype = 0;
      MPI_Datatype exahype::records::ADERDGCellDescription::FullDatatype = 0;
      
      
      void exahype::records::ADERDGCellDescription::initDatatype() {
         {
            ADERDGCellDescription dummyADERDGCellDescription[2];
            
            #ifdef MPI2
            const int Attributes = 28;
            #else
            const int Attributes = 29;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //hasToHoldDataForMasterWorkerCommunication
               , MPI_INT		 //parentCellLevel
               , MPI_DOUBLE		 //parentOffset
               , MPI_CXX_BOOL		 //hasVirtualChildren
               , MPI_INT		 //type
               , MPI_INT		 //refinementEvent
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_DOUBLE		 //previousCorrectorTimeStamp
               , MPI_DOUBLE		 //previousCorrectorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStamp
               , MPI_DOUBLE		 //predictorTimeStepSize
               , MPI_DOUBLE		 //predictorTimeStamp
               , MPI_INT		 //augmentationStatus
               , MPI_INT		 //previousAugmentationStatus
               , MPI_INT		 //communicationStatus
               , MPI_INT		 //refinementStatus
               , MPI_INT		 //previousRefinementStatus
               , MPI_INT		 //iterationsToCureTroubledCell
               , MPI_INT		 //compressionState
               , MPI_INT		 //bytesPerDoFInPreviousSolution
               , MPI_INT		 //bytesPerDoFInSolution
               , MPI_INT		 //bytesPerDoFInUpdate
               , MPI_INT		 //bytesPerDoFInExtrapolatedPredictor
               , MPI_INT		 //bytesPerDoFInFluctuation
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , 1		 //hasToHoldDataForMasterWorkerCommunication
               , 1		 //parentCellLevel
               , DIMENSIONS		 //parentOffset
               , 1		 //hasVirtualChildren
               , 1		 //type
               , 1		 //refinementEvent
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //previousCorrectorTimeStamp
               , 1		 //previousCorrectorTimeStepSize
               , 1		 //correctorTimeStepSize
               , 1		 //correctorTimeStamp
               , 1		 //predictorTimeStepSize
               , 1		 //predictorTimeStamp
               , 1		 //augmentationStatus
               , 1		 //previousAugmentationStatus
               , 1		 //communicationStatus
               , 1		 //refinementStatus
               , 1		 //previousRefinementStatus
               , 1		 //iterationsToCureTroubledCell
               , 1		 //compressionState
               , 1		 //bytesPerDoFInPreviousSolution
               , 1		 //bytesPerDoFInSolution
               , 1		 //bytesPerDoFInUpdate
               , 1		 //bytesPerDoFInExtrapolatedPredictor
               , 1		 //bytesPerDoFInFluctuation
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[23] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[24] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[24] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[25] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[25] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[26] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[26] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[27] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[27] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescription)), i, disp[i], Attributes, sizeof(ADERDGCellDescription));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]))), 		&disp[28] );
            disp[28] -= base;
            disp[28] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescription::Datatype );
            MPI_Type_commit( &ADERDGCellDescription::Datatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::Datatype);
            MPI_Type_commit( &ADERDGCellDescription::Datatype );
            #endif
            
         }
         {
            ADERDGCellDescription dummyADERDGCellDescription[2];
            
            #ifdef MPI2
            const int Attributes = 53;
            #else
            const int Attributes = 54;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //neighbourMergePerformed
               , MPI_CXX_BOOL		 //adjacentToRemoteRank
               , MPI_CXX_BOOL		 //hasToHoldDataForMasterWorkerCommunication
               , MPI_INT		 //faceDataExchangeCounter
               , MPI_INT		 //parentIndex
               , MPI_INT		 //parentCellLevel
               , MPI_DOUBLE		 //parentOffset
               , MPI_CXX_BOOL		 //hasVirtualChildren
               , MPI_INT		 //type
               , MPI_INT		 //refinementEvent
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_DOUBLE		 //previousCorrectorTimeStamp
               , MPI_DOUBLE		 //previousCorrectorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStamp
               , MPI_DOUBLE		 //predictorTimeStepSize
               , MPI_DOUBLE		 //predictorTimeStamp
               , MPI_INT		 //solution
               , MPI_INT		 //solutionAverages
               , MPI_INT		 //solutionCompressed
               , MPI_INT		 //previousSolution
               , MPI_INT		 //previousSolutionAverages
               , MPI_INT		 //previousSolutionCompressed
               , MPI_INT		 //update
               , MPI_INT		 //updateAverages
               , MPI_INT		 //updateCompressed
               , MPI_INT		 //extrapolatedPredictor
               , MPI_INT		 //extrapolatedPredictorAverages
               , MPI_INT		 //extrapolatedPredictorCompressed
               , MPI_INT		 //fluctuation
               , MPI_INT		 //fluctuationAverages
               , MPI_INT		 //fluctuationCompressed
               , MPI_INT		 //solutionMin
               , MPI_INT		 //solutionMax
               , MPI_INT		 //facewiseAugmentationStatus
               , MPI_INT		 //augmentationStatus
               , MPI_INT		 //previousAugmentationStatus
               , MPI_INT		 //facewiseCommunicationStatus
               , MPI_INT		 //communicationStatus
               , MPI_INT		 //facewiseRefinementStatus
               , MPI_CXX_BOOL		 //refinementFlag
               , MPI_INT		 //refinementStatus
               , MPI_INT		 //previousRefinementStatus
               , MPI_INT		 //iterationsToCureTroubledCell
               , MPI_INT		 //compressionState
               , MPI_INT		 //bytesPerDoFInPreviousSolution
               , MPI_INT		 //bytesPerDoFInSolution
               , MPI_INT		 //bytesPerDoFInUpdate
               , MPI_INT		 //bytesPerDoFInExtrapolatedPredictor
               , MPI_INT		 //bytesPerDoFInFluctuation
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
               , 1		 //adjacentToRemoteRank
               , 1		 //hasToHoldDataForMasterWorkerCommunication
               , DIMENSIONS_TIMES_TWO		 //faceDataExchangeCounter
               , 1		 //parentIndex
               , 1		 //parentCellLevel
               , DIMENSIONS		 //parentOffset
               , 1		 //hasVirtualChildren
               , 1		 //type
               , 1		 //refinementEvent
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //previousCorrectorTimeStamp
               , 1		 //previousCorrectorTimeStepSize
               , 1		 //correctorTimeStepSize
               , 1		 //correctorTimeStamp
               , 1		 //predictorTimeStepSize
               , 1		 //predictorTimeStamp
               , 1		 //solution
               , 1		 //solutionAverages
               , 1		 //solutionCompressed
               , 1		 //previousSolution
               , 1		 //previousSolutionAverages
               , 1		 //previousSolutionCompressed
               , 1		 //update
               , 1		 //updateAverages
               , 1		 //updateCompressed
               , 1		 //extrapolatedPredictor
               , 1		 //extrapolatedPredictorAverages
               , 1		 //extrapolatedPredictorCompressed
               , 1		 //fluctuation
               , 1		 //fluctuationAverages
               , 1		 //fluctuationCompressed
               , 1		 //solutionMin
               , 1		 //solutionMax
               , DIMENSIONS_TIMES_TWO		 //facewiseAugmentationStatus
               , 1		 //augmentationStatus
               , 1		 //previousAugmentationStatus
               , DIMENSIONS_TIMES_TWO		 //facewiseCommunicationStatus
               , 1		 //communicationStatus
               , DIMENSIONS_TIMES_TWO		 //facewiseRefinementStatus
               , 1		 //refinementFlag
               , 1		 //refinementStatus
               , 1		 //previousRefinementStatus
               , 1		 //iterationsToCureTroubledCell
               , 1		 //compressionState
               , 1		 //bytesPerDoFInPreviousSolution
               , 1		 //bytesPerDoFInSolution
               , 1		 //bytesPerDoFInUpdate
               , 1		 //bytesPerDoFInExtrapolatedPredictor
               , 1		 //bytesPerDoFInFluctuation
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._adjacentToRemoteRank))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._adjacentToRemoteRank))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[23] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[24] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[24] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[25] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[25] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[26] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[26] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[27] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[27] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[28] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[28] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[29] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[29] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[30] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[30] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[31] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[31] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[32] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[32] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[33] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[33] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[34] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[34] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[35] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[35] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[36] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[36] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[37] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[37] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[38] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[38] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[39] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[39] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[40] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[40] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[41] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[41] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[42] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[42] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[43] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[43] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[44] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[44] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[45] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[45] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[46] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[46] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[47] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[47] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[48] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[48] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[49] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[49] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[50] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[50] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[51] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[51] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[52] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[52] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescription)), i, disp[i], Attributes, sizeof(ADERDGCellDescription));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]))), 		&disp[53] );
            disp[53] -= base;
            disp[53] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescription::FullDatatype );
            MPI_Type_commit( &ADERDGCellDescription::FullDatatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::FullDatatype);
            MPI_Type_commit( &ADERDGCellDescription::FullDatatype );
            #endif
            
         }
         
      }
      
      
      void exahype::records::ADERDGCellDescription::shutdownDatatype() {
         MPI_Type_free( &ADERDGCellDescription::Datatype );
         MPI_Type_free( &ADERDGCellDescription::FullDatatype );
         
      }
      
      void exahype::records::ADERDGCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
         // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::ADERDGCellDescription " 
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result); 
         _log.error( "send(int)",msg.str() ); 
       } 
    } 
    break; 
   case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    {
      MPI_Request* sendRequestHandle = new MPI_Request(); 
      int          flag = 0; 
       int          result; 
       clock_t      timeOutWarning   = -1; 
       clock_t      timeOutShutdown  = -1; 
       bool         triggeredTimeoutWarning = false;  
       result = MPI_Isend(  
         this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination,  
         tag, tarch::parallel::Node::getInstance().getCommunicator(), 
         sendRequestHandle  
       ); 
       if  (result!=MPI_SUCCESS) {  
         std::ostringstream msg;  
         msg << "was not able to send message exahype::records::ADERDGCellDescription "  
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result);  
         _log.error( "send(int)",msg.str() );  
       }  
       result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
       while (!flag) { 
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
         result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
           std::ostringstream msg; 
           msg << "testing for finished send task for exahype::records::ADERDGCellDescription " 
               << toString() 
               << " sent to node " << destination 
               << " failed: " << tarch::parallel::MPIReturnValueToString(result); 
           _log.error("send(int)", msg.str() ); 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
           (clock()>timeOutWarning) && 
           (!triggeredTimeoutWarning) 
         ) { 
           tarch::parallel::Node::getInstance().writeTimeOutWarning( 
             "exahype::records::ADERDGCellDescription", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::ADERDGCellDescription", 
             "send(int)", destination,tag,1 
           ); 
         } 
 	       tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
       } 
       delete sendRequestHandle; 
     }  
     break; 
   case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    assertionMsg(false,"should not be called"); 
    break; 
} 
 // ============================= 
// end injected snippet/aspect 
// ============================= 

         
      }
      
      
      
      void exahype::records::ADERDGCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
         // ============================= 
// start injected snippet/aspect 
// ============================= 
MPI_Status status; 
switch (mode) { 
  case ExchangeMode::Blocking: 
    { 
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescription from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    } 
    break; 
  case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    { 
      int          flag = 0; 
      int          result; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      MPI_Request* sendRequestHandle = new MPI_Request(); 
       result = MPI_Irecv( 
        this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, 
        tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle 
      ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescription from node " 
             << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
      result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      delete sendRequestHandle; 
    }    break; 
  case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    {
      int flag; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      int result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
       if (result!=MPI_SUCCESS) { 
        std::ostringstream msg; 
        msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: " 
            << tarch::parallel::MPIReturnValueToString(result); 
        _log.error("receive(int)", msg.str() ); 
      } 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescription from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    }
    break; 
  } 
// =========================== 
// end injected snippet/aspect 
// =========================== 

         
      }
      
      
      
      bool exahype::records::ADERDGCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
         MPI_Status status;
         int  flag        = 0;
         MPI_Iprobe(
            MPI_ANY_SOURCE, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
         );
         if (flag) {
            int  messageCounter;
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               MPI_Get_count(&status, Datatype, &messageCounter);
            }
            else {
               MPI_Get_count(&status, FullDatatype, &messageCounter);
            }
            return messageCounter > 0;
         }
         else return false;
         
      }
      
      
   #endif
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords() {
      if ((DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _solverNumber(solverNumber),
   _hasToHoldDataForMasterWorkerCommunication(hasToHoldDataForMasterWorkerCommunication),
   _faceDataExchangeCounter(faceDataExchangeCounter),
   _parentIndex(parentIndex),
   _parentCellLevel(parentCellLevel),
   _parentOffset(parentOffset),
   _hasVirtualChildren(hasVirtualChildren),
   _level(level),
   _offset(offset),
   _size(size),
   _previousCorrectorTimeStamp(previousCorrectorTimeStamp),
   _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
   _correctorTimeStepSize(correctorTimeStepSize),
   _correctorTimeStamp(correctorTimeStamp),
   _predictorTimeStepSize(predictorTimeStepSize),
   _predictorTimeStamp(predictorTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _update(update),
   _updateAverages(updateAverages),
   _updateCompressed(updateCompressed),
   _extrapolatedPredictor(extrapolatedPredictor),
   _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
   _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
   _fluctuation(fluctuation),
   _fluctuationAverages(fluctuationAverages),
   _fluctuationCompressed(fluctuationCompressed),
   _solutionMin(solutionMin),
   _solutionMax(solutionMax),
   _facewiseAugmentationStatus(facewiseAugmentationStatus),
   _augmentationStatus(augmentationStatus),
   _previousAugmentationStatus(previousAugmentationStatus),
   _facewiseCommunicationStatus(facewiseCommunicationStatus),
   _communicationStatus(communicationStatus),
   _facewiseRefinementStatus(facewiseRefinementStatus),
   _refinementFlag(refinementFlag),
   _refinementStatus(refinementStatus),
   _previousRefinementStatus(previousRefinementStatus),
   _iterationsToCureTroubledCell(iterationsToCureTroubledCell) {
      setNeighbourMergePerformed(neighbourMergePerformed);
      setAdjacentToRemoteRank(adjacentToRemoteRank);
      setType(type);
      setRefinementEvent(refinementEvent);
      setCompressionState(compressionState);
      setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
      setBytesPerDoFInSolution(bytesPerDoFInSolution);
      setBytesPerDoFInUpdate(bytesPerDoFInUpdate);
      setBytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor);
      setBytesPerDoFInFluctuation(bytesPerDoFInFluctuation);
      if ((DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked() {
      if ((DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords.getNeighbourMergePerformed(), persistentRecords.getAdjacentToRemoteRank(), persistentRecords._hasToHoldDataForMasterWorkerCommunication, persistentRecords._faceDataExchangeCounter, persistentRecords._parentIndex, persistentRecords._parentCellLevel, persistentRecords._parentOffset, persistentRecords._hasVirtualChildren, persistentRecords.getType(), persistentRecords.getRefinementEvent(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._previousCorrectorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._facewiseAugmentationStatus, persistentRecords._augmentationStatus, persistentRecords._previousAugmentationStatus, persistentRecords._facewiseCommunicationStatus, persistentRecords._communicationStatus, persistentRecords._facewiseRefinementStatus, persistentRecords._refinementFlag, persistentRecords._refinementStatus, persistentRecords._previousRefinementStatus, persistentRecords._iterationsToCureTroubledCell, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInUpdate(), persistentRecords.getBytesPerDoFInExtrapolatedPredictor(), persistentRecords.getBytesPerDoFInFluctuation()) {
      if ((DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _persistentRecords(solverNumber, neighbourMergePerformed, adjacentToRemoteRank, hasToHoldDataForMasterWorkerCommunication, faceDataExchangeCounter, parentIndex, parentCellLevel, parentOffset, hasVirtualChildren, type, refinementEvent, level, offset, size, previousCorrectorTimeStamp, previousCorrectorTimeStepSize, correctorTimeStepSize, correctorTimeStamp, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, facewiseAugmentationStatus, augmentationStatus, previousAugmentationStatus, facewiseCommunicationStatus, communicationStatus, facewiseRefinementStatus, refinementFlag, refinementStatus, previousRefinementStatus, iterationsToCureTroubledCell, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
      if ((DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::~ADERDGCellDescriptionPacked() { }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const Type& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getTypeMapping() {
      return exahype::records::ADERDGCellDescription::getTypeMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const RefinementEvent& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getRefinementEventMapping() {
      return exahype::records::ADERDGCellDescription::getRefinementEventMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const CompressionState& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getCompressionStateMapping() {
      return exahype::records::ADERDGCellDescription::getCompressionStateMapping();
   }
   
   
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::ADERDGCellDescriptionPacked::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "adjacentToRemoteRank:" << getAdjacentToRemoteRank();
      out << ",";
      out << "hasToHoldDataForMasterWorkerCommunication:" << getHasToHoldDataForMasterWorkerCommunication();
      out << ",";
      out << "faceDataExchangeCounter:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFaceDataExchangeCounter(i) << ",";
   }
   out << getFaceDataExchangeCounter(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "parentCellLevel:" << getParentCellLevel();
      out << ",";
      out << "parentOffset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getParentOffset(i) << ",";
   }
   out << getParentOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "hasVirtualChildren:" << getHasVirtualChildren();
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out << ",";
      out << "level:" << getLevel();
      out << ",";
      out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
      out << ",";
      out << "previousCorrectorTimeStamp:" << getPreviousCorrectorTimeStamp();
      out << ",";
      out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStamp:" << getCorrectorTimeStamp();
      out << ",";
      out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
      out << ",";
      out << "predictorTimeStamp:" << getPredictorTimeStamp();
      out << ",";
      out << "solution:" << getSolution();
      out << ",";
      out << "solutionAverages:" << getSolutionAverages();
      out << ",";
      out << "solutionCompressed:" << getSolutionCompressed();
      out << ",";
      out << "previousSolution:" << getPreviousSolution();
      out << ",";
      out << "previousSolutionAverages:" << getPreviousSolutionAverages();
      out << ",";
      out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
      out << ",";
      out << "update:" << getUpdate();
      out << ",";
      out << "updateAverages:" << getUpdateAverages();
      out << ",";
      out << "updateCompressed:" << getUpdateCompressed();
      out << ",";
      out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
      out << ",";
      out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
      out << ",";
      out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
      out << ",";
      out << "fluctuation:" << getFluctuation();
      out << ",";
      out << "fluctuationAverages:" << getFluctuationAverages();
      out << ",";
      out << "fluctuationCompressed:" << getFluctuationCompressed();
      out << ",";
      out << "solutionMin:" << getSolutionMin();
      out << ",";
      out << "solutionMax:" << getSolutionMax();
      out << ",";
      out << "facewiseAugmentationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseAugmentationStatus(i) << ",";
   }
   out << getFacewiseAugmentationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "augmentationStatus:" << getAugmentationStatus();
      out << ",";
      out << "previousAugmentationStatus:" << getPreviousAugmentationStatus();
      out << ",";
      out << "facewiseCommunicationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseCommunicationStatus(i) << ",";
   }
   out << getFacewiseCommunicationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "communicationStatus:" << getCommunicationStatus();
      out << ",";
      out << "facewiseRefinementStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseRefinementStatus(i) << ",";
   }
   out << getFacewiseRefinementStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "refinementFlag:" << getRefinementFlag();
      out << ",";
      out << "refinementStatus:" << getRefinementStatus();
      out << ",";
      out << "previousRefinementStatus:" << getPreviousRefinementStatus();
      out << ",";
      out << "iterationsToCureTroubledCell:" << getIterationsToCureTroubledCell();
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
      out << ",";
      out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
      out << ",";
      out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
      out <<  ")";
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords exahype::records::ADERDGCellDescriptionPacked::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::ADERDGCellDescription exahype::records::ADERDGCellDescriptionPacked::convert() const{
      return ADERDGCellDescription(
         getSolverNumber(),
         getNeighbourMergePerformed(),
         getAdjacentToRemoteRank(),
         getHasToHoldDataForMasterWorkerCommunication(),
         getFaceDataExchangeCounter(),
         getParentIndex(),
         getParentCellLevel(),
         getParentOffset(),
         getHasVirtualChildren(),
         getType(),
         getRefinementEvent(),
         getLevel(),
         getOffset(),
         getSize(),
         getPreviousCorrectorTimeStamp(),
         getPreviousCorrectorTimeStepSize(),
         getCorrectorTimeStepSize(),
         getCorrectorTimeStamp(),
         getPredictorTimeStepSize(),
         getPredictorTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getUpdate(),
         getUpdateAverages(),
         getUpdateCompressed(),
         getExtrapolatedPredictor(),
         getExtrapolatedPredictorAverages(),
         getExtrapolatedPredictorCompressed(),
         getFluctuation(),
         getFluctuationAverages(),
         getFluctuationCompressed(),
         getSolutionMin(),
         getSolutionMax(),
         getFacewiseAugmentationStatus(),
         getAugmentationStatus(),
         getPreviousAugmentationStatus(),
         getFacewiseCommunicationStatus(),
         getCommunicationStatus(),
         getFacewiseRefinementStatus(),
         getRefinementFlag(),
         getRefinementStatus(),
         getPreviousRefinementStatus(),
         getIterationsToCureTroubledCell(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInUpdate(),
         getBytesPerDoFInExtrapolatedPredictor(),
         getBytesPerDoFInFluctuation()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::ADERDGCellDescriptionPacked::_log( "exahype::records::ADERDGCellDescriptionPacked" );
      
      MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::Datatype = 0;
      MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::FullDatatype = 0;
      
      
      void exahype::records::ADERDGCellDescriptionPacked::initDatatype() {
         {
            ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
            
            #ifdef MPI2
            const int Attributes = 21;
            #else
            const int Attributes = 22;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //hasToHoldDataForMasterWorkerCommunication
               , MPI_INT		 //parentCellLevel
               , MPI_DOUBLE		 //parentOffset
               , MPI_CXX_BOOL		 //hasVirtualChildren
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_DOUBLE		 //previousCorrectorTimeStamp
               , MPI_DOUBLE		 //previousCorrectorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStamp
               , MPI_DOUBLE		 //predictorTimeStepSize
               , MPI_DOUBLE		 //predictorTimeStamp
               , MPI_INT		 //augmentationStatus
               , MPI_INT		 //previousAugmentationStatus
               , MPI_INT		 //communicationStatus
               , MPI_INT		 //refinementStatus
               , MPI_INT		 //previousRefinementStatus
               , MPI_INT		 //iterationsToCureTroubledCell
               , MPI_INT		 //_packedRecords0
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , 1		 //hasToHoldDataForMasterWorkerCommunication
               , 1		 //parentCellLevel
               , DIMENSIONS		 //parentOffset
               , 1		 //hasVirtualChildren
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //previousCorrectorTimeStamp
               , 1		 //previousCorrectorTimeStepSize
               , 1		 //correctorTimeStepSize
               , 1		 //correctorTimeStamp
               , 1		 //predictorTimeStepSize
               , 1		 //predictorTimeStamp
               , 1		 //augmentationStatus
               , 1		 //previousAugmentationStatus
               , 1		 //communicationStatus
               , 1		 //refinementStatus
               , 1		 //previousRefinementStatus
               , 1		 //iterationsToCureTroubledCell
               , 1		 //_packedRecords0
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[20] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescriptionPacked)), i, disp[i], Attributes, sizeof(ADERDGCellDescriptionPacked));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]))), 		&disp[21] );
            disp[21] -= base;
            disp[21] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescriptionPacked::Datatype );
            MPI_Type_commit( &ADERDGCellDescriptionPacked::Datatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::Datatype);
            MPI_Type_commit( &ADERDGCellDescriptionPacked::Datatype );
            #endif
            
         }
         {
            ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
            
            #ifdef MPI2
            const int Attributes = 44;
            #else
            const int Attributes = 45;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //hasToHoldDataForMasterWorkerCommunication
               , MPI_INT		 //faceDataExchangeCounter
               , MPI_INT		 //parentIndex
               , MPI_INT		 //parentCellLevel
               , MPI_DOUBLE		 //parentOffset
               , MPI_CXX_BOOL		 //hasVirtualChildren
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_DOUBLE		 //previousCorrectorTimeStamp
               , MPI_DOUBLE		 //previousCorrectorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStepSize
               , MPI_DOUBLE		 //correctorTimeStamp
               , MPI_DOUBLE		 //predictorTimeStepSize
               , MPI_DOUBLE		 //predictorTimeStamp
               , MPI_INT		 //solution
               , MPI_INT		 //solutionAverages
               , MPI_INT		 //solutionCompressed
               , MPI_INT		 //previousSolution
               , MPI_INT		 //previousSolutionAverages
               , MPI_INT		 //previousSolutionCompressed
               , MPI_INT		 //update
               , MPI_INT		 //updateAverages
               , MPI_INT		 //updateCompressed
               , MPI_INT		 //extrapolatedPredictor
               , MPI_INT		 //extrapolatedPredictorAverages
               , MPI_INT		 //extrapolatedPredictorCompressed
               , MPI_INT		 //fluctuation
               , MPI_INT		 //fluctuationAverages
               , MPI_INT		 //fluctuationCompressed
               , MPI_INT		 //solutionMin
               , MPI_INT		 //solutionMax
               , MPI_INT		 //facewiseAugmentationStatus
               , MPI_INT		 //augmentationStatus
               , MPI_INT		 //previousAugmentationStatus
               , MPI_INT		 //facewiseCommunicationStatus
               , MPI_INT		 //communicationStatus
               , MPI_INT		 //facewiseRefinementStatus
               , MPI_CXX_BOOL		 //refinementFlag
               , MPI_INT		 //refinementStatus
               , MPI_INT		 //previousRefinementStatus
               , MPI_INT		 //iterationsToCureTroubledCell
               , MPI_INT		 //_packedRecords0
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , 1		 //hasToHoldDataForMasterWorkerCommunication
               , DIMENSIONS_TIMES_TWO		 //faceDataExchangeCounter
               , 1		 //parentIndex
               , 1		 //parentCellLevel
               , DIMENSIONS		 //parentOffset
               , 1		 //hasVirtualChildren
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //previousCorrectorTimeStamp
               , 1		 //previousCorrectorTimeStepSize
               , 1		 //correctorTimeStepSize
               , 1		 //correctorTimeStamp
               , 1		 //predictorTimeStepSize
               , 1		 //predictorTimeStamp
               , 1		 //solution
               , 1		 //solutionAverages
               , 1		 //solutionCompressed
               , 1		 //previousSolution
               , 1		 //previousSolutionAverages
               , 1		 //previousSolutionCompressed
               , 1		 //update
               , 1		 //updateAverages
               , 1		 //updateCompressed
               , 1		 //extrapolatedPredictor
               , 1		 //extrapolatedPredictorAverages
               , 1		 //extrapolatedPredictorCompressed
               , 1		 //fluctuation
               , 1		 //fluctuationAverages
               , 1		 //fluctuationCompressed
               , 1		 //solutionMin
               , 1		 //solutionMax
               , DIMENSIONS_TIMES_TWO		 //facewiseAugmentationStatus
               , 1		 //augmentationStatus
               , 1		 //previousAugmentationStatus
               , DIMENSIONS_TIMES_TWO		 //facewiseCommunicationStatus
               , 1		 //communicationStatus
               , DIMENSIONS_TIMES_TWO		 //facewiseRefinementStatus
               , 1		 //refinementFlag
               , 1		 //refinementStatus
               , 1		 //previousRefinementStatus
               , 1		 //iterationsToCureTroubledCell
               , 1		 //_packedRecords0
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[23] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[24] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[24] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[25] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[25] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[26] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[26] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[27] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[27] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[28] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[28] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[29] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[29] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[30] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[30] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[31] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[31] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[32] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[32] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[33] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[33] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[34] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[34] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[35] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[35] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[36] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[36] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[37] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[37] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[38] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[38] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[39] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[39] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[40] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[40] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[41] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[41] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[42] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[42] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[43] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[43] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescriptionPacked)), i, disp[i], Attributes, sizeof(ADERDGCellDescriptionPacked));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]))), 		&disp[44] );
            disp[44] -= base;
            disp[44] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescriptionPacked::FullDatatype );
            MPI_Type_commit( &ADERDGCellDescriptionPacked::FullDatatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::FullDatatype);
            MPI_Type_commit( &ADERDGCellDescriptionPacked::FullDatatype );
            #endif
            
         }
         
      }
      
      
      void exahype::records::ADERDGCellDescriptionPacked::shutdownDatatype() {
         MPI_Type_free( &ADERDGCellDescriptionPacked::Datatype );
         MPI_Type_free( &ADERDGCellDescriptionPacked::FullDatatype );
         
      }
      
      void exahype::records::ADERDGCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
         // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked " 
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result); 
         _log.error( "send(int)",msg.str() ); 
       } 
    } 
    break; 
   case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    {
      MPI_Request* sendRequestHandle = new MPI_Request(); 
      int          flag = 0; 
       int          result; 
       clock_t      timeOutWarning   = -1; 
       clock_t      timeOutShutdown  = -1; 
       bool         triggeredTimeoutWarning = false;  
       result = MPI_Isend(  
         this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination,  
         tag, tarch::parallel::Node::getInstance().getCommunicator(), 
         sendRequestHandle  
       ); 
       if  (result!=MPI_SUCCESS) {  
         std::ostringstream msg;  
         msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked "  
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result);  
         _log.error( "send(int)",msg.str() );  
       }  
       result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
       while (!flag) { 
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
         result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
           std::ostringstream msg; 
           msg << "testing for finished send task for exahype::records::ADERDGCellDescriptionPacked " 
               << toString() 
               << " sent to node " << destination 
               << " failed: " << tarch::parallel::MPIReturnValueToString(result); 
           _log.error("send(int)", msg.str() ); 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
           (clock()>timeOutWarning) && 
           (!triggeredTimeoutWarning) 
         ) { 
           tarch::parallel::Node::getInstance().writeTimeOutWarning( 
             "exahype::records::ADERDGCellDescriptionPacked", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::ADERDGCellDescriptionPacked", 
             "send(int)", destination,tag,1 
           ); 
         } 
 	       tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
       } 
       delete sendRequestHandle; 
     }  
     break; 
   case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    assertionMsg(false,"should not be called"); 
    break; 
} 
 // ============================= 
// end injected snippet/aspect 
// ============================= 

         
      }
      
      
      
      void exahype::records::ADERDGCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
         // ============================= 
// start injected snippet/aspect 
// ============================= 
MPI_Status status; 
switch (mode) { 
  case ExchangeMode::Blocking: 
    { 
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    } 
    break; 
  case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    { 
      int          flag = 0; 
      int          result; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      MPI_Request* sendRequestHandle = new MPI_Request(); 
       result = MPI_Irecv( 
        this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, 
        tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle 
      ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node " 
             << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
      result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      delete sendRequestHandle; 
    }    break; 
  case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    {
      int flag; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      int result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
       if (result!=MPI_SUCCESS) { 
        std::ostringstream msg; 
        msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: " 
            << tarch::parallel::MPIReturnValueToString(result); 
        _log.error("receive(int)", msg.str() ); 
      } 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    }
    break; 
  } 
// =========================== 
// end injected snippet/aspect 
// =========================== 

         
      }
      
      
      
      bool exahype::records::ADERDGCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
         MPI_Status status;
         int  flag        = 0;
         MPI_Iprobe(
            MPI_ANY_SOURCE, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
         );
         if (flag) {
            int  messageCounter;
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               MPI_Get_count(&status, Datatype, &messageCounter);
            }
            else {
               MPI_Get_count(&status, FullDatatype, &messageCounter);
            }
            return messageCounter > 0;
         }
         else return false;
         
      }
      
      
   #endif
   
   
   #elif !defined(Parallel)
      exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords() {
         
      }
      
      
      exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
      _solverNumber(solverNumber),
      _neighbourMergePerformed(neighbourMergePerformed),
      _parentIndex(parentIndex),
      _parentCellLevel(parentCellLevel),
      _parentOffset(parentOffset),
      _hasVirtualChildren(hasVirtualChildren),
      _type(type),
      _refinementEvent(refinementEvent),
      _level(level),
      _offset(offset),
      _size(size),
      _previousCorrectorTimeStamp(previousCorrectorTimeStamp),
      _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
      _correctorTimeStepSize(correctorTimeStepSize),
      _correctorTimeStamp(correctorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _predictorTimeStamp(predictorTimeStamp),
      _solution(solution),
      _solutionAverages(solutionAverages),
      _solutionCompressed(solutionCompressed),
      _previousSolution(previousSolution),
      _previousSolutionAverages(previousSolutionAverages),
      _previousSolutionCompressed(previousSolutionCompressed),
      _update(update),
      _updateAverages(updateAverages),
      _updateCompressed(updateCompressed),
      _extrapolatedPredictor(extrapolatedPredictor),
      _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
      _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
      _fluctuation(fluctuation),
      _fluctuationAverages(fluctuationAverages),
      _fluctuationCompressed(fluctuationCompressed),
      _solutionMin(solutionMin),
      _solutionMax(solutionMax),
      _facewiseAugmentationStatus(facewiseAugmentationStatus),
      _augmentationStatus(augmentationStatus),
      _previousAugmentationStatus(previousAugmentationStatus),
      _facewiseCommunicationStatus(facewiseCommunicationStatus),
      _communicationStatus(communicationStatus),
      _facewiseRefinementStatus(facewiseRefinementStatus),
      _refinementFlag(refinementFlag),
      _refinementStatus(refinementStatus),
      _previousRefinementStatus(previousRefinementStatus),
      _iterationsToCureTroubledCell(iterationsToCureTroubledCell),
      _compressionState(compressionState),
      _bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
      _bytesPerDoFInSolution(bytesPerDoFInSolution),
      _bytesPerDoFInUpdate(bytesPerDoFInUpdate),
      _bytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor),
      _bytesPerDoFInFluctuation(bytesPerDoFInFluctuation) {
         
      }
      
      exahype::records::ADERDGCellDescription::ADERDGCellDescription() {
         
      }
      
      
      exahype::records::ADERDGCellDescription::ADERDGCellDescription(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords._parentIndex, persistentRecords._parentCellLevel, persistentRecords._parentOffset, persistentRecords._hasVirtualChildren, persistentRecords._type, persistentRecords._refinementEvent, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._previousCorrectorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._facewiseAugmentationStatus, persistentRecords._augmentationStatus, persistentRecords._previousAugmentationStatus, persistentRecords._facewiseCommunicationStatus, persistentRecords._communicationStatus, persistentRecords._facewiseRefinementStatus, persistentRecords._refinementFlag, persistentRecords._refinementStatus, persistentRecords._previousRefinementStatus, persistentRecords._iterationsToCureTroubledCell, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInUpdate, persistentRecords._bytesPerDoFInExtrapolatedPredictor, persistentRecords._bytesPerDoFInFluctuation) {
         
      }
      
      
      exahype::records::ADERDGCellDescription::ADERDGCellDescription(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
      _persistentRecords(solverNumber, neighbourMergePerformed, parentIndex, parentCellLevel, parentOffset, hasVirtualChildren, type, refinementEvent, level, offset, size, previousCorrectorTimeStamp, previousCorrectorTimeStepSize, correctorTimeStepSize, correctorTimeStamp, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, facewiseAugmentationStatus, augmentationStatus, previousAugmentationStatus, facewiseCommunicationStatus, communicationStatus, facewiseRefinementStatus, refinementFlag, refinementStatus, previousRefinementStatus, iterationsToCureTroubledCell, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
         
      }
      
      
      exahype::records::ADERDGCellDescription::~ADERDGCellDescription() { }
      
      std::string exahype::records::ADERDGCellDescription::toString(const CompressionState& param) {
         switch (param) {
            case Uncompressed: return "Uncompressed";
            case CurrentlyProcessed: return "CurrentlyProcessed";
            case Compressed: return "Compressed";
         }
         return "undefined";
      }
      
      std::string exahype::records::ADERDGCellDescription::getCompressionStateMapping() {
         return "CompressionState(Uncompressed=0,CurrentlyProcessed=1,Compressed=2)";
      }
      std::string exahype::records::ADERDGCellDescription::toString(const RefinementEvent& param) {
         switch (param) {
            case None: return "None";
            case ErasingChildrenRequested: return "ErasingChildrenRequested";
            case ErasingChildren: return "ErasingChildren";
            case ChangeChildrenToVirtualChildrenRequested: return "ChangeChildrenToVirtualChildrenRequested";
            case ChangeChildrenToVirtualChildren: return "ChangeChildrenToVirtualChildren";
            case RefiningRequested: return "RefiningRequested";
            case Refining: return "Refining";
            case Prolongating: return "Prolongating";
            case ErasingVirtualChildrenRequested: return "ErasingVirtualChildrenRequested";
            case ErasingVirtualChildren: return "ErasingVirtualChildren";
            case VirtualRefiningRequested: return "VirtualRefiningRequested";
            case VirtualRefining: return "VirtualRefining";
            case Erasing: return "Erasing";
            case ChangeToVirtualCell: return "ChangeToVirtualCell";
            case ErasingVirtualCell: return "ErasingVirtualCell";
         }
         return "undefined";
      }
      
      std::string exahype::records::ADERDGCellDescription::getRefinementEventMapping() {
         return "RefinementEvent(None=0,ErasingChildrenRequested=1,ErasingChildren=2,ChangeChildrenToVirtualChildrenRequested=3,ChangeChildrenToVirtualChildren=4,RefiningRequested=5,Refining=6,Prolongating=7,ErasingVirtualChildrenRequested=8,ErasingVirtualChildren=9,VirtualRefiningRequested=10,VirtualRefining=11,Erasing=12,ChangeToVirtualCell=13,ErasingVirtualCell=14)";
      }
      std::string exahype::records::ADERDGCellDescription::toString(const Type& param) {
         switch (param) {
            case Erased: return "Erased";
            case Ancestor: return "Ancestor";
            case Cell: return "Cell";
            case Descendant: return "Descendant";
         }
         return "undefined";
      }
      
      std::string exahype::records::ADERDGCellDescription::getTypeMapping() {
         return "Type(Erased=0,Ancestor=1,Cell=2,Descendant=3)";
      }
      
      
      std::string exahype::records::ADERDGCellDescription::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::ADERDGCellDescription::toString (std::ostream& out) const {
         out << "("; 
         out << "solverNumber:" << getSolverNumber();
         out << ",";
         out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "parentIndex:" << getParentIndex();
         out << ",";
         out << "parentCellLevel:" << getParentCellLevel();
         out << ",";
         out << "parentOffset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getParentOffset(i) << ",";
   }
   out << getParentOffset(DIMENSIONS-1) << "]";
         out << ",";
         out << "hasVirtualChildren:" << getHasVirtualChildren();
         out << ",";
         out << "type:" << toString(getType());
         out << ",";
         out << "refinementEvent:" << toString(getRefinementEvent());
         out << ",";
         out << "level:" << getLevel();
         out << ",";
         out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
         out << ",";
         out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
         out << ",";
         out << "previousCorrectorTimeStamp:" << getPreviousCorrectorTimeStamp();
         out << ",";
         out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
         out << ",";
         out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
         out << ",";
         out << "correctorTimeStamp:" << getCorrectorTimeStamp();
         out << ",";
         out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
         out << ",";
         out << "predictorTimeStamp:" << getPredictorTimeStamp();
         out << ",";
         out << "solution:" << getSolution();
         out << ",";
         out << "solutionAverages:" << getSolutionAverages();
         out << ",";
         out << "solutionCompressed:" << getSolutionCompressed();
         out << ",";
         out << "previousSolution:" << getPreviousSolution();
         out << ",";
         out << "previousSolutionAverages:" << getPreviousSolutionAverages();
         out << ",";
         out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
         out << ",";
         out << "update:" << getUpdate();
         out << ",";
         out << "updateAverages:" << getUpdateAverages();
         out << ",";
         out << "updateCompressed:" << getUpdateCompressed();
         out << ",";
         out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
         out << ",";
         out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
         out << ",";
         out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
         out << ",";
         out << "fluctuation:" << getFluctuation();
         out << ",";
         out << "fluctuationAverages:" << getFluctuationAverages();
         out << ",";
         out << "fluctuationCompressed:" << getFluctuationCompressed();
         out << ",";
         out << "solutionMin:" << getSolutionMin();
         out << ",";
         out << "solutionMax:" << getSolutionMax();
         out << ",";
         out << "facewiseAugmentationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseAugmentationStatus(i) << ",";
   }
   out << getFacewiseAugmentationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "augmentationStatus:" << getAugmentationStatus();
         out << ",";
         out << "previousAugmentationStatus:" << getPreviousAugmentationStatus();
         out << ",";
         out << "facewiseCommunicationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseCommunicationStatus(i) << ",";
   }
   out << getFacewiseCommunicationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "communicationStatus:" << getCommunicationStatus();
         out << ",";
         out << "facewiseRefinementStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseRefinementStatus(i) << ",";
   }
   out << getFacewiseRefinementStatus(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "refinementFlag:" << getRefinementFlag();
         out << ",";
         out << "refinementStatus:" << getRefinementStatus();
         out << ",";
         out << "previousRefinementStatus:" << getPreviousRefinementStatus();
         out << ",";
         out << "iterationsToCureTroubledCell:" << getIterationsToCureTroubledCell();
         out << ",";
         out << "compressionState:" << toString(getCompressionState());
         out << ",";
         out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
         out << ",";
         out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
         out << ",";
         out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
         out << ",";
         out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
         out << ",";
         out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
         out <<  ")";
      }
      
      
      exahype::records::ADERDGCellDescription::PersistentRecords exahype::records::ADERDGCellDescription::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::ADERDGCellDescriptionPacked exahype::records::ADERDGCellDescription::convert() const{
         return ADERDGCellDescriptionPacked(
            getSolverNumber(),
            getNeighbourMergePerformed(),
            getParentIndex(),
            getParentCellLevel(),
            getParentOffset(),
            getHasVirtualChildren(),
            getType(),
            getRefinementEvent(),
            getLevel(),
            getOffset(),
            getSize(),
            getPreviousCorrectorTimeStamp(),
            getPreviousCorrectorTimeStepSize(),
            getCorrectorTimeStepSize(),
            getCorrectorTimeStamp(),
            getPredictorTimeStepSize(),
            getPredictorTimeStamp(),
            getSolution(),
            getSolutionAverages(),
            getSolutionCompressed(),
            getPreviousSolution(),
            getPreviousSolutionAverages(),
            getPreviousSolutionCompressed(),
            getUpdate(),
            getUpdateAverages(),
            getUpdateCompressed(),
            getExtrapolatedPredictor(),
            getExtrapolatedPredictorAverages(),
            getExtrapolatedPredictorCompressed(),
            getFluctuation(),
            getFluctuationAverages(),
            getFluctuationCompressed(),
            getSolutionMin(),
            getSolutionMax(),
            getFacewiseAugmentationStatus(),
            getAugmentationStatus(),
            getPreviousAugmentationStatus(),
            getFacewiseCommunicationStatus(),
            getCommunicationStatus(),
            getFacewiseRefinementStatus(),
            getRefinementFlag(),
            getRefinementStatus(),
            getPreviousRefinementStatus(),
            getIterationsToCureTroubledCell(),
            getCompressionState(),
            getBytesPerDoFInPreviousSolution(),
            getBytesPerDoFInSolution(),
            getBytesPerDoFInUpdate(),
            getBytesPerDoFInExtrapolatedPredictor(),
            getBytesPerDoFInFluctuation()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::ADERDGCellDescription::_log( "exahype::records::ADERDGCellDescription" );
         
         MPI_Datatype exahype::records::ADERDGCellDescription::Datatype = 0;
         MPI_Datatype exahype::records::ADERDGCellDescription::FullDatatype = 0;
         
         
         void exahype::records::ADERDGCellDescription::initDatatype() {
            {
               ADERDGCellDescription dummyADERDGCellDescription[2];
               
               #ifdef MPI2
               const int Attributes = 27;
               #else
               const int Attributes = 28;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_INT		 //parentCellLevel
                  , MPI_DOUBLE		 //parentOffset
                  , MPI_CXX_BOOL		 //hasVirtualChildren
                  , MPI_INT		 //type
                  , MPI_INT		 //refinementEvent
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_DOUBLE		 //previousCorrectorTimeStamp
                  , MPI_DOUBLE		 //previousCorrectorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStamp
                  , MPI_DOUBLE		 //predictorTimeStepSize
                  , MPI_DOUBLE		 //predictorTimeStamp
                  , MPI_INT		 //augmentationStatus
                  , MPI_INT		 //previousAugmentationStatus
                  , MPI_INT		 //communicationStatus
                  , MPI_INT		 //refinementStatus
                  , MPI_INT		 //previousRefinementStatus
                  , MPI_INT		 //iterationsToCureTroubledCell
                  , MPI_INT		 //compressionState
                  , MPI_INT		 //bytesPerDoFInPreviousSolution
                  , MPI_INT		 //bytesPerDoFInSolution
                  , MPI_INT		 //bytesPerDoFInUpdate
                  , MPI_INT		 //bytesPerDoFInExtrapolatedPredictor
                  , MPI_INT		 //bytesPerDoFInFluctuation
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , 1		 //parentCellLevel
                  , DIMENSIONS		 //parentOffset
                  , 1		 //hasVirtualChildren
                  , 1		 //type
                  , 1		 //refinementEvent
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //previousCorrectorTimeStamp
                  , 1		 //previousCorrectorTimeStepSize
                  , 1		 //correctorTimeStepSize
                  , 1		 //correctorTimeStamp
                  , 1		 //predictorTimeStepSize
                  , 1		 //predictorTimeStamp
                  , 1		 //augmentationStatus
                  , 1		 //previousAugmentationStatus
                  , 1		 //communicationStatus
                  , 1		 //refinementStatus
                  , 1		 //previousRefinementStatus
                  , 1		 //iterationsToCureTroubledCell
                  , 1		 //compressionState
                  , 1		 //bytesPerDoFInPreviousSolution
                  , 1		 //bytesPerDoFInSolution
                  , 1		 //bytesPerDoFInUpdate
                  , 1		 //bytesPerDoFInExtrapolatedPredictor
                  , 1		 //bytesPerDoFInFluctuation
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[24] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[25] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[25] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[26] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[26] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescription)), i, disp[i], Attributes, sizeof(ADERDGCellDescription));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]))), 		&disp[27] );
               disp[27] -= base;
               disp[27] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescription::Datatype );
               MPI_Type_commit( &ADERDGCellDescription::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::Datatype);
               MPI_Type_commit( &ADERDGCellDescription::Datatype );
               #endif
               
            }
            {
               ADERDGCellDescription dummyADERDGCellDescription[2];
               
               #ifdef MPI2
               const int Attributes = 50;
               #else
               const int Attributes = 51;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_CXX_BOOL		 //neighbourMergePerformed
                  , MPI_INT		 //parentIndex
                  , MPI_INT		 //parentCellLevel
                  , MPI_DOUBLE		 //parentOffset
                  , MPI_CXX_BOOL		 //hasVirtualChildren
                  , MPI_INT		 //type
                  , MPI_INT		 //refinementEvent
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_DOUBLE		 //previousCorrectorTimeStamp
                  , MPI_DOUBLE		 //previousCorrectorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStamp
                  , MPI_DOUBLE		 //predictorTimeStepSize
                  , MPI_DOUBLE		 //predictorTimeStamp
                  , MPI_INT		 //solution
                  , MPI_INT		 //solutionAverages
                  , MPI_INT		 //solutionCompressed
                  , MPI_INT		 //previousSolution
                  , MPI_INT		 //previousSolutionAverages
                  , MPI_INT		 //previousSolutionCompressed
                  , MPI_INT		 //update
                  , MPI_INT		 //updateAverages
                  , MPI_INT		 //updateCompressed
                  , MPI_INT		 //extrapolatedPredictor
                  , MPI_INT		 //extrapolatedPredictorAverages
                  , MPI_INT		 //extrapolatedPredictorCompressed
                  , MPI_INT		 //fluctuation
                  , MPI_INT		 //fluctuationAverages
                  , MPI_INT		 //fluctuationCompressed
                  , MPI_INT		 //solutionMin
                  , MPI_INT		 //solutionMax
                  , MPI_INT		 //facewiseAugmentationStatus
                  , MPI_INT		 //augmentationStatus
                  , MPI_INT		 //previousAugmentationStatus
                  , MPI_INT		 //facewiseCommunicationStatus
                  , MPI_INT		 //communicationStatus
                  , MPI_INT		 //facewiseRefinementStatus
                  , MPI_CXX_BOOL		 //refinementFlag
                  , MPI_INT		 //refinementStatus
                  , MPI_INT		 //previousRefinementStatus
                  , MPI_INT		 //iterationsToCureTroubledCell
                  , MPI_INT		 //compressionState
                  , MPI_INT		 //bytesPerDoFInPreviousSolution
                  , MPI_INT		 //bytesPerDoFInSolution
                  , MPI_INT		 //bytesPerDoFInUpdate
                  , MPI_INT		 //bytesPerDoFInExtrapolatedPredictor
                  , MPI_INT		 //bytesPerDoFInFluctuation
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
                  , 1		 //parentIndex
                  , 1		 //parentCellLevel
                  , DIMENSIONS		 //parentOffset
                  , 1		 //hasVirtualChildren
                  , 1		 //type
                  , 1		 //refinementEvent
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //previousCorrectorTimeStamp
                  , 1		 //previousCorrectorTimeStepSize
                  , 1		 //correctorTimeStepSize
                  , 1		 //correctorTimeStamp
                  , 1		 //predictorTimeStepSize
                  , 1		 //predictorTimeStamp
                  , 1		 //solution
                  , 1		 //solutionAverages
                  , 1		 //solutionCompressed
                  , 1		 //previousSolution
                  , 1		 //previousSolutionAverages
                  , 1		 //previousSolutionCompressed
                  , 1		 //update
                  , 1		 //updateAverages
                  , 1		 //updateCompressed
                  , 1		 //extrapolatedPredictor
                  , 1		 //extrapolatedPredictorAverages
                  , 1		 //extrapolatedPredictorCompressed
                  , 1		 //fluctuation
                  , 1		 //fluctuationAverages
                  , 1		 //fluctuationCompressed
                  , 1		 //solutionMin
                  , 1		 //solutionMax
                  , DIMENSIONS_TIMES_TWO		 //facewiseAugmentationStatus
                  , 1		 //augmentationStatus
                  , 1		 //previousAugmentationStatus
                  , DIMENSIONS_TIMES_TWO		 //facewiseCommunicationStatus
                  , 1		 //communicationStatus
                  , DIMENSIONS_TIMES_TWO		 //facewiseRefinementStatus
                  , 1		 //refinementFlag
                  , 1		 //refinementStatus
                  , 1		 //previousRefinementStatus
                  , 1		 //iterationsToCureTroubledCell
                  , 1		 //compressionState
                  , 1		 //bytesPerDoFInPreviousSolution
                  , 1		 //bytesPerDoFInSolution
                  , 1		 //bytesPerDoFInUpdate
                  , 1		 //bytesPerDoFInExtrapolatedPredictor
                  , 1		 //bytesPerDoFInFluctuation
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentCellLevel))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentOffset[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasVirtualChildren))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[24] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[25] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[25] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[26] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[26] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[27] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[27] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[28] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[28] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[29] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[29] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[30] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[30] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[31] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[31] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[32] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[32] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[33] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[33] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[34] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[34] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[35] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[35] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[36] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousAugmentationStatus))), 		&disp[36] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[37] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[37] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[38] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[38] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[39] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[39] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[40] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[40] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[41] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[41] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[42] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[42] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[43] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[43] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[44] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[44] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[45] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[45] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[46] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[46] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[47] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[47] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[48] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[48] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[49] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[49] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescription)), i, disp[i], Attributes, sizeof(ADERDGCellDescription));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]))), 		&disp[50] );
               disp[50] -= base;
               disp[50] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescription::FullDatatype );
               MPI_Type_commit( &ADERDGCellDescription::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::FullDatatype);
               MPI_Type_commit( &ADERDGCellDescription::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::ADERDGCellDescription::shutdownDatatype() {
            MPI_Type_free( &ADERDGCellDescription::Datatype );
            MPI_Type_free( &ADERDGCellDescription::FullDatatype );
            
         }
         
         void exahype::records::ADERDGCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
            // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::ADERDGCellDescription " 
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result); 
         _log.error( "send(int)",msg.str() ); 
       } 
    } 
    break; 
   case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    {
      MPI_Request* sendRequestHandle = new MPI_Request(); 
      int          flag = 0; 
       int          result; 
       clock_t      timeOutWarning   = -1; 
       clock_t      timeOutShutdown  = -1; 
       bool         triggeredTimeoutWarning = false;  
       result = MPI_Isend(  
         this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination,  
         tag, tarch::parallel::Node::getInstance().getCommunicator(), 
         sendRequestHandle  
       ); 
       if  (result!=MPI_SUCCESS) {  
         std::ostringstream msg;  
         msg << "was not able to send message exahype::records::ADERDGCellDescription "  
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result);  
         _log.error( "send(int)",msg.str() );  
       }  
       result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
       while (!flag) { 
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
         result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
           std::ostringstream msg; 
           msg << "testing for finished send task for exahype::records::ADERDGCellDescription " 
               << toString() 
               << " sent to node " << destination 
               << " failed: " << tarch::parallel::MPIReturnValueToString(result); 
           _log.error("send(int)", msg.str() ); 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
           (clock()>timeOutWarning) && 
           (!triggeredTimeoutWarning) 
         ) { 
           tarch::parallel::Node::getInstance().writeTimeOutWarning( 
             "exahype::records::ADERDGCellDescription", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::ADERDGCellDescription", 
             "send(int)", destination,tag,1 
           ); 
         } 
 	       tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
       } 
       delete sendRequestHandle; 
     }  
     break; 
   case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    assertionMsg(false,"should not be called"); 
    break; 
} 
 // ============================= 
// end injected snippet/aspect 
// ============================= 

            
         }
         
         
         
         void exahype::records::ADERDGCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
            // ============================= 
// start injected snippet/aspect 
// ============================= 
MPI_Status status; 
switch (mode) { 
  case ExchangeMode::Blocking: 
    { 
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescription from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    } 
    break; 
  case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    { 
      int          flag = 0; 
      int          result; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      MPI_Request* sendRequestHandle = new MPI_Request(); 
       result = MPI_Irecv( 
        this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, 
        tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle 
      ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescription from node " 
             << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
      result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      delete sendRequestHandle; 
    }    break; 
  case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    {
      int flag; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      int result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
       if (result!=MPI_SUCCESS) { 
        std::ostringstream msg; 
        msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: " 
            << tarch::parallel::MPIReturnValueToString(result); 
        _log.error("receive(int)", msg.str() ); 
      } 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescription from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    }
    break; 
  } 
// =========================== 
// end injected snippet/aspect 
// =========================== 

            
         }
         
         
         
         bool exahype::records::ADERDGCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         
      #endif
      
      
      exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords() {
         if ((DIMENSIONS_TIMES_TWO+23 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((DIMENSIONS_TIMES_TWO+23 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
      _solverNumber(solverNumber),
      _parentIndex(parentIndex),
      _parentCellLevel(parentCellLevel),
      _parentOffset(parentOffset),
      _hasVirtualChildren(hasVirtualChildren),
      _level(level),
      _offset(offset),
      _size(size),
      _previousCorrectorTimeStamp(previousCorrectorTimeStamp),
      _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
      _correctorTimeStepSize(correctorTimeStepSize),
      _correctorTimeStamp(correctorTimeStamp),
      _predictorTimeStepSize(predictorTimeStepSize),
      _predictorTimeStamp(predictorTimeStamp),
      _solution(solution),
      _solutionAverages(solutionAverages),
      _solutionCompressed(solutionCompressed),
      _previousSolution(previousSolution),
      _previousSolutionAverages(previousSolutionAverages),
      _previousSolutionCompressed(previousSolutionCompressed),
      _update(update),
      _updateAverages(updateAverages),
      _updateCompressed(updateCompressed),
      _extrapolatedPredictor(extrapolatedPredictor),
      _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
      _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
      _fluctuation(fluctuation),
      _fluctuationAverages(fluctuationAverages),
      _fluctuationCompressed(fluctuationCompressed),
      _solutionMin(solutionMin),
      _solutionMax(solutionMax),
      _facewiseAugmentationStatus(facewiseAugmentationStatus),
      _augmentationStatus(augmentationStatus),
      _previousAugmentationStatus(previousAugmentationStatus),
      _facewiseCommunicationStatus(facewiseCommunicationStatus),
      _communicationStatus(communicationStatus),
      _facewiseRefinementStatus(facewiseRefinementStatus),
      _refinementFlag(refinementFlag),
      _refinementStatus(refinementStatus),
      _previousRefinementStatus(previousRefinementStatus),
      _iterationsToCureTroubledCell(iterationsToCureTroubledCell) {
         setNeighbourMergePerformed(neighbourMergePerformed);
         setType(type);
         setRefinementEvent(refinementEvent);
         setCompressionState(compressionState);
         setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
         setBytesPerDoFInSolution(bytesPerDoFInSolution);
         setBytesPerDoFInUpdate(bytesPerDoFInUpdate);
         setBytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor);
         setBytesPerDoFInFluctuation(bytesPerDoFInFluctuation);
         if ((DIMENSIONS_TIMES_TWO+23 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((DIMENSIONS_TIMES_TWO+23 < (8 * sizeof(int))));
         
      }
      
      exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked() {
         if ((DIMENSIONS_TIMES_TWO+23 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((DIMENSIONS_TIMES_TWO+23 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._solverNumber, persistentRecords.getNeighbourMergePerformed(), persistentRecords._parentIndex, persistentRecords._parentCellLevel, persistentRecords._parentOffset, persistentRecords._hasVirtualChildren, persistentRecords.getType(), persistentRecords.getRefinementEvent(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._previousCorrectorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._facewiseAugmentationStatus, persistentRecords._augmentationStatus, persistentRecords._previousAugmentationStatus, persistentRecords._facewiseCommunicationStatus, persistentRecords._communicationStatus, persistentRecords._facewiseRefinementStatus, persistentRecords._refinementFlag, persistentRecords._refinementStatus, persistentRecords._previousRefinementStatus, persistentRecords._iterationsToCureTroubledCell, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInUpdate(), persistentRecords.getBytesPerDoFInExtrapolatedPredictor(), persistentRecords.getBytesPerDoFInFluctuation()) {
         if ((DIMENSIONS_TIMES_TWO+23 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((DIMENSIONS_TIMES_TWO+23 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const int& parentIndex, const int& parentCellLevel, const tarch::la::Vector<DIMENSIONS,double>& parentOffset, const bool& hasVirtualChildren, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousCorrectorTimeStamp, const double& previousCorrectorTimeStepSize, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const int& previousAugmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const bool& refinementFlag, const int& refinementStatus, const int& previousRefinementStatus, const int& iterationsToCureTroubledCell, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
      _persistentRecords(solverNumber, neighbourMergePerformed, parentIndex, parentCellLevel, parentOffset, hasVirtualChildren, type, refinementEvent, level, offset, size, previousCorrectorTimeStamp, previousCorrectorTimeStepSize, correctorTimeStepSize, correctorTimeStamp, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, facewiseAugmentationStatus, augmentationStatus, previousAugmentationStatus, facewiseCommunicationStatus, communicationStatus, facewiseRefinementStatus, refinementFlag, refinementStatus, previousRefinementStatus, iterationsToCureTroubledCell, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
         if ((DIMENSIONS_TIMES_TWO+23 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((DIMENSIONS_TIMES_TWO+23 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::ADERDGCellDescriptionPacked::~ADERDGCellDescriptionPacked() { }
      
      std::string exahype::records::ADERDGCellDescriptionPacked::toString(const Type& param) {
         return exahype::records::ADERDGCellDescription::toString(param);
      }
      
      std::string exahype::records::ADERDGCellDescriptionPacked::getTypeMapping() {
         return exahype::records::ADERDGCellDescription::getTypeMapping();
      }
      
      std::string exahype::records::ADERDGCellDescriptionPacked::toString(const RefinementEvent& param) {
         return exahype::records::ADERDGCellDescription::toString(param);
      }
      
      std::string exahype::records::ADERDGCellDescriptionPacked::getRefinementEventMapping() {
         return exahype::records::ADERDGCellDescription::getRefinementEventMapping();
      }
      
      std::string exahype::records::ADERDGCellDescriptionPacked::toString(const CompressionState& param) {
         return exahype::records::ADERDGCellDescription::toString(param);
      }
      
      std::string exahype::records::ADERDGCellDescriptionPacked::getCompressionStateMapping() {
         return exahype::records::ADERDGCellDescription::getCompressionStateMapping();
      }
      
      
      
      std::string exahype::records::ADERDGCellDescriptionPacked::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::ADERDGCellDescriptionPacked::toString (std::ostream& out) const {
         out << "("; 
         out << "solverNumber:" << getSolverNumber();
         out << ",";
         out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "parentIndex:" << getParentIndex();
         out << ",";
         out << "parentCellLevel:" << getParentCellLevel();
         out << ",";
         out << "parentOffset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getParentOffset(i) << ",";
   }
   out << getParentOffset(DIMENSIONS-1) << "]";
         out << ",";
         out << "hasVirtualChildren:" << getHasVirtualChildren();
         out << ",";
         out << "type:" << toString(getType());
         out << ",";
         out << "refinementEvent:" << toString(getRefinementEvent());
         out << ",";
         out << "level:" << getLevel();
         out << ",";
         out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
         out << ",";
         out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
         out << ",";
         out << "previousCorrectorTimeStamp:" << getPreviousCorrectorTimeStamp();
         out << ",";
         out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
         out << ",";
         out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
         out << ",";
         out << "correctorTimeStamp:" << getCorrectorTimeStamp();
         out << ",";
         out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
         out << ",";
         out << "predictorTimeStamp:" << getPredictorTimeStamp();
         out << ",";
         out << "solution:" << getSolution();
         out << ",";
         out << "solutionAverages:" << getSolutionAverages();
         out << ",";
         out << "solutionCompressed:" << getSolutionCompressed();
         out << ",";
         out << "previousSolution:" << getPreviousSolution();
         out << ",";
         out << "previousSolutionAverages:" << getPreviousSolutionAverages();
         out << ",";
         out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
         out << ",";
         out << "update:" << getUpdate();
         out << ",";
         out << "updateAverages:" << getUpdateAverages();
         out << ",";
         out << "updateCompressed:" << getUpdateCompressed();
         out << ",";
         out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
         out << ",";
         out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
         out << ",";
         out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
         out << ",";
         out << "fluctuation:" << getFluctuation();
         out << ",";
         out << "fluctuationAverages:" << getFluctuationAverages();
         out << ",";
         out << "fluctuationCompressed:" << getFluctuationCompressed();
         out << ",";
         out << "solutionMin:" << getSolutionMin();
         out << ",";
         out << "solutionMax:" << getSolutionMax();
         out << ",";
         out << "facewiseAugmentationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseAugmentationStatus(i) << ",";
   }
   out << getFacewiseAugmentationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "augmentationStatus:" << getAugmentationStatus();
         out << ",";
         out << "previousAugmentationStatus:" << getPreviousAugmentationStatus();
         out << ",";
         out << "facewiseCommunicationStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseCommunicationStatus(i) << ",";
   }
   out << getFacewiseCommunicationStatus(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "communicationStatus:" << getCommunicationStatus();
         out << ",";
         out << "facewiseRefinementStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFacewiseRefinementStatus(i) << ",";
   }
   out << getFacewiseRefinementStatus(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "refinementFlag:" << getRefinementFlag();
         out << ",";
         out << "refinementStatus:" << getRefinementStatus();
         out << ",";
         out << "previousRefinementStatus:" << getPreviousRefinementStatus();
         out << ",";
         out << "iterationsToCureTroubledCell:" << getIterationsToCureTroubledCell();
         out << ",";
         out << "compressionState:" << toString(getCompressionState());
         out << ",";
         out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
         out << ",";
         out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
         out << ",";
         out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
         out << ",";
         out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
         out << ",";
         out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
         out <<  ")";
      }
      
      
      exahype::records::ADERDGCellDescriptionPacked::PersistentRecords exahype::records::ADERDGCellDescriptionPacked::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::ADERDGCellDescription exahype::records::ADERDGCellDescriptionPacked::convert() const{
         return ADERDGCellDescription(
            getSolverNumber(),
            getNeighbourMergePerformed(),
            getParentIndex(),
            getParentCellLevel(),
            getParentOffset(),
            getHasVirtualChildren(),
            getType(),
            getRefinementEvent(),
            getLevel(),
            getOffset(),
            getSize(),
            getPreviousCorrectorTimeStamp(),
            getPreviousCorrectorTimeStepSize(),
            getCorrectorTimeStepSize(),
            getCorrectorTimeStamp(),
            getPredictorTimeStepSize(),
            getPredictorTimeStamp(),
            getSolution(),
            getSolutionAverages(),
            getSolutionCompressed(),
            getPreviousSolution(),
            getPreviousSolutionAverages(),
            getPreviousSolutionCompressed(),
            getUpdate(),
            getUpdateAverages(),
            getUpdateCompressed(),
            getExtrapolatedPredictor(),
            getExtrapolatedPredictorAverages(),
            getExtrapolatedPredictorCompressed(),
            getFluctuation(),
            getFluctuationAverages(),
            getFluctuationCompressed(),
            getSolutionMin(),
            getSolutionMax(),
            getFacewiseAugmentationStatus(),
            getAugmentationStatus(),
            getPreviousAugmentationStatus(),
            getFacewiseCommunicationStatus(),
            getCommunicationStatus(),
            getFacewiseRefinementStatus(),
            getRefinementFlag(),
            getRefinementStatus(),
            getPreviousRefinementStatus(),
            getIterationsToCureTroubledCell(),
            getCompressionState(),
            getBytesPerDoFInPreviousSolution(),
            getBytesPerDoFInSolution(),
            getBytesPerDoFInUpdate(),
            getBytesPerDoFInExtrapolatedPredictor(),
            getBytesPerDoFInFluctuation()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::ADERDGCellDescriptionPacked::_log( "exahype::records::ADERDGCellDescriptionPacked" );
         
         MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::Datatype = 0;
         MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::FullDatatype = 0;
         
         
         void exahype::records::ADERDGCellDescriptionPacked::initDatatype() {
            {
               ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
               
               #ifdef MPI2
               const int Attributes = 20;
               #else
               const int Attributes = 21;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_INT		 //parentCellLevel
                  , MPI_DOUBLE		 //parentOffset
                  , MPI_CXX_BOOL		 //hasVirtualChildren
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_DOUBLE		 //previousCorrectorTimeStamp
                  , MPI_DOUBLE		 //previousCorrectorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStamp
                  , MPI_DOUBLE		 //predictorTimeStepSize
                  , MPI_DOUBLE		 //predictorTimeStamp
                  , MPI_INT		 //augmentationStatus
                  , MPI_INT		 //previousAugmentationStatus
                  , MPI_INT		 //communicationStatus
                  , MPI_INT		 //refinementStatus
                  , MPI_INT		 //previousRefinementStatus
                  , MPI_INT		 //iterationsToCureTroubledCell
                  , MPI_INT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , 1		 //parentCellLevel
                  , DIMENSIONS		 //parentOffset
                  , 1		 //hasVirtualChildren
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //previousCorrectorTimeStamp
                  , 1		 //previousCorrectorTimeStepSize
                  , 1		 //correctorTimeStepSize
                  , 1		 //correctorTimeStamp
                  , 1		 //predictorTimeStepSize
                  , 1		 //predictorTimeStamp
                  , 1		 //augmentationStatus
                  , 1		 //previousAugmentationStatus
                  , 1		 //communicationStatus
                  , 1		 //refinementStatus
                  , 1		 //previousRefinementStatus
                  , 1		 //iterationsToCureTroubledCell
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[19] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescriptionPacked)), i, disp[i], Attributes, sizeof(ADERDGCellDescriptionPacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]))), 		&disp[20] );
               disp[20] -= base;
               disp[20] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescriptionPacked::Datatype );
               MPI_Type_commit( &ADERDGCellDescriptionPacked::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::Datatype);
               MPI_Type_commit( &ADERDGCellDescriptionPacked::Datatype );
               #endif
               
            }
            {
               ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
               
               #ifdef MPI2
               const int Attributes = 42;
               #else
               const int Attributes = 43;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_INT		 //parentIndex
                  , MPI_INT		 //parentCellLevel
                  , MPI_DOUBLE		 //parentOffset
                  , MPI_CXX_BOOL		 //hasVirtualChildren
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_DOUBLE		 //previousCorrectorTimeStamp
                  , MPI_DOUBLE		 //previousCorrectorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStepSize
                  , MPI_DOUBLE		 //correctorTimeStamp
                  , MPI_DOUBLE		 //predictorTimeStepSize
                  , MPI_DOUBLE		 //predictorTimeStamp
                  , MPI_INT		 //solution
                  , MPI_INT		 //solutionAverages
                  , MPI_INT		 //solutionCompressed
                  , MPI_INT		 //previousSolution
                  , MPI_INT		 //previousSolutionAverages
                  , MPI_INT		 //previousSolutionCompressed
                  , MPI_INT		 //update
                  , MPI_INT		 //updateAverages
                  , MPI_INT		 //updateCompressed
                  , MPI_INT		 //extrapolatedPredictor
                  , MPI_INT		 //extrapolatedPredictorAverages
                  , MPI_INT		 //extrapolatedPredictorCompressed
                  , MPI_INT		 //fluctuation
                  , MPI_INT		 //fluctuationAverages
                  , MPI_INT		 //fluctuationCompressed
                  , MPI_INT		 //solutionMin
                  , MPI_INT		 //solutionMax
                  , MPI_INT		 //facewiseAugmentationStatus
                  , MPI_INT		 //augmentationStatus
                  , MPI_INT		 //previousAugmentationStatus
                  , MPI_INT		 //facewiseCommunicationStatus
                  , MPI_INT		 //communicationStatus
                  , MPI_INT		 //facewiseRefinementStatus
                  , MPI_CXX_BOOL		 //refinementFlag
                  , MPI_INT		 //refinementStatus
                  , MPI_INT		 //previousRefinementStatus
                  , MPI_INT		 //iterationsToCureTroubledCell
                  , MPI_INT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , 1		 //parentIndex
                  , 1		 //parentCellLevel
                  , DIMENSIONS		 //parentOffset
                  , 1		 //hasVirtualChildren
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //previousCorrectorTimeStamp
                  , 1		 //previousCorrectorTimeStepSize
                  , 1		 //correctorTimeStepSize
                  , 1		 //correctorTimeStamp
                  , 1		 //predictorTimeStepSize
                  , 1		 //predictorTimeStamp
                  , 1		 //solution
                  , 1		 //solutionAverages
                  , 1		 //solutionCompressed
                  , 1		 //previousSolution
                  , 1		 //previousSolutionAverages
                  , 1		 //previousSolutionCompressed
                  , 1		 //update
                  , 1		 //updateAverages
                  , 1		 //updateCompressed
                  , 1		 //extrapolatedPredictor
                  , 1		 //extrapolatedPredictorAverages
                  , 1		 //extrapolatedPredictorCompressed
                  , 1		 //fluctuation
                  , 1		 //fluctuationAverages
                  , 1		 //fluctuationCompressed
                  , 1		 //solutionMin
                  , 1		 //solutionMax
                  , DIMENSIONS_TIMES_TWO		 //facewiseAugmentationStatus
                  , 1		 //augmentationStatus
                  , 1		 //previousAugmentationStatus
                  , DIMENSIONS_TIMES_TWO		 //facewiseCommunicationStatus
                  , 1		 //communicationStatus
                  , DIMENSIONS_TIMES_TWO		 //facewiseRefinementStatus
                  , 1		 //refinementFlag
                  , 1		 //refinementStatus
                  , 1		 //previousRefinementStatus
                  , 1		 //iterationsToCureTroubledCell
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentCellLevel))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentOffset[0]))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasVirtualChildren))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStamp))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[24] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[25] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[25] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[26] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[26] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[27] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[27] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[28] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[28] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[29] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[29] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[30] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[30] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[31] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[31] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[32] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[32] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[33] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousAugmentationStatus))), 		&disp[33] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[34] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[34] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[35] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[35] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[36] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[36] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[37] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[37] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[38] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[38] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[39] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[39] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[40] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._iterationsToCureTroubledCell))), 		&disp[40] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[41] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[41] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(ADERDGCellDescriptionPacked)), i, disp[i], Attributes, sizeof(ADERDGCellDescriptionPacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]))), 		&disp[42] );
               disp[42] -= base;
               disp[42] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &ADERDGCellDescriptionPacked::FullDatatype );
               MPI_Type_commit( &ADERDGCellDescriptionPacked::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::FullDatatype);
               MPI_Type_commit( &ADERDGCellDescriptionPacked::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::ADERDGCellDescriptionPacked::shutdownDatatype() {
            MPI_Type_free( &ADERDGCellDescriptionPacked::Datatype );
            MPI_Type_free( &ADERDGCellDescriptionPacked::FullDatatype );
            
         }
         
         void exahype::records::ADERDGCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
            // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked " 
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result); 
         _log.error( "send(int)",msg.str() ); 
       } 
    } 
    break; 
   case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    {
      MPI_Request* sendRequestHandle = new MPI_Request(); 
      int          flag = 0; 
       int          result; 
       clock_t      timeOutWarning   = -1; 
       clock_t      timeOutShutdown  = -1; 
       bool         triggeredTimeoutWarning = false;  
       result = MPI_Isend(  
         this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination,  
         tag, tarch::parallel::Node::getInstance().getCommunicator(), 
         sendRequestHandle  
       ); 
       if  (result!=MPI_SUCCESS) {  
         std::ostringstream msg;  
         msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked "  
             << toString() 
             << " to node " << destination 
             << ": " << tarch::parallel::MPIReturnValueToString(result);  
         _log.error( "send(int)",msg.str() );  
       }  
       result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
       while (!flag) { 
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
         result = MPI_Test( sendRequestHandle, &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
           std::ostringstream msg; 
           msg << "testing for finished send task for exahype::records::ADERDGCellDescriptionPacked " 
               << toString() 
               << " sent to node " << destination 
               << " failed: " << tarch::parallel::MPIReturnValueToString(result); 
           _log.error("send(int)", msg.str() ); 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
           (clock()>timeOutWarning) && 
           (!triggeredTimeoutWarning) 
         ) { 
           tarch::parallel::Node::getInstance().writeTimeOutWarning( 
             "exahype::records::ADERDGCellDescriptionPacked", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::ADERDGCellDescriptionPacked", 
             "send(int)", destination,tag,1 
           ); 
         } 
 	       tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
       } 
       delete sendRequestHandle; 
     }  
     break; 
   case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    assertionMsg(false,"should not be called"); 
    break; 
} 
 // ============================= 
// end injected snippet/aspect 
// ============================= 

            
         }
         
         
         
         void exahype::records::ADERDGCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
            // ============================= 
// start injected snippet/aspect 
// ============================= 
MPI_Status status; 
switch (mode) { 
  case ExchangeMode::Blocking: 
    { 
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    } 
    break; 
  case ExchangeMode::NonblockingWithPollingLoopOverTests: 
    { 
      int          flag = 0; 
      int          result; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      MPI_Request* sendRequestHandle = new MPI_Request(); 
       result = MPI_Irecv( 
        this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, 
        tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle 
      ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node " 
             << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
      result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      delete sendRequestHandle; 
    }    break; 
  case ExchangeMode::LoopOverProbeWithBlockingReceive: 
    {
      int flag; 
      clock_t      timeOutWarning   = -1; 
      clock_t      timeOutShutdown  = -1; 
      bool         triggeredTimeoutWarning = false; 
      int result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
       if (result!=MPI_SUCCESS) { 
        std::ostringstream msg; 
        msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: " 
            << tarch::parallel::MPIReturnValueToString(result); 
        _log.error("receive(int)", msg.str() ); 
      } 
      while (!flag) { 
        if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp(); 
        if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp(); 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() && 
          (clock()>timeOutWarning) && 
          (!triggeredTimeoutWarning) 
        ) { 
          tarch::parallel::Node::getInstance().writeTimeOutWarning( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::ADERDGCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node " 
            << source << ": " << tarch::parallel::MPIReturnValueToString(result); 
        _log.error( "receive(int)", msg.str() ); 
      } 
    }
    break; 
  } 
// =========================== 
// end injected snippet/aspect 
// =========================== 

            
         }
         
         
         
         bool exahype::records::ADERDGCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         
      #endif
      
      
      
   
#endif


