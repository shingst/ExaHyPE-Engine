#include "exahype/records/ADERDGCellDescription.h"

exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords() {
   
}


exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed, const int& parentIndex, const Type& type, const Type& parentType, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousTimeStamp, const double& previousTimeStepSize, const double& timeStepSize, const double& timeStamp, const int& solutionIndex, const int& solutionAveragesIndex, const int& solutionCompressedIndex, void* solution, void* solutionAverages, void* solutionCompressed, const int& previousSolutionIndex, const int& previousSolutionAveragesIndex, const int& previousSolutionCompressedIndex, void* previousSolution, void* previousSolutionAverages, void* previousSolutionCompressed, const int& updateIndex, const int& updateAveragesIndex, const int& updateCompressedIndex, void* update, void* updateAverages, void* updateCompressed, const int& extrapolatedPredictorIndex, const int& extrapolatedPredictorAveragesIndex, const int& extrapolatedPredictorCompressedIndex, void* extrapolatedPredictor, void* extrapolatedPredictorAverages, void* extrapolatedPredictorCompressed, const int& extrapolatedPredictorGradientIndex, void* extrapolatedPredictorGradient, const int& fluctuationIndex, const int& fluctuationAveragesIndex, const int& fluctuationCompressedIndex, void* fluctuation, void* fluctuationAverages, void* fluctuationCompressed, const int& solutionMinIndex, const int& solutionMaxIndex, void* solutionMin, void* solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const int& refinementStatus, const int& previousRefinementStatus, const bool& refinementFlag, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation, const Creation& creation):
_solverNumber(solverNumber),
_neighbourMergePerformed(neighbourMergePerformed),
_parentIndex(parentIndex),
_type(type),
_parentType(parentType),
_level(level),
_offset(offset),
_size(size),
_previousTimeStamp(previousTimeStamp),
_previousTimeStepSize(previousTimeStepSize),
_timeStepSize(timeStepSize),
_timeStamp(timeStamp),
_solutionIndex(solutionIndex),
_solutionAveragesIndex(solutionAveragesIndex),
_solutionCompressedIndex(solutionCompressedIndex),
_solution(solution),
_solutionAverages(solutionAverages),
_solutionCompressed(solutionCompressed),
_previousSolutionIndex(previousSolutionIndex),
_previousSolutionAveragesIndex(previousSolutionAveragesIndex),
_previousSolutionCompressedIndex(previousSolutionCompressedIndex),
_previousSolution(previousSolution),
_previousSolutionAverages(previousSolutionAverages),
_previousSolutionCompressed(previousSolutionCompressed),
_updateIndex(updateIndex),
_updateAveragesIndex(updateAveragesIndex),
_updateCompressedIndex(updateCompressedIndex),
_update(update),
_updateAverages(updateAverages),
_updateCompressed(updateCompressed),
_extrapolatedPredictorIndex(extrapolatedPredictorIndex),
_extrapolatedPredictorAveragesIndex(extrapolatedPredictorAveragesIndex),
_extrapolatedPredictorCompressedIndex(extrapolatedPredictorCompressedIndex),
_extrapolatedPredictor(extrapolatedPredictor),
_extrapolatedPredictorAverages(extrapolatedPredictorAverages),
_extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
_extrapolatedPredictorGradientIndex(extrapolatedPredictorGradientIndex),
_extrapolatedPredictorGradient(extrapolatedPredictorGradient),
_fluctuationIndex(fluctuationIndex),
_fluctuationAveragesIndex(fluctuationAveragesIndex),
_fluctuationCompressedIndex(fluctuationCompressedIndex),
_fluctuation(fluctuation),
_fluctuationAverages(fluctuationAverages),
_fluctuationCompressed(fluctuationCompressed),
_solutionMinIndex(solutionMinIndex),
_solutionMaxIndex(solutionMaxIndex),
_solutionMin(solutionMin),
_solutionMax(solutionMax),
_facewiseAugmentationStatus(facewiseAugmentationStatus),
_augmentationStatus(augmentationStatus),
_facewiseCommunicationStatus(facewiseCommunicationStatus),
_communicationStatus(communicationStatus),
_facewiseRefinementStatus(facewiseRefinementStatus),
_refinementStatus(refinementStatus),
_previousRefinementStatus(previousRefinementStatus),
_refinementFlag(refinementFlag),
_compressionState(compressionState),
_bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
_bytesPerDoFInSolution(bytesPerDoFInSolution),
_bytesPerDoFInUpdate(bytesPerDoFInUpdate),
_bytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor),
_bytesPerDoFInFluctuation(bytesPerDoFInFluctuation),
_creation(creation) {
   
}

exahype::records::ADERDGCellDescription::ADERDGCellDescription() {
   
}


exahype::records::ADERDGCellDescription::ADERDGCellDescription(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords._parentIndex, persistentRecords._type, persistentRecords._parentType, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._previousTimeStamp, persistentRecords._previousTimeStepSize, persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._solutionIndex, persistentRecords._solutionAveragesIndex, persistentRecords._solutionCompressedIndex, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolutionIndex, persistentRecords._previousSolutionAveragesIndex, persistentRecords._previousSolutionCompressedIndex, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._updateIndex, persistentRecords._updateAveragesIndex, persistentRecords._updateCompressedIndex, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictorIndex, persistentRecords._extrapolatedPredictorAveragesIndex, persistentRecords._extrapolatedPredictorCompressedIndex, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._extrapolatedPredictorGradientIndex, persistentRecords._extrapolatedPredictorGradient, persistentRecords._fluctuationIndex, persistentRecords._fluctuationAveragesIndex, persistentRecords._fluctuationCompressedIndex, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMinIndex, persistentRecords._solutionMaxIndex, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._facewiseAugmentationStatus, persistentRecords._augmentationStatus, persistentRecords._facewiseCommunicationStatus, persistentRecords._communicationStatus, persistentRecords._facewiseRefinementStatus, persistentRecords._refinementStatus, persistentRecords._previousRefinementStatus, persistentRecords._refinementFlag, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInUpdate, persistentRecords._bytesPerDoFInExtrapolatedPredictor, persistentRecords._bytesPerDoFInFluctuation, persistentRecords._creation) {
   
}


exahype::records::ADERDGCellDescription::ADERDGCellDescription(const int& solverNumber, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed, const int& parentIndex, const Type& type, const Type& parentType, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousTimeStamp, const double& previousTimeStepSize, const double& timeStepSize, const double& timeStamp, const int& solutionIndex, const int& solutionAveragesIndex, const int& solutionCompressedIndex, void* solution, void* solutionAverages, void* solutionCompressed, const int& previousSolutionIndex, const int& previousSolutionAveragesIndex, const int& previousSolutionCompressedIndex, void* previousSolution, void* previousSolutionAverages, void* previousSolutionCompressed, const int& updateIndex, const int& updateAveragesIndex, const int& updateCompressedIndex, void* update, void* updateAverages, void* updateCompressed, const int& extrapolatedPredictorIndex, const int& extrapolatedPredictorAveragesIndex, const int& extrapolatedPredictorCompressedIndex, void* extrapolatedPredictor, void* extrapolatedPredictorAverages, void* extrapolatedPredictorCompressed, const int& extrapolatedPredictorGradientIndex, void* extrapolatedPredictorGradient, const int& fluctuationIndex, const int& fluctuationAveragesIndex, const int& fluctuationCompressedIndex, void* fluctuation, void* fluctuationAverages, void* fluctuationCompressed, const int& solutionMinIndex, const int& solutionMaxIndex, void* solutionMin, void* solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const int& refinementStatus, const int& previousRefinementStatus, const bool& refinementFlag, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation, const Creation& creation):
_persistentRecords(solverNumber, neighbourMergePerformed, parentIndex, type, parentType, level, offset, size, previousTimeStamp, previousTimeStepSize, timeStepSize, timeStamp, solutionIndex, solutionAveragesIndex, solutionCompressedIndex, solution, solutionAverages, solutionCompressed, previousSolutionIndex, previousSolutionAveragesIndex, previousSolutionCompressedIndex, previousSolution, previousSolutionAverages, previousSolutionCompressed, updateIndex, updateAveragesIndex, updateCompressedIndex, update, updateAverages, updateCompressed, extrapolatedPredictorIndex, extrapolatedPredictorAveragesIndex, extrapolatedPredictorCompressedIndex, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, extrapolatedPredictorGradientIndex, extrapolatedPredictorGradient, fluctuationIndex, fluctuationAveragesIndex, fluctuationCompressedIndex, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMinIndex, solutionMaxIndex, solutionMin, solutionMax, facewiseAugmentationStatus, augmentationStatus, facewiseCommunicationStatus, communicationStatus, facewiseRefinementStatus, refinementStatus, previousRefinementStatus, refinementFlag, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation, creation) {
   
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
std::string exahype::records::ADERDGCellDescription::toString(const Creation& param) {
   switch (param) {
      case NotSpecified: return "NotSpecified";
      case UniformRefinement: return "UniformRefinement";
      case AdaptiveRefinement: return "AdaptiveRefinement";
      case AdaptiveCoarsening: return "AdaptiveCoarsening";
      case ReceivedDueToForkOrJoin: return "ReceivedDueToForkOrJoin";
      case ReceivedFromWorker: return "ReceivedFromWorker";
   }
   return "undefined";
}

std::string exahype::records::ADERDGCellDescription::getCreationMapping() {
   return "Creation(NotSpecified=0,UniformRefinement=1,AdaptiveRefinement=2,AdaptiveCoarsening=3,ReceivedDueToForkOrJoin=4,ReceivedFromWorker=5)";
}
std::string exahype::records::ADERDGCellDescription::toString(const Type& param) {
   switch (param) {
      case Leaf: return "Leaf";
      case LeafChecked: return "LeafChecked";
      case LeafInitiatesRefining: return "LeafInitiatesRefining";
      case LeafRefines: return "LeafRefines";
      case LeafProlongates: return "LeafProlongates";
      case Parent: return "Parent";
      case ParentChecked: return "ParentChecked";
      case ParentRequestsCoarseningA: return "ParentRequestsCoarseningA";
      case ParentRequestsCoarseningB: return "ParentRequestsCoarseningB";
      case ParentCoarsens: return "ParentCoarsens";
      case Virtual: return "Virtual";
      case Erased: return "Erased";
   }
   return "undefined";
}

std::string exahype::records::ADERDGCellDescription::getTypeMapping() {
   return "Type(Leaf=0,LeafChecked=1,LeafInitiatesRefining=2,LeafRefines=3,LeafProlongates=4,Parent=5,ParentChecked=6,ParentRequestsCoarseningA=7,ParentRequestsCoarseningB=8,ParentCoarsens=9,Virtual=10,Erased=11)";
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
   out << "type:" << toString(getType());
   out << ",";
   out << "parentType:" << toString(getParentType());
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
   out << "previousTimeStamp:" << getPreviousTimeStamp();
   out << ",";
   out << "previousTimeStepSize:" << getPreviousTimeStepSize();
   out << ",";
   out << "timeStepSize:" << getTimeStepSize();
   out << ",";
   out << "timeStamp:" << getTimeStamp();
   out << ",";
   out << "solutionIndex:" << getSolutionIndex();
   out << ",";
   out << "solutionAveragesIndex:" << getSolutionAveragesIndex();
   out << ",";
   out << "solutionCompressedIndex:" << getSolutionCompressedIndex();
   out << ",";
   out << "solution:" << getSolution();
   out << ",";
   out << "solutionAverages:" << getSolutionAverages();
   out << ",";
   out << "solutionCompressed:" << getSolutionCompressed();
   out << ",";
   out << "previousSolutionIndex:" << getPreviousSolutionIndex();
   out << ",";
   out << "previousSolutionAveragesIndex:" << getPreviousSolutionAveragesIndex();
   out << ",";
   out << "previousSolutionCompressedIndex:" << getPreviousSolutionCompressedIndex();
   out << ",";
   out << "previousSolution:" << getPreviousSolution();
   out << ",";
   out << "previousSolutionAverages:" << getPreviousSolutionAverages();
   out << ",";
   out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
   out << ",";
   out << "updateIndex:" << getUpdateIndex();
   out << ",";
   out << "updateAveragesIndex:" << getUpdateAveragesIndex();
   out << ",";
   out << "updateCompressedIndex:" << getUpdateCompressedIndex();
   out << ",";
   out << "update:" << getUpdate();
   out << ",";
   out << "updateAverages:" << getUpdateAverages();
   out << ",";
   out << "updateCompressed:" << getUpdateCompressed();
   out << ",";
   out << "extrapolatedPredictorIndex:" << getExtrapolatedPredictorIndex();
   out << ",";
   out << "extrapolatedPredictorAveragesIndex:" << getExtrapolatedPredictorAveragesIndex();
   out << ",";
   out << "extrapolatedPredictorCompressedIndex:" << getExtrapolatedPredictorCompressedIndex();
   out << ",";
   out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
   out << ",";
   out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
   out << ",";
   out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
   out << ",";
   out << "extrapolatedPredictorGradientIndex:" << getExtrapolatedPredictorGradientIndex();
   out << ",";
   out << "extrapolatedPredictorGradient:" << getExtrapolatedPredictorGradient();
   out << ",";
   out << "fluctuationIndex:" << getFluctuationIndex();
   out << ",";
   out << "fluctuationAveragesIndex:" << getFluctuationAveragesIndex();
   out << ",";
   out << "fluctuationCompressedIndex:" << getFluctuationCompressedIndex();
   out << ",";
   out << "fluctuation:" << getFluctuation();
   out << ",";
   out << "fluctuationAverages:" << getFluctuationAverages();
   out << ",";
   out << "fluctuationCompressed:" << getFluctuationCompressed();
   out << ",";
   out << "solutionMinIndex:" << getSolutionMinIndex();
   out << ",";
   out << "solutionMaxIndex:" << getSolutionMaxIndex();
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
   out << "refinementStatus:" << getRefinementStatus();
   out << ",";
   out << "previousRefinementStatus:" << getPreviousRefinementStatus();
   out << ",";
   out << "refinementFlag:" << getRefinementFlag();
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
   out << ",";
   out << "creation:" << toString(getCreation());
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
      getType(),
      getParentType(),
      getLevel(),
      getOffset(),
      getSize(),
      getPreviousTimeStamp(),
      getPreviousTimeStepSize(),
      getTimeStepSize(),
      getTimeStamp(),
      getSolutionIndex(),
      getSolutionAveragesIndex(),
      getSolutionCompressedIndex(),
      getSolution(),
      getSolutionAverages(),
      getSolutionCompressed(),
      getPreviousSolutionIndex(),
      getPreviousSolutionAveragesIndex(),
      getPreviousSolutionCompressedIndex(),
      getPreviousSolution(),
      getPreviousSolutionAverages(),
      getPreviousSolutionCompressed(),
      getUpdateIndex(),
      getUpdateAveragesIndex(),
      getUpdateCompressedIndex(),
      getUpdate(),
      getUpdateAverages(),
      getUpdateCompressed(),
      getExtrapolatedPredictorIndex(),
      getExtrapolatedPredictorAveragesIndex(),
      getExtrapolatedPredictorCompressedIndex(),
      getExtrapolatedPredictor(),
      getExtrapolatedPredictorAverages(),
      getExtrapolatedPredictorCompressed(),
      getExtrapolatedPredictorGradientIndex(),
      getExtrapolatedPredictorGradient(),
      getFluctuationIndex(),
      getFluctuationAveragesIndex(),
      getFluctuationCompressedIndex(),
      getFluctuation(),
      getFluctuationAverages(),
      getFluctuationCompressed(),
      getSolutionMinIndex(),
      getSolutionMaxIndex(),
      getSolutionMin(),
      getSolutionMax(),
      getFacewiseAugmentationStatus(),
      getAugmentationStatus(),
      getFacewiseCommunicationStatus(),
      getCommunicationStatus(),
      getFacewiseRefinementStatus(),
      getRefinementStatus(),
      getPreviousRefinementStatus(),
      getRefinementFlag(),
      getCompressionState(),
      getBytesPerDoFInPreviousSolution(),
      getBytesPerDoFInSolution(),
      getBytesPerDoFInUpdate(),
      getBytesPerDoFInExtrapolatedPredictor(),
      getBytesPerDoFInFluctuation(),
      getCreation()
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
         const int Attributes = 21;
         #else
         const int Attributes = 22;
         #endif
         MPI_Datatype subtypes[Attributes] = {
              MPI_INT		 //solverNumber
            , MPI_INT		 //type
            , MPI_INT		 //parentType
            , MPI_INT		 //level
            , MPI_DOUBLE		 //offset
            , MPI_DOUBLE		 //size
            , MPI_DOUBLE		 //previousTimeStamp
            , MPI_DOUBLE		 //previousTimeStepSize
            , MPI_DOUBLE		 //timeStepSize
            , MPI_DOUBLE		 //timeStamp
            , MPI_INT		 //augmentationStatus
            , MPI_INT		 //communicationStatus
            , MPI_INT		 //refinementStatus
            , MPI_INT		 //previousRefinementStatus
            , MPI_CXX_BOOL		 //refinementFlag
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
            , 1		 //type
            , 1		 //parentType
            , 1		 //level
            , DIMENSIONS		 //offset
            , DIMENSIONS		 //size
            , 1		 //previousTimeStamp
            , 1		 //previousTimeStepSize
            , 1		 //timeStepSize
            , 1		 //timeStamp
            , 1		 //augmentationStatus
            , 1		 //communicationStatus
            , 1		 //refinementStatus
            , 1		 //previousRefinementStatus
            , 1		 //refinementFlag
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
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[1] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[1] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentType))), 		&disp[2] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentType))), 		&disp[2] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[3] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[3] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[4] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[4] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[5] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[5] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[7] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[7] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[8] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[8] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStamp))), 		&disp[9] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStamp))), 		&disp[9] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[10] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[10] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[11] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[11] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[12] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[12] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[13] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[13] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[14] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[14] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[15] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[15] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[16] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[16] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[17] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[17] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[18] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[18] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[19] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[19] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[20] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[20] );
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
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]))), 		&disp[21] );
         disp[21] -= base;
         disp[21] += disp[0];
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
         const int Attributes = 63;
         #else
         const int Attributes = 64;
         #endif
         MPI_Datatype subtypes[Attributes] = {
              MPI_INT		 //solverNumber
            , MPI_BYTE		 //neighbourMergePerformed
            , MPI_INT		 //parentIndex
            , MPI_INT		 //type
            , MPI_INT		 //parentType
            , MPI_INT		 //level
            , MPI_DOUBLE		 //offset
            , MPI_DOUBLE		 //size
            , MPI_DOUBLE		 //previousTimeStamp
            , MPI_DOUBLE		 //previousTimeStepSize
            , MPI_DOUBLE		 //timeStepSize
            , MPI_DOUBLE		 //timeStamp
            , MPI_INT		 //solutionIndex
            , MPI_INT		 //solutionAveragesIndex
            , MPI_INT		 //solutionCompressedIndex
            , MPI_UNSIGNED_LONG		 //solution
            , MPI_UNSIGNED_LONG		 //solutionAverages
            , MPI_UNSIGNED_LONG		 //solutionCompressed
            , MPI_INT		 //previousSolutionIndex
            , MPI_INT		 //previousSolutionAveragesIndex
            , MPI_INT		 //previousSolutionCompressedIndex
            , MPI_UNSIGNED_LONG		 //previousSolution
            , MPI_UNSIGNED_LONG		 //previousSolutionAverages
            , MPI_UNSIGNED_LONG		 //previousSolutionCompressed
            , MPI_INT		 //updateIndex
            , MPI_INT		 //updateAveragesIndex
            , MPI_INT		 //updateCompressedIndex
            , MPI_UNSIGNED_LONG		 //update
            , MPI_UNSIGNED_LONG		 //updateAverages
            , MPI_UNSIGNED_LONG		 //updateCompressed
            , MPI_INT		 //extrapolatedPredictorIndex
            , MPI_INT		 //extrapolatedPredictorAveragesIndex
            , MPI_INT		 //extrapolatedPredictorCompressedIndex
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictor
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictorAverages
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictorCompressed
            , MPI_INT		 //extrapolatedPredictorGradientIndex
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictorGradient
            , MPI_INT		 //fluctuationIndex
            , MPI_INT		 //fluctuationAveragesIndex
            , MPI_INT		 //fluctuationCompressedIndex
            , MPI_UNSIGNED_LONG		 //fluctuation
            , MPI_UNSIGNED_LONG		 //fluctuationAverages
            , MPI_UNSIGNED_LONG		 //fluctuationCompressed
            , MPI_INT		 //solutionMinIndex
            , MPI_INT		 //solutionMaxIndex
            , MPI_UNSIGNED_LONG		 //solutionMin
            , MPI_UNSIGNED_LONG		 //solutionMax
            , MPI_INT		 //facewiseAugmentationStatus
            , MPI_INT		 //augmentationStatus
            , MPI_INT		 //facewiseCommunicationStatus
            , MPI_INT		 //communicationStatus
            , MPI_INT		 //facewiseRefinementStatus
            , MPI_INT		 //refinementStatus
            , MPI_INT		 //previousRefinementStatus
            , MPI_CXX_BOOL		 //refinementFlag
            , MPI_INT		 //compressionState
            , MPI_INT		 //bytesPerDoFInPreviousSolution
            , MPI_INT		 //bytesPerDoFInSolution
            , MPI_INT		 //bytesPerDoFInUpdate
            , MPI_INT		 //bytesPerDoFInExtrapolatedPredictor
            , MPI_INT		 //bytesPerDoFInFluctuation
            , MPI_INT		 //creation
            #ifndef MPI2
            , MPI_UB
            #endif
            
         };
         
         int blocklen[Attributes] = {
              1		 //solverNumber
            , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
            , 1		 //parentIndex
            , 1		 //type
            , 1		 //parentType
            , 1		 //level
            , DIMENSIONS		 //offset
            , DIMENSIONS		 //size
            , 1		 //previousTimeStamp
            , 1		 //previousTimeStepSize
            , 1		 //timeStepSize
            , 1		 //timeStamp
            , 1		 //solutionIndex
            , 1		 //solutionAveragesIndex
            , 1		 //solutionCompressedIndex
            , 1		 //solution
            , 1		 //solutionAverages
            , 1		 //solutionCompressed
            , 1		 //previousSolutionIndex
            , 1		 //previousSolutionAveragesIndex
            , 1		 //previousSolutionCompressedIndex
            , 1		 //previousSolution
            , 1		 //previousSolutionAverages
            , 1		 //previousSolutionCompressed
            , 1		 //updateIndex
            , 1		 //updateAveragesIndex
            , 1		 //updateCompressedIndex
            , 1		 //update
            , 1		 //updateAverages
            , 1		 //updateCompressed
            , 1		 //extrapolatedPredictorIndex
            , 1		 //extrapolatedPredictorAveragesIndex
            , 1		 //extrapolatedPredictorCompressedIndex
            , 1		 //extrapolatedPredictor
            , 1		 //extrapolatedPredictorAverages
            , 1		 //extrapolatedPredictorCompressed
            , 1		 //extrapolatedPredictorGradientIndex
            , 1		 //extrapolatedPredictorGradient
            , 1		 //fluctuationIndex
            , 1		 //fluctuationAveragesIndex
            , 1		 //fluctuationCompressedIndex
            , 1		 //fluctuation
            , 1		 //fluctuationAverages
            , 1		 //fluctuationCompressed
            , 1		 //solutionMinIndex
            , 1		 //solutionMaxIndex
            , 1		 //solutionMin
            , 1		 //solutionMax
            , DIMENSIONS_TIMES_TWO		 //facewiseAugmentationStatus
            , 1		 //augmentationStatus
            , DIMENSIONS_TIMES_TWO		 //facewiseCommunicationStatus
            , 1		 //communicationStatus
            , DIMENSIONS_TIMES_TWO		 //facewiseRefinementStatus
            , 1		 //refinementStatus
            , 1		 //previousRefinementStatus
            , 1		 //refinementFlag
            , 1		 //compressionState
            , 1		 //bytesPerDoFInPreviousSolution
            , 1		 //bytesPerDoFInSolution
            , 1		 //bytesPerDoFInUpdate
            , 1		 //bytesPerDoFInExtrapolatedPredictor
            , 1		 //bytesPerDoFInFluctuation
            , 1		 //creation
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
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._neighbourMergePerformed[0]))), 		&disp[1] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._neighbourMergePerformed[0]))), 		&disp[1] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[2] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[2] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[3] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[3] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentType))), 		&disp[4] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentType))), 		&disp[4] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[5] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[5] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[6] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[6] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[7] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[7] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[8] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[8] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[9] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[9] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[10] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[10] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStamp))), 		&disp[11] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._timeStamp))), 		&disp[11] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionIndex))), 		&disp[12] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionIndex))), 		&disp[12] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAveragesIndex))), 		&disp[13] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAveragesIndex))), 		&disp[13] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressedIndex))), 		&disp[14] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressedIndex))), 		&disp[14] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[15] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[15] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[16] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[16] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[17] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[17] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionIndex))), 		&disp[18] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionIndex))), 		&disp[18] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAveragesIndex))), 		&disp[19] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAveragesIndex))), 		&disp[19] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressedIndex))), 		&disp[20] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressedIndex))), 		&disp[20] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[21] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[21] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[22] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[22] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[23] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[23] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateIndex))), 		&disp[24] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateIndex))), 		&disp[24] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAveragesIndex))), 		&disp[25] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAveragesIndex))), 		&disp[25] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressedIndex))), 		&disp[26] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressedIndex))), 		&disp[26] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[27] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[27] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[28] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[28] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[29] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[29] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorIndex))), 		&disp[30] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorIndex))), 		&disp[30] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAveragesIndex))), 		&disp[31] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAveragesIndex))), 		&disp[31] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressedIndex))), 		&disp[32] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressedIndex))), 		&disp[32] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[33] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[33] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[34] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[34] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[35] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[35] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorGradientIndex))), 		&disp[36] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorGradientIndex))), 		&disp[36] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorGradient))), 		&disp[37] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorGradient))), 		&disp[37] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationIndex))), 		&disp[38] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationIndex))), 		&disp[38] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAveragesIndex))), 		&disp[39] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAveragesIndex))), 		&disp[39] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressedIndex))), 		&disp[40] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressedIndex))), 		&disp[40] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[41] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[41] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[42] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[42] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[43] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[43] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMinIndex))), 		&disp[44] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMinIndex))), 		&disp[44] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMaxIndex))), 		&disp[45] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMaxIndex))), 		&disp[45] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[46] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[46] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[47] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[47] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[48] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[48] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[49] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._augmentationStatus))), 		&disp[49] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[50] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[50] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[51] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._communicationStatus))), 		&disp[51] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[52] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[52] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[53] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementStatus))), 		&disp[53] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[54] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousRefinementStatus))), 		&disp[54] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[55] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementFlag))), 		&disp[55] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[56] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[56] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[57] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[57] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[58] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[58] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[59] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[59] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[60] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[60] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[61] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[61] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._creation))), 		&disp[62] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._creation))), 		&disp[62] );
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
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]))), 		&disp[63] );
         disp[63] -= base;
         disp[63] += disp[0];
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
   if ((25 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((25 < (8 * sizeof(int))));
   
}


exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed, const int& parentIndex, const Type& type, const Type& parentType, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousTimeStamp, const double& previousTimeStepSize, const double& timeStepSize, const double& timeStamp, const int& solutionIndex, const int& solutionAveragesIndex, const int& solutionCompressedIndex, void* solution, void* solutionAverages, void* solutionCompressed, const int& previousSolutionIndex, const int& previousSolutionAveragesIndex, const int& previousSolutionCompressedIndex, void* previousSolution, void* previousSolutionAverages, void* previousSolutionCompressed, const int& updateIndex, const int& updateAveragesIndex, const int& updateCompressedIndex, void* update, void* updateAverages, void* updateCompressed, const int& extrapolatedPredictorIndex, const int& extrapolatedPredictorAveragesIndex, const int& extrapolatedPredictorCompressedIndex, void* extrapolatedPredictor, void* extrapolatedPredictorAverages, void* extrapolatedPredictorCompressed, const int& extrapolatedPredictorGradientIndex, void* extrapolatedPredictorGradient, const int& fluctuationIndex, const int& fluctuationAveragesIndex, const int& fluctuationCompressedIndex, void* fluctuation, void* fluctuationAverages, void* fluctuationCompressed, const int& solutionMinIndex, const int& solutionMaxIndex, void* solutionMin, void* solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const int& refinementStatus, const int& previousRefinementStatus, const bool& refinementFlag, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation, const Creation& creation):
_solverNumber(solverNumber),
_parentIndex(parentIndex),
_level(level),
_offset(offset),
_size(size),
_previousTimeStamp(previousTimeStamp),
_previousTimeStepSize(previousTimeStepSize),
_timeStepSize(timeStepSize),
_timeStamp(timeStamp),
_solutionIndex(solutionIndex),
_solutionAveragesIndex(solutionAveragesIndex),
_solutionCompressedIndex(solutionCompressedIndex),
_solution(solution),
_solutionAverages(solutionAverages),
_solutionCompressed(solutionCompressed),
_previousSolutionIndex(previousSolutionIndex),
_previousSolutionAveragesIndex(previousSolutionAveragesIndex),
_previousSolutionCompressedIndex(previousSolutionCompressedIndex),
_previousSolution(previousSolution),
_previousSolutionAverages(previousSolutionAverages),
_previousSolutionCompressed(previousSolutionCompressed),
_updateIndex(updateIndex),
_updateAveragesIndex(updateAveragesIndex),
_updateCompressedIndex(updateCompressedIndex),
_update(update),
_updateAverages(updateAverages),
_updateCompressed(updateCompressed),
_extrapolatedPredictorIndex(extrapolatedPredictorIndex),
_extrapolatedPredictorAveragesIndex(extrapolatedPredictorAveragesIndex),
_extrapolatedPredictorCompressedIndex(extrapolatedPredictorCompressedIndex),
_extrapolatedPredictor(extrapolatedPredictor),
_extrapolatedPredictorAverages(extrapolatedPredictorAverages),
_extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
_extrapolatedPredictorGradientIndex(extrapolatedPredictorGradientIndex),
_extrapolatedPredictorGradient(extrapolatedPredictorGradient),
_fluctuationIndex(fluctuationIndex),
_fluctuationAveragesIndex(fluctuationAveragesIndex),
_fluctuationCompressedIndex(fluctuationCompressedIndex),
_fluctuation(fluctuation),
_fluctuationAverages(fluctuationAverages),
_fluctuationCompressed(fluctuationCompressed),
_solutionMinIndex(solutionMinIndex),
_solutionMaxIndex(solutionMaxIndex),
_solutionMin(solutionMin),
_solutionMax(solutionMax),
_facewiseAugmentationStatus(facewiseAugmentationStatus),
_augmentationStatus(augmentationStatus),
_facewiseCommunicationStatus(facewiseCommunicationStatus),
_communicationStatus(communicationStatus),
_facewiseRefinementStatus(facewiseRefinementStatus),
_refinementStatus(refinementStatus),
_previousRefinementStatus(previousRefinementStatus),
_refinementFlag(refinementFlag),
_creation(creation) {
   setNeighbourMergePerformed(neighbourMergePerformed);
   setType(type);
   setParentType(parentType);
   setCompressionState(compressionState);
   setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
   setBytesPerDoFInSolution(bytesPerDoFInSolution);
   setBytesPerDoFInUpdate(bytesPerDoFInUpdate);
   setBytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor);
   setBytesPerDoFInFluctuation(bytesPerDoFInFluctuation);
   if ((25 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((25 < (8 * sizeof(int))));
   
}

exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked() {
   if ((25 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((25 < (8 * sizeof(int))));
   
}


exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._solverNumber, persistentRecords.getNeighbourMergePerformed(), persistentRecords._parentIndex, persistentRecords.getType(), persistentRecords.getParentType(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._previousTimeStamp, persistentRecords._previousTimeStepSize, persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._solutionIndex, persistentRecords._solutionAveragesIndex, persistentRecords._solutionCompressedIndex, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolutionIndex, persistentRecords._previousSolutionAveragesIndex, persistentRecords._previousSolutionCompressedIndex, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._updateIndex, persistentRecords._updateAveragesIndex, persistentRecords._updateCompressedIndex, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictorIndex, persistentRecords._extrapolatedPredictorAveragesIndex, persistentRecords._extrapolatedPredictorCompressedIndex, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._extrapolatedPredictorGradientIndex, persistentRecords._extrapolatedPredictorGradient, persistentRecords._fluctuationIndex, persistentRecords._fluctuationAveragesIndex, persistentRecords._fluctuationCompressedIndex, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMinIndex, persistentRecords._solutionMaxIndex, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._facewiseAugmentationStatus, persistentRecords._augmentationStatus, persistentRecords._facewiseCommunicationStatus, persistentRecords._communicationStatus, persistentRecords._facewiseRefinementStatus, persistentRecords._refinementStatus, persistentRecords._previousRefinementStatus, persistentRecords._refinementFlag, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInUpdate(), persistentRecords.getBytesPerDoFInExtrapolatedPredictor(), persistentRecords.getBytesPerDoFInFluctuation(), persistentRecords._creation) {
   if ((25 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((25 < (8 * sizeof(int))));
   
}


exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const int& solverNumber, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,signed char>& neighbourMergePerformed, const int& parentIndex, const Type& type, const Type& parentType, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& previousTimeStamp, const double& previousTimeStepSize, const double& timeStepSize, const double& timeStamp, const int& solutionIndex, const int& solutionAveragesIndex, const int& solutionCompressedIndex, void* solution, void* solutionAverages, void* solutionCompressed, const int& previousSolutionIndex, const int& previousSolutionAveragesIndex, const int& previousSolutionCompressedIndex, void* previousSolution, void* previousSolutionAverages, void* previousSolutionCompressed, const int& updateIndex, const int& updateAveragesIndex, const int& updateCompressedIndex, void* update, void* updateAverages, void* updateCompressed, const int& extrapolatedPredictorIndex, const int& extrapolatedPredictorAveragesIndex, const int& extrapolatedPredictorCompressedIndex, void* extrapolatedPredictor, void* extrapolatedPredictorAverages, void* extrapolatedPredictorCompressed, const int& extrapolatedPredictorGradientIndex, void* extrapolatedPredictorGradient, const int& fluctuationIndex, const int& fluctuationAveragesIndex, const int& fluctuationCompressedIndex, void* fluctuation, void* fluctuationAverages, void* fluctuationCompressed, const int& solutionMinIndex, const int& solutionMaxIndex, void* solutionMin, void* solutionMax, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseAugmentationStatus, const int& augmentationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseCommunicationStatus, const int& communicationStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& facewiseRefinementStatus, const int& refinementStatus, const int& previousRefinementStatus, const bool& refinementFlag, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation, const Creation& creation):
_persistentRecords(solverNumber, neighbourMergePerformed, parentIndex, type, parentType, level, offset, size, previousTimeStamp, previousTimeStepSize, timeStepSize, timeStamp, solutionIndex, solutionAveragesIndex, solutionCompressedIndex, solution, solutionAverages, solutionCompressed, previousSolutionIndex, previousSolutionAveragesIndex, previousSolutionCompressedIndex, previousSolution, previousSolutionAverages, previousSolutionCompressed, updateIndex, updateAveragesIndex, updateCompressedIndex, update, updateAverages, updateCompressed, extrapolatedPredictorIndex, extrapolatedPredictorAveragesIndex, extrapolatedPredictorCompressedIndex, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, extrapolatedPredictorGradientIndex, extrapolatedPredictorGradient, fluctuationIndex, fluctuationAveragesIndex, fluctuationCompressedIndex, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMinIndex, solutionMaxIndex, solutionMin, solutionMax, facewiseAugmentationStatus, augmentationStatus, facewiseCommunicationStatus, communicationStatus, facewiseRefinementStatus, refinementStatus, previousRefinementStatus, refinementFlag, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation, creation) {
   if ((25 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((25 < (8 * sizeof(int))));
   
}


exahype::records::ADERDGCellDescriptionPacked::~ADERDGCellDescriptionPacked() { }

std::string exahype::records::ADERDGCellDescriptionPacked::toString(const Type& param) {
   return exahype::records::ADERDGCellDescription::toString(param);
}

std::string exahype::records::ADERDGCellDescriptionPacked::getTypeMapping() {
   return exahype::records::ADERDGCellDescription::getTypeMapping();
}

std::string exahype::records::ADERDGCellDescriptionPacked::toString(const CompressionState& param) {
   return exahype::records::ADERDGCellDescription::toString(param);
}

std::string exahype::records::ADERDGCellDescriptionPacked::getCompressionStateMapping() {
   return exahype::records::ADERDGCellDescription::getCompressionStateMapping();
}

std::string exahype::records::ADERDGCellDescriptionPacked::toString(const Creation& param) {
   return exahype::records::ADERDGCellDescription::toString(param);
}

std::string exahype::records::ADERDGCellDescriptionPacked::getCreationMapping() {
   return exahype::records::ADERDGCellDescription::getCreationMapping();
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
   out << "type:" << toString(getType());
   out << ",";
   out << "parentType:" << toString(getParentType());
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
   out << "previousTimeStamp:" << getPreviousTimeStamp();
   out << ",";
   out << "previousTimeStepSize:" << getPreviousTimeStepSize();
   out << ",";
   out << "timeStepSize:" << getTimeStepSize();
   out << ",";
   out << "timeStamp:" << getTimeStamp();
   out << ",";
   out << "solutionIndex:" << getSolutionIndex();
   out << ",";
   out << "solutionAveragesIndex:" << getSolutionAveragesIndex();
   out << ",";
   out << "solutionCompressedIndex:" << getSolutionCompressedIndex();
   out << ",";
   out << "solution:" << getSolution();
   out << ",";
   out << "solutionAverages:" << getSolutionAverages();
   out << ",";
   out << "solutionCompressed:" << getSolutionCompressed();
   out << ",";
   out << "previousSolutionIndex:" << getPreviousSolutionIndex();
   out << ",";
   out << "previousSolutionAveragesIndex:" << getPreviousSolutionAveragesIndex();
   out << ",";
   out << "previousSolutionCompressedIndex:" << getPreviousSolutionCompressedIndex();
   out << ",";
   out << "previousSolution:" << getPreviousSolution();
   out << ",";
   out << "previousSolutionAverages:" << getPreviousSolutionAverages();
   out << ",";
   out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
   out << ",";
   out << "updateIndex:" << getUpdateIndex();
   out << ",";
   out << "updateAveragesIndex:" << getUpdateAveragesIndex();
   out << ",";
   out << "updateCompressedIndex:" << getUpdateCompressedIndex();
   out << ",";
   out << "update:" << getUpdate();
   out << ",";
   out << "updateAverages:" << getUpdateAverages();
   out << ",";
   out << "updateCompressed:" << getUpdateCompressed();
   out << ",";
   out << "extrapolatedPredictorIndex:" << getExtrapolatedPredictorIndex();
   out << ",";
   out << "extrapolatedPredictorAveragesIndex:" << getExtrapolatedPredictorAveragesIndex();
   out << ",";
   out << "extrapolatedPredictorCompressedIndex:" << getExtrapolatedPredictorCompressedIndex();
   out << ",";
   out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
   out << ",";
   out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
   out << ",";
   out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
   out << ",";
   out << "extrapolatedPredictorGradientIndex:" << getExtrapolatedPredictorGradientIndex();
   out << ",";
   out << "extrapolatedPredictorGradient:" << getExtrapolatedPredictorGradient();
   out << ",";
   out << "fluctuationIndex:" << getFluctuationIndex();
   out << ",";
   out << "fluctuationAveragesIndex:" << getFluctuationAveragesIndex();
   out << ",";
   out << "fluctuationCompressedIndex:" << getFluctuationCompressedIndex();
   out << ",";
   out << "fluctuation:" << getFluctuation();
   out << ",";
   out << "fluctuationAverages:" << getFluctuationAverages();
   out << ",";
   out << "fluctuationCompressed:" << getFluctuationCompressed();
   out << ",";
   out << "solutionMinIndex:" << getSolutionMinIndex();
   out << ",";
   out << "solutionMaxIndex:" << getSolutionMaxIndex();
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
   out << "refinementStatus:" << getRefinementStatus();
   out << ",";
   out << "previousRefinementStatus:" << getPreviousRefinementStatus();
   out << ",";
   out << "refinementFlag:" << getRefinementFlag();
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
   out << ",";
   out << "creation:" << toString(getCreation());
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
      getType(),
      getParentType(),
      getLevel(),
      getOffset(),
      getSize(),
      getPreviousTimeStamp(),
      getPreviousTimeStepSize(),
      getTimeStepSize(),
      getTimeStamp(),
      getSolutionIndex(),
      getSolutionAveragesIndex(),
      getSolutionCompressedIndex(),
      getSolution(),
      getSolutionAverages(),
      getSolutionCompressed(),
      getPreviousSolutionIndex(),
      getPreviousSolutionAveragesIndex(),
      getPreviousSolutionCompressedIndex(),
      getPreviousSolution(),
      getPreviousSolutionAverages(),
      getPreviousSolutionCompressed(),
      getUpdateIndex(),
      getUpdateAveragesIndex(),
      getUpdateCompressedIndex(),
      getUpdate(),
      getUpdateAverages(),
      getUpdateCompressed(),
      getExtrapolatedPredictorIndex(),
      getExtrapolatedPredictorAveragesIndex(),
      getExtrapolatedPredictorCompressedIndex(),
      getExtrapolatedPredictor(),
      getExtrapolatedPredictorAverages(),
      getExtrapolatedPredictorCompressed(),
      getExtrapolatedPredictorGradientIndex(),
      getExtrapolatedPredictorGradient(),
      getFluctuationIndex(),
      getFluctuationAveragesIndex(),
      getFluctuationCompressedIndex(),
      getFluctuation(),
      getFluctuationAverages(),
      getFluctuationCompressed(),
      getSolutionMinIndex(),
      getSolutionMaxIndex(),
      getSolutionMin(),
      getSolutionMax(),
      getFacewiseAugmentationStatus(),
      getAugmentationStatus(),
      getFacewiseCommunicationStatus(),
      getCommunicationStatus(),
      getFacewiseRefinementStatus(),
      getRefinementStatus(),
      getPreviousRefinementStatus(),
      getRefinementFlag(),
      getCompressionState(),
      getBytesPerDoFInPreviousSolution(),
      getBytesPerDoFInSolution(),
      getBytesPerDoFInUpdate(),
      getBytesPerDoFInExtrapolatedPredictor(),
      getBytesPerDoFInFluctuation(),
      getCreation()
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
         const int Attributes = 14;
         #else
         const int Attributes = 15;
         #endif
         MPI_Datatype subtypes[Attributes] = {
              MPI_INT		 //solverNumber
            , MPI_INT		 //level
            , MPI_DOUBLE		 //offset
            , MPI_DOUBLE		 //size
            , MPI_DOUBLE		 //previousTimeStamp
            , MPI_DOUBLE		 //previousTimeStepSize
            , MPI_DOUBLE		 //timeStepSize
            , MPI_DOUBLE		 //timeStamp
            , MPI_INT		 //augmentationStatus
            , MPI_INT		 //communicationStatus
            , MPI_INT		 //refinementStatus
            , MPI_INT		 //previousRefinementStatus
            , MPI_CXX_BOOL		 //refinementFlag
            , MPI_INT		 //_packedRecords0
            #ifndef MPI2
            , MPI_UB
            #endif
            
         };
         
         int blocklen[Attributes] = {
              1		 //solverNumber
            , 1		 //level
            , DIMENSIONS		 //offset
            , DIMENSIONS		 //size
            , 1		 //previousTimeStamp
            , 1		 //previousTimeStepSize
            , 1		 //timeStepSize
            , 1		 //timeStamp
            , 1		 //augmentationStatus
            , 1		 //communicationStatus
            , 1		 //refinementStatus
            , 1		 //previousRefinementStatus
            , 1		 //refinementFlag
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
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[1] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[1] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[2] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[2] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[3] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[3] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[4] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[4] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[5] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[5] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[6] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[6] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[7] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[7] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[8] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[8] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[9] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[9] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[10] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[10] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[11] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[11] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[12] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[12] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[13] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[13] );
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
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]))), 		&disp[14] );
         disp[14] -= base;
         disp[14] += disp[0];
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
         const int Attributes = 56;
         #else
         const int Attributes = 57;
         #endif
         MPI_Datatype subtypes[Attributes] = {
              MPI_INT		 //solverNumber
            , MPI_BYTE		 //neighbourMergePerformed
            , MPI_INT		 //parentIndex
            , MPI_INT		 //level
            , MPI_DOUBLE		 //offset
            , MPI_DOUBLE		 //size
            , MPI_DOUBLE		 //previousTimeStamp
            , MPI_DOUBLE		 //previousTimeStepSize
            , MPI_DOUBLE		 //timeStepSize
            , MPI_DOUBLE		 //timeStamp
            , MPI_INT		 //solutionIndex
            , MPI_INT		 //solutionAveragesIndex
            , MPI_INT		 //solutionCompressedIndex
            , MPI_UNSIGNED_LONG		 //solution
            , MPI_UNSIGNED_LONG		 //solutionAverages
            , MPI_UNSIGNED_LONG		 //solutionCompressed
            , MPI_INT		 //previousSolutionIndex
            , MPI_INT		 //previousSolutionAveragesIndex
            , MPI_INT		 //previousSolutionCompressedIndex
            , MPI_UNSIGNED_LONG		 //previousSolution
            , MPI_UNSIGNED_LONG		 //previousSolutionAverages
            , MPI_UNSIGNED_LONG		 //previousSolutionCompressed
            , MPI_INT		 //updateIndex
            , MPI_INT		 //updateAveragesIndex
            , MPI_INT		 //updateCompressedIndex
            , MPI_UNSIGNED_LONG		 //update
            , MPI_UNSIGNED_LONG		 //updateAverages
            , MPI_UNSIGNED_LONG		 //updateCompressed
            , MPI_INT		 //extrapolatedPredictorIndex
            , MPI_INT		 //extrapolatedPredictorAveragesIndex
            , MPI_INT		 //extrapolatedPredictorCompressedIndex
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictor
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictorAverages
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictorCompressed
            , MPI_INT		 //extrapolatedPredictorGradientIndex
            , MPI_UNSIGNED_LONG		 //extrapolatedPredictorGradient
            , MPI_INT		 //fluctuationIndex
            , MPI_INT		 //fluctuationAveragesIndex
            , MPI_INT		 //fluctuationCompressedIndex
            , MPI_UNSIGNED_LONG		 //fluctuation
            , MPI_UNSIGNED_LONG		 //fluctuationAverages
            , MPI_UNSIGNED_LONG		 //fluctuationCompressed
            , MPI_INT		 //solutionMinIndex
            , MPI_INT		 //solutionMaxIndex
            , MPI_UNSIGNED_LONG		 //solutionMin
            , MPI_UNSIGNED_LONG		 //solutionMax
            , MPI_INT		 //facewiseAugmentationStatus
            , MPI_INT		 //augmentationStatus
            , MPI_INT		 //facewiseCommunicationStatus
            , MPI_INT		 //communicationStatus
            , MPI_INT		 //facewiseRefinementStatus
            , MPI_INT		 //refinementStatus
            , MPI_INT		 //previousRefinementStatus
            , MPI_CXX_BOOL		 //refinementFlag
            , MPI_INT		 //creation
            , MPI_INT		 //_packedRecords0
            #ifndef MPI2
            , MPI_UB
            #endif
            
         };
         
         int blocklen[Attributes] = {
              1		 //solverNumber
            , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
            , 1		 //parentIndex
            , 1		 //level
            , DIMENSIONS		 //offset
            , DIMENSIONS		 //size
            , 1		 //previousTimeStamp
            , 1		 //previousTimeStepSize
            , 1		 //timeStepSize
            , 1		 //timeStamp
            , 1		 //solutionIndex
            , 1		 //solutionAveragesIndex
            , 1		 //solutionCompressedIndex
            , 1		 //solution
            , 1		 //solutionAverages
            , 1		 //solutionCompressed
            , 1		 //previousSolutionIndex
            , 1		 //previousSolutionAveragesIndex
            , 1		 //previousSolutionCompressedIndex
            , 1		 //previousSolution
            , 1		 //previousSolutionAverages
            , 1		 //previousSolutionCompressed
            , 1		 //updateIndex
            , 1		 //updateAveragesIndex
            , 1		 //updateCompressedIndex
            , 1		 //update
            , 1		 //updateAverages
            , 1		 //updateCompressed
            , 1		 //extrapolatedPredictorIndex
            , 1		 //extrapolatedPredictorAveragesIndex
            , 1		 //extrapolatedPredictorCompressedIndex
            , 1		 //extrapolatedPredictor
            , 1		 //extrapolatedPredictorAverages
            , 1		 //extrapolatedPredictorCompressed
            , 1		 //extrapolatedPredictorGradientIndex
            , 1		 //extrapolatedPredictorGradient
            , 1		 //fluctuationIndex
            , 1		 //fluctuationAveragesIndex
            , 1		 //fluctuationCompressedIndex
            , 1		 //fluctuation
            , 1		 //fluctuationAverages
            , 1		 //fluctuationCompressed
            , 1		 //solutionMinIndex
            , 1		 //solutionMaxIndex
            , 1		 //solutionMin
            , 1		 //solutionMax
            , DIMENSIONS_TIMES_TWO		 //facewiseAugmentationStatus
            , 1		 //augmentationStatus
            , DIMENSIONS_TIMES_TWO		 //facewiseCommunicationStatus
            , 1		 //communicationStatus
            , DIMENSIONS_TIMES_TWO		 //facewiseRefinementStatus
            , 1		 //refinementStatus
            , 1		 //previousRefinementStatus
            , 1		 //refinementFlag
            , 1		 //creation
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
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed[0]))), 		&disp[1] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed[0]))), 		&disp[1] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[2] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[2] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[3] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[3] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[4] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[4] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[5] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[5] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[7] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[7] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[8] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[8] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[9] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[9] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionIndex))), 		&disp[10] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionIndex))), 		&disp[10] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAveragesIndex))), 		&disp[11] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAveragesIndex))), 		&disp[11] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressedIndex))), 		&disp[12] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressedIndex))), 		&disp[12] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[13] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[13] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[14] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[14] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[15] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[15] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionIndex))), 		&disp[16] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionIndex))), 		&disp[16] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAveragesIndex))), 		&disp[17] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAveragesIndex))), 		&disp[17] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressedIndex))), 		&disp[18] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressedIndex))), 		&disp[18] );
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
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateIndex))), 		&disp[22] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateIndex))), 		&disp[22] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAveragesIndex))), 		&disp[23] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAveragesIndex))), 		&disp[23] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressedIndex))), 		&disp[24] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressedIndex))), 		&disp[24] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[25] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[25] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[26] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[26] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[27] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[27] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorIndex))), 		&disp[28] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorIndex))), 		&disp[28] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAveragesIndex))), 		&disp[29] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAveragesIndex))), 		&disp[29] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressedIndex))), 		&disp[30] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressedIndex))), 		&disp[30] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[31] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[31] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[32] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[32] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[33] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[33] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorGradientIndex))), 		&disp[34] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorGradientIndex))), 		&disp[34] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorGradient))), 		&disp[35] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorGradient))), 		&disp[35] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationIndex))), 		&disp[36] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationIndex))), 		&disp[36] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAveragesIndex))), 		&disp[37] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAveragesIndex))), 		&disp[37] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressedIndex))), 		&disp[38] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressedIndex))), 		&disp[38] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[39] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[39] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[40] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[40] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[41] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[41] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMinIndex))), 		&disp[42] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMinIndex))), 		&disp[42] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMaxIndex))), 		&disp[43] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMaxIndex))), 		&disp[43] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[44] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[44] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[45] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[45] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[46] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseAugmentationStatus[0]))), 		&disp[46] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[47] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._augmentationStatus))), 		&disp[47] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[48] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseCommunicationStatus[0]))), 		&disp[48] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[49] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._communicationStatus))), 		&disp[49] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[50] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._facewiseRefinementStatus[0]))), 		&disp[50] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[51] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementStatus))), 		&disp[51] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[52] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousRefinementStatus))), 		&disp[52] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[53] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._refinementFlag))), 		&disp[53] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._creation))), 		&disp[54] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._creation))), 		&disp[54] );
         #endif
         #ifdef MPI2
         MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[55] );
         #else
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[55] );
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
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]))), 		&disp[56] );
         disp[56] -= base;
         disp[56] += disp[0];
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



