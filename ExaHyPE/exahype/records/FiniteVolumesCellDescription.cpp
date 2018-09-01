#include "exahype/records/FiniteVolumesCellDescription.h"

#if defined(Parallel)
   exahype::records::FiniteVolumesCellDescription::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::FiniteVolumesCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const bool& oneRemoteBoundaryNeighbourIsOfTypeCell, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
   _solverNumber(solverNumber),
   _neighbourMergePerformed(neighbourMergePerformed),
   _hasCompletedTimeStep(hasCompletedTimeStep),
   _timeStepSize(timeStepSize),
   _timeStamp(timeStamp),
   _previousTimeStepSize(previousTimeStepSize),
   _previousTimeStamp(previousTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _extrapolatedSolution(extrapolatedSolution),
   _extrapolatedSolutionAverages(extrapolatedSolutionAverages),
   _extrapolatedSolutionCompressed(extrapolatedSolutionCompressed),
   _compressionState(compressionState),
   _bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
   _bytesPerDoFInSolution(bytesPerDoFInSolution),
   _bytesPerDoFInExtrapolatedSolution(bytesPerDoFInExtrapolatedSolution),
   _level(level),
   _offset(offset),
   _size(size),
   _oneRemoteBoundaryNeighbourIsOfTypeCell(oneRemoteBoundaryNeighbourIsOfTypeCell),
   _faceDataExchangeCounter(faceDataExchangeCounter),
   _type(type),
   _parentIndex(parentIndex),
   _refinementEvent(refinementEvent) {
      
   }
   
   exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription() {
      
   }
   
   
   exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords._hasCompletedTimeStep, persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._previousTimeStepSize, persistentRecords._previousTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._extrapolatedSolution, persistentRecords._extrapolatedSolutionAverages, persistentRecords._extrapolatedSolutionCompressed, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInExtrapolatedSolution, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell, persistentRecords._faceDataExchangeCounter, persistentRecords._type, persistentRecords._parentIndex, persistentRecords._refinementEvent) {
      
   }
   
   
   exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const bool& oneRemoteBoundaryNeighbourIsOfTypeCell, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
   _persistentRecords(solverNumber, neighbourMergePerformed, hasCompletedTimeStep, timeStepSize, timeStamp, previousTimeStepSize, previousTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, extrapolatedSolution, extrapolatedSolutionAverages, extrapolatedSolutionCompressed, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInExtrapolatedSolution, level, offset, size, oneRemoteBoundaryNeighbourIsOfTypeCell, faceDataExchangeCounter, type, parentIndex, refinementEvent) {
      
   }
   
   
   exahype::records::FiniteVolumesCellDescription::~FiniteVolumesCellDescription() { }
   
   std::string exahype::records::FiniteVolumesCellDescription::toString(const CompressionState& param) {
      switch (param) {
         case Uncompressed: return "Uncompressed";
         case CurrentlyProcessed: return "CurrentlyProcessed";
         case Compressed: return "Compressed";
      }
      return "undefined";
   }
   
   std::string exahype::records::FiniteVolumesCellDescription::getCompressionStateMapping() {
      return "CompressionState(Uncompressed=0,CurrentlyProcessed=1,Compressed=2)";
   }
   std::string exahype::records::FiniteVolumesCellDescription::toString(const RefinementEvent& param) {
      switch (param) {
         case None: return "None";
         case ErasingChildrenRequested: return "ErasingChildrenRequested";
         case ErasingChildren: return "ErasingChildren";
         case ChangeChildrenToDescendantsRequested: return "ChangeChildrenToDescendantsRequested";
         case ChangeChildrenToDescendants: return "ChangeChildrenToDescendants";
         case RefiningRequested: return "RefiningRequested";
         case Refining: return "Refining";
         case DeaugmentingChildrenRequestedTriggered: return "DeaugmentingChildrenRequestedTriggered";
         case DeaugmentingChildrenRequested: return "DeaugmentingChildrenRequested";
         case DeaugmentingChildren: return "DeaugmentingChildren";
         case AugmentingRequested: return "AugmentingRequested";
         case Augmenting: return "Augmenting";
      }
      return "undefined";
   }
   
   std::string exahype::records::FiniteVolumesCellDescription::getRefinementEventMapping() {
      return "RefinementEvent(None=0,ErasingChildrenRequested=1,ErasingChildren=2,ChangeChildrenToDescendantsRequested=3,ChangeChildrenToDescendants=4,RefiningRequested=5,Refining=6,DeaugmentingChildrenRequestedTriggered=7,DeaugmentingChildrenRequested=8,DeaugmentingChildren=9,AugmentingRequested=10,Augmenting=11)";
   }
   std::string exahype::records::FiniteVolumesCellDescription::toString(const Type& param) {
      switch (param) {
         case Erased: return "Erased";
         case Ancestor: return "Ancestor";
         case Cell: return "Cell";
         case Descendant: return "Descendant";
      }
      return "undefined";
   }
   
   std::string exahype::records::FiniteVolumesCellDescription::getTypeMapping() {
      return "Type(Erased=0,Ancestor=1,Cell=2,Descendant=3)";
   }
   
   
   std::string exahype::records::FiniteVolumesCellDescription::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::FiniteVolumesCellDescription::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "hasCompletedTimeStep:" << getHasCompletedTimeStep();
      out << ",";
      out << "timeStepSize:" << getTimeStepSize();
      out << ",";
      out << "timeStamp:" << getTimeStamp();
      out << ",";
      out << "previousTimeStepSize:" << getPreviousTimeStepSize();
      out << ",";
      out << "previousTimeStamp:" << getPreviousTimeStamp();
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
      out << "extrapolatedSolution:" << getExtrapolatedSolution();
      out << ",";
      out << "extrapolatedSolutionAverages:" << getExtrapolatedSolutionAverages();
      out << ",";
      out << "extrapolatedSolutionCompressed:" << getExtrapolatedSolutionCompressed();
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInExtrapolatedSolution:" << getBytesPerDoFInExtrapolatedSolution();
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
      out << "oneRemoteBoundaryNeighbourIsOfTypeCell:" << getOneRemoteBoundaryNeighbourIsOfTypeCell();
      out << ",";
      out << "faceDataExchangeCounter:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFaceDataExchangeCounter(i) << ",";
   }
   out << getFaceDataExchangeCounter(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out <<  ")";
   }
   
   
   exahype::records::FiniteVolumesCellDescription::PersistentRecords exahype::records::FiniteVolumesCellDescription::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::FiniteVolumesCellDescriptionPacked exahype::records::FiniteVolumesCellDescription::convert() const{
      return FiniteVolumesCellDescriptionPacked(
         getSolverNumber(),
         getNeighbourMergePerformed(),
         getHasCompletedTimeStep(),
         getTimeStepSize(),
         getTimeStamp(),
         getPreviousTimeStepSize(),
         getPreviousTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getExtrapolatedSolution(),
         getExtrapolatedSolutionAverages(),
         getExtrapolatedSolutionCompressed(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInExtrapolatedSolution(),
         getLevel(),
         getOffset(),
         getSize(),
         getOneRemoteBoundaryNeighbourIsOfTypeCell(),
         getFaceDataExchangeCounter(),
         getType(),
         getParentIndex(),
         getRefinementEvent()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::FiniteVolumesCellDescription::_log( "exahype::records::FiniteVolumesCellDescription" );
      
      MPI_Datatype exahype::records::FiniteVolumesCellDescription::Datatype = 0;
      MPI_Datatype exahype::records::FiniteVolumesCellDescription::FullDatatype = 0;
      
      
      void exahype::records::FiniteVolumesCellDescription::initDatatype() {
         {
            FiniteVolumesCellDescription dummyFiniteVolumesCellDescription[2];
            
            #ifdef MPI2
            const int Attributes = 27;
            #else
            const int Attributes = 28;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //neighbourMergePerformed
               , MPI_DOUBLE		 //timeStepSize
               , MPI_DOUBLE		 //timeStamp
               , MPI_DOUBLE		 //previousTimeStepSize
               , MPI_DOUBLE		 //previousTimeStamp
               , MPI_INT		 //solution
               , MPI_INT		 //solutionAverages
               , MPI_INT		 //solutionCompressed
               , MPI_INT		 //previousSolution
               , MPI_INT		 //previousSolutionAverages
               , MPI_INT		 //previousSolutionCompressed
               , MPI_INT		 //extrapolatedSolution
               , MPI_INT		 //extrapolatedSolutionAverages
               , MPI_INT		 //extrapolatedSolutionCompressed
               , MPI_INT		 //compressionState
               , MPI_INT		 //bytesPerDoFInPreviousSolution
               , MPI_INT		 //bytesPerDoFInSolution
               , MPI_INT		 //bytesPerDoFInExtrapolatedSolution
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_CXX_BOOL		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , MPI_INT		 //faceDataExchangeCounter
               , MPI_INT		 //type
               , MPI_INT		 //parentIndex
               , MPI_INT		 //refinementEvent
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
               , 1		 //timeStepSize
               , 1		 //timeStamp
               , 1		 //previousTimeStepSize
               , 1		 //previousTimeStamp
               , 1		 //solution
               , 1		 //solutionAverages
               , 1		 //solutionCompressed
               , 1		 //previousSolution
               , 1		 //previousSolutionAverages
               , 1		 //previousSolutionCompressed
               , 1		 //extrapolatedSolution
               , 1		 //extrapolatedSolutionAverages
               , 1		 //extrapolatedSolutionCompressed
               , 1		 //compressionState
               , 1		 //bytesPerDoFInPreviousSolution
               , 1		 //bytesPerDoFInSolution
               , 1		 //bytesPerDoFInExtrapolatedSolution
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , DIMENSIONS_TIMES_TWO		 //faceDataExchangeCounter
               , 1		 //type
               , 1		 //parentIndex
               , 1		 //refinementEvent
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[23] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[24] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[24] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[25] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[25] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[26] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[26] );
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
               assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescription)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescription));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[1]))), 		&disp[27] );
            disp[27] -= base;
            disp[27] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescription::Datatype );
            MPI_Type_commit( &FiniteVolumesCellDescription::Datatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescription::Datatype);
            MPI_Type_commit( &FiniteVolumesCellDescription::Datatype );
            #endif
            
         }
         {
            FiniteVolumesCellDescription dummyFiniteVolumesCellDescription[2];
            
            #ifdef MPI2
            const int Attributes = 28;
            #else
            const int Attributes = 29;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //neighbourMergePerformed
               , MPI_CXX_BOOL		 //hasCompletedTimeStep
               , MPI_DOUBLE		 //timeStepSize
               , MPI_DOUBLE		 //timeStamp
               , MPI_DOUBLE		 //previousTimeStepSize
               , MPI_DOUBLE		 //previousTimeStamp
               , MPI_INT		 //solution
               , MPI_INT		 //solutionAverages
               , MPI_INT		 //solutionCompressed
               , MPI_INT		 //previousSolution
               , MPI_INT		 //previousSolutionAverages
               , MPI_INT		 //previousSolutionCompressed
               , MPI_INT		 //extrapolatedSolution
               , MPI_INT		 //extrapolatedSolutionAverages
               , MPI_INT		 //extrapolatedSolutionCompressed
               , MPI_INT		 //compressionState
               , MPI_INT		 //bytesPerDoFInPreviousSolution
               , MPI_INT		 //bytesPerDoFInSolution
               , MPI_INT		 //bytesPerDoFInExtrapolatedSolution
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_CXX_BOOL		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , MPI_INT		 //faceDataExchangeCounter
               , MPI_INT		 //type
               , MPI_INT		 //parentIndex
               , MPI_INT		 //refinementEvent
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
               , 1		 //hasCompletedTimeStep
               , 1		 //timeStepSize
               , 1		 //timeStamp
               , 1		 //previousTimeStepSize
               , 1		 //previousTimeStamp
               , 1		 //solution
               , 1		 //solutionAverages
               , 1		 //solutionCompressed
               , 1		 //previousSolution
               , 1		 //previousSolutionAverages
               , 1		 //previousSolutionCompressed
               , 1		 //extrapolatedSolution
               , 1		 //extrapolatedSolutionAverages
               , 1		 //extrapolatedSolutionCompressed
               , 1		 //compressionState
               , 1		 //bytesPerDoFInPreviousSolution
               , 1		 //bytesPerDoFInSolution
               , 1		 //bytesPerDoFInExtrapolatedSolution
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , DIMENSIONS_TIMES_TWO		 //faceDataExchangeCounter
               , 1		 //type
               , 1		 //parentIndex
               , 1		 //refinementEvent
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._hasCompletedTimeStep))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._hasCompletedTimeStep))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[23] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[24] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[24] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[25] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[25] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[26] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[26] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[27] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[27] );
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
               assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescription)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescription));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[1]))), 		&disp[28] );
            disp[28] -= base;
            disp[28] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescription::FullDatatype );
            MPI_Type_commit( &FiniteVolumesCellDescription::FullDatatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescription::FullDatatype);
            MPI_Type_commit( &FiniteVolumesCellDescription::FullDatatype );
            #endif
            
         }
         
      }
      
      
      void exahype::records::FiniteVolumesCellDescription::shutdownDatatype() {
         MPI_Type_free( &FiniteVolumesCellDescription::Datatype );
         MPI_Type_free( &FiniteVolumesCellDescription::FullDatatype );
         
      }
      
      void exahype::records::FiniteVolumesCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
         // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescription " 
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
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescription "  
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
           msg << "testing for finished send task for exahype::records::FiniteVolumesCellDescription " 
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
             "exahype::records::FiniteVolumesCellDescription", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::FiniteVolumesCellDescription", 
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
      
      
      
      void exahype::records::FiniteVolumesCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node " 
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node " 
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
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: " 
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
        msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: " 
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
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node " 
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
      
      
      
      bool exahype::records::FiniteVolumesCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   
   exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords::PersistentRecords() {
      if ((12 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((12 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const bool& oneRemoteBoundaryNeighbourIsOfTypeCell, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
   _solverNumber(solverNumber),
   _neighbourMergePerformed(neighbourMergePerformed),
   _timeStepSize(timeStepSize),
   _timeStamp(timeStamp),
   _previousTimeStepSize(previousTimeStepSize),
   _previousTimeStamp(previousTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _extrapolatedSolution(extrapolatedSolution),
   _extrapolatedSolutionAverages(extrapolatedSolutionAverages),
   _extrapolatedSolutionCompressed(extrapolatedSolutionCompressed),
   _level(level),
   _offset(offset),
   _size(size),
   _oneRemoteBoundaryNeighbourIsOfTypeCell(oneRemoteBoundaryNeighbourIsOfTypeCell),
   _faceDataExchangeCounter(faceDataExchangeCounter),
   _type(type),
   _parentIndex(parentIndex),
   _refinementEvent(refinementEvent) {
      setHasCompletedTimeStep(hasCompletedTimeStep);
      setCompressionState(compressionState);
      setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
      setBytesPerDoFInSolution(bytesPerDoFInSolution);
      setBytesPerDoFInExtrapolatedSolution(bytesPerDoFInExtrapolatedSolution);
      if ((12 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((12 < (8 * sizeof(int))));
      
   }
   
   exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked() {
      if ((12 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((12 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords.getHasCompletedTimeStep(), persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._previousTimeStepSize, persistentRecords._previousTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._extrapolatedSolution, persistentRecords._extrapolatedSolutionAverages, persistentRecords._extrapolatedSolutionCompressed, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInExtrapolatedSolution(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell, persistentRecords._faceDataExchangeCounter, persistentRecords._type, persistentRecords._parentIndex, persistentRecords._refinementEvent) {
      if ((12 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((12 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const bool& oneRemoteBoundaryNeighbourIsOfTypeCell, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
   _persistentRecords(solverNumber, neighbourMergePerformed, hasCompletedTimeStep, timeStepSize, timeStamp, previousTimeStepSize, previousTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, extrapolatedSolution, extrapolatedSolutionAverages, extrapolatedSolutionCompressed, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInExtrapolatedSolution, level, offset, size, oneRemoteBoundaryNeighbourIsOfTypeCell, faceDataExchangeCounter, type, parentIndex, refinementEvent) {
      if ((12 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((12 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::FiniteVolumesCellDescriptionPacked::~FiniteVolumesCellDescriptionPacked() { }
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const CompressionState& param) {
      return exahype::records::FiniteVolumesCellDescription::toString(param);
   }
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::getCompressionStateMapping() {
      return exahype::records::FiniteVolumesCellDescription::getCompressionStateMapping();
   }
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const Type& param) {
      return exahype::records::FiniteVolumesCellDescription::toString(param);
   }
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::getTypeMapping() {
      return exahype::records::FiniteVolumesCellDescription::getTypeMapping();
   }
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const RefinementEvent& param) {
      return exahype::records::FiniteVolumesCellDescription::toString(param);
   }
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::getRefinementEventMapping() {
      return exahype::records::FiniteVolumesCellDescription::getRefinementEventMapping();
   }
   
   
   
   std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::FiniteVolumesCellDescriptionPacked::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "hasCompletedTimeStep:" << getHasCompletedTimeStep();
      out << ",";
      out << "timeStepSize:" << getTimeStepSize();
      out << ",";
      out << "timeStamp:" << getTimeStamp();
      out << ",";
      out << "previousTimeStepSize:" << getPreviousTimeStepSize();
      out << ",";
      out << "previousTimeStamp:" << getPreviousTimeStamp();
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
      out << "extrapolatedSolution:" << getExtrapolatedSolution();
      out << ",";
      out << "extrapolatedSolutionAverages:" << getExtrapolatedSolutionAverages();
      out << ",";
      out << "extrapolatedSolutionCompressed:" << getExtrapolatedSolutionCompressed();
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInExtrapolatedSolution:" << getBytesPerDoFInExtrapolatedSolution();
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
      out << "oneRemoteBoundaryNeighbourIsOfTypeCell:" << getOneRemoteBoundaryNeighbourIsOfTypeCell();
      out << ",";
      out << "faceDataExchangeCounter:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFaceDataExchangeCounter(i) << ",";
   }
   out << getFaceDataExchangeCounter(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out <<  ")";
   }
   
   
   exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords exahype::records::FiniteVolumesCellDescriptionPacked::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::FiniteVolumesCellDescription exahype::records::FiniteVolumesCellDescriptionPacked::convert() const{
      return FiniteVolumesCellDescription(
         getSolverNumber(),
         getNeighbourMergePerformed(),
         getHasCompletedTimeStep(),
         getTimeStepSize(),
         getTimeStamp(),
         getPreviousTimeStepSize(),
         getPreviousTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getExtrapolatedSolution(),
         getExtrapolatedSolutionAverages(),
         getExtrapolatedSolutionCompressed(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInExtrapolatedSolution(),
         getLevel(),
         getOffset(),
         getSize(),
         getOneRemoteBoundaryNeighbourIsOfTypeCell(),
         getFaceDataExchangeCounter(),
         getType(),
         getParentIndex(),
         getRefinementEvent()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::FiniteVolumesCellDescriptionPacked::_log( "exahype::records::FiniteVolumesCellDescriptionPacked" );
      
      MPI_Datatype exahype::records::FiniteVolumesCellDescriptionPacked::Datatype = 0;
      MPI_Datatype exahype::records::FiniteVolumesCellDescriptionPacked::FullDatatype = 0;
      
      
      void exahype::records::FiniteVolumesCellDescriptionPacked::initDatatype() {
         {
            FiniteVolumesCellDescriptionPacked dummyFiniteVolumesCellDescriptionPacked[2];
            
            #ifdef MPI2
            const int Attributes = 24;
            #else
            const int Attributes = 25;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //neighbourMergePerformed
               , MPI_DOUBLE		 //timeStepSize
               , MPI_DOUBLE		 //timeStamp
               , MPI_DOUBLE		 //previousTimeStepSize
               , MPI_DOUBLE		 //previousTimeStamp
               , MPI_INT		 //solution
               , MPI_INT		 //solutionAverages
               , MPI_INT		 //solutionCompressed
               , MPI_INT		 //previousSolution
               , MPI_INT		 //previousSolutionAverages
               , MPI_INT		 //previousSolutionCompressed
               , MPI_INT		 //extrapolatedSolution
               , MPI_INT		 //extrapolatedSolutionAverages
               , MPI_INT		 //extrapolatedSolutionCompressed
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_CXX_BOOL		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , MPI_INT		 //faceDataExchangeCounter
               , MPI_INT		 //type
               , MPI_INT		 //parentIndex
               , MPI_INT		 //refinementEvent
               , MPI_INT		 //_packedRecords0
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
               , 1		 //timeStepSize
               , 1		 //timeStamp
               , 1		 //previousTimeStepSize
               , 1		 //previousTimeStamp
               , 1		 //solution
               , 1		 //solutionAverages
               , 1		 //solutionCompressed
               , 1		 //previousSolution
               , 1		 //previousSolutionAverages
               , 1		 //previousSolutionCompressed
               , 1		 //extrapolatedSolution
               , 1		 //extrapolatedSolutionAverages
               , 1		 //extrapolatedSolutionCompressed
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , DIMENSIONS_TIMES_TWO		 //faceDataExchangeCounter
               , 1		 //type
               , 1		 //parentIndex
               , 1		 //refinementEvent
               , 1		 //_packedRecords0
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[23] );
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
               assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescriptionPacked)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescriptionPacked));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[1]))), 		&disp[24] );
            disp[24] -= base;
            disp[24] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescriptionPacked::Datatype );
            MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::Datatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescriptionPacked::Datatype);
            MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::Datatype );
            #endif
            
         }
         {
            FiniteVolumesCellDescriptionPacked dummyFiniteVolumesCellDescriptionPacked[2];
            
            #ifdef MPI2
            const int Attributes = 24;
            #else
            const int Attributes = 25;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //solverNumber
               , MPI_CXX_BOOL		 //neighbourMergePerformed
               , MPI_DOUBLE		 //timeStepSize
               , MPI_DOUBLE		 //timeStamp
               , MPI_DOUBLE		 //previousTimeStepSize
               , MPI_DOUBLE		 //previousTimeStamp
               , MPI_INT		 //solution
               , MPI_INT		 //solutionAverages
               , MPI_INT		 //solutionCompressed
               , MPI_INT		 //previousSolution
               , MPI_INT		 //previousSolutionAverages
               , MPI_INT		 //previousSolutionCompressed
               , MPI_INT		 //extrapolatedSolution
               , MPI_INT		 //extrapolatedSolutionAverages
               , MPI_INT		 //extrapolatedSolutionCompressed
               , MPI_INT		 //level
               , MPI_DOUBLE		 //offset
               , MPI_DOUBLE		 //size
               , MPI_CXX_BOOL		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , MPI_INT		 //faceDataExchangeCounter
               , MPI_INT		 //type
               , MPI_INT		 //parentIndex
               , MPI_INT		 //refinementEvent
               , MPI_INT		 //_packedRecords0
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //solverNumber
               , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
               , 1		 //timeStepSize
               , 1		 //timeStamp
               , 1		 //previousTimeStepSize
               , 1		 //previousTimeStamp
               , 1		 //solution
               , 1		 //solutionAverages
               , 1		 //solutionCompressed
               , 1		 //previousSolution
               , 1		 //previousSolutionAverages
               , 1		 //previousSolutionCompressed
               , 1		 //extrapolatedSolution
               , 1		 //extrapolatedSolutionAverages
               , 1		 //extrapolatedSolutionCompressed
               , 1		 //level
               , DIMENSIONS		 //offset
               , DIMENSIONS		 //size
               , 1		 //oneRemoteBoundaryNeighbourIsOfTypeCell
               , DIMENSIONS_TIMES_TWO		 //faceDataExchangeCounter
               , 1		 //type
               , 1		 //parentIndex
               , 1		 //refinementEvent
               , 1		 //_packedRecords0
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[18] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._oneRemoteBoundaryNeighbourIsOfTypeCell))), 		&disp[18] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[19] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[19] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[20] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[20] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[21] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[21] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[22] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[22] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[23] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[23] );
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
               assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescriptionPacked)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescriptionPacked));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[1]))), 		&disp[24] );
            disp[24] -= base;
            disp[24] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescriptionPacked::FullDatatype );
            MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::FullDatatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescriptionPacked::FullDatatype);
            MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::FullDatatype );
            #endif
            
         }
         
      }
      
      
      void exahype::records::FiniteVolumesCellDescriptionPacked::shutdownDatatype() {
         MPI_Type_free( &FiniteVolumesCellDescriptionPacked::Datatype );
         MPI_Type_free( &FiniteVolumesCellDescriptionPacked::FullDatatype );
         
      }
      
      void exahype::records::FiniteVolumesCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
         // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescriptionPacked " 
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
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescriptionPacked "  
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
           msg << "testing for finished send task for exahype::records::FiniteVolumesCellDescriptionPacked " 
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
             "exahype::records::FiniteVolumesCellDescriptionPacked", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::FiniteVolumesCellDescriptionPacked", 
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
      
      
      
      void exahype::records::FiniteVolumesCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node " 
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node " 
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
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: " 
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
        msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: " 
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
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node " 
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
      
      
      
      bool exahype::records::FiniteVolumesCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
      exahype::records::FiniteVolumesCellDescription::PersistentRecords::PersistentRecords() {
         
      }
      
      
      exahype::records::FiniteVolumesCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
      _solverNumber(solverNumber),
      _neighbourMergePerformed(neighbourMergePerformed),
      _hasCompletedTimeStep(hasCompletedTimeStep),
      _timeStepSize(timeStepSize),
      _timeStamp(timeStamp),
      _previousTimeStepSize(previousTimeStepSize),
      _previousTimeStamp(previousTimeStamp),
      _solution(solution),
      _solutionAverages(solutionAverages),
      _solutionCompressed(solutionCompressed),
      _previousSolution(previousSolution),
      _previousSolutionAverages(previousSolutionAverages),
      _previousSolutionCompressed(previousSolutionCompressed),
      _extrapolatedSolution(extrapolatedSolution),
      _extrapolatedSolutionAverages(extrapolatedSolutionAverages),
      _extrapolatedSolutionCompressed(extrapolatedSolutionCompressed),
      _compressionState(compressionState),
      _bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
      _bytesPerDoFInSolution(bytesPerDoFInSolution),
      _bytesPerDoFInExtrapolatedSolution(bytesPerDoFInExtrapolatedSolution),
      _level(level),
      _offset(offset),
      _size(size),
      _type(type),
      _parentIndex(parentIndex),
      _refinementEvent(refinementEvent) {
         
      }
      
      exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription() {
         
      }
      
      
      exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords._hasCompletedTimeStep, persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._previousTimeStepSize, persistentRecords._previousTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._extrapolatedSolution, persistentRecords._extrapolatedSolutionAverages, persistentRecords._extrapolatedSolutionCompressed, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInExtrapolatedSolution, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._type, persistentRecords._parentIndex, persistentRecords._refinementEvent) {
         
      }
      
      
      exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
      _persistentRecords(solverNumber, neighbourMergePerformed, hasCompletedTimeStep, timeStepSize, timeStamp, previousTimeStepSize, previousTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, extrapolatedSolution, extrapolatedSolutionAverages, extrapolatedSolutionCompressed, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInExtrapolatedSolution, level, offset, size, type, parentIndex, refinementEvent) {
         
      }
      
      
      exahype::records::FiniteVolumesCellDescription::~FiniteVolumesCellDescription() { }
      
      std::string exahype::records::FiniteVolumesCellDescription::toString(const CompressionState& param) {
         switch (param) {
            case Uncompressed: return "Uncompressed";
            case CurrentlyProcessed: return "CurrentlyProcessed";
            case Compressed: return "Compressed";
         }
         return "undefined";
      }
      
      std::string exahype::records::FiniteVolumesCellDescription::getCompressionStateMapping() {
         return "CompressionState(Uncompressed=0,CurrentlyProcessed=1,Compressed=2)";
      }
      std::string exahype::records::FiniteVolumesCellDescription::toString(const RefinementEvent& param) {
         switch (param) {
            case None: return "None";
            case ErasingChildrenRequested: return "ErasingChildrenRequested";
            case ErasingChildren: return "ErasingChildren";
            case ChangeChildrenToDescendantsRequested: return "ChangeChildrenToDescendantsRequested";
            case ChangeChildrenToDescendants: return "ChangeChildrenToDescendants";
            case RefiningRequested: return "RefiningRequested";
            case Refining: return "Refining";
            case DeaugmentingChildrenRequestedTriggered: return "DeaugmentingChildrenRequestedTriggered";
            case DeaugmentingChildrenRequested: return "DeaugmentingChildrenRequested";
            case DeaugmentingChildren: return "DeaugmentingChildren";
            case AugmentingRequested: return "AugmentingRequested";
            case Augmenting: return "Augmenting";
         }
         return "undefined";
      }
      
      std::string exahype::records::FiniteVolumesCellDescription::getRefinementEventMapping() {
         return "RefinementEvent(None=0,ErasingChildrenRequested=1,ErasingChildren=2,ChangeChildrenToDescendantsRequested=3,ChangeChildrenToDescendants=4,RefiningRequested=5,Refining=6,DeaugmentingChildrenRequestedTriggered=7,DeaugmentingChildrenRequested=8,DeaugmentingChildren=9,AugmentingRequested=10,Augmenting=11)";
      }
      std::string exahype::records::FiniteVolumesCellDescription::toString(const Type& param) {
         switch (param) {
            case Erased: return "Erased";
            case Ancestor: return "Ancestor";
            case Cell: return "Cell";
            case Descendant: return "Descendant";
         }
         return "undefined";
      }
      
      std::string exahype::records::FiniteVolumesCellDescription::getTypeMapping() {
         return "Type(Erased=0,Ancestor=1,Cell=2,Descendant=3)";
      }
      
      
      std::string exahype::records::FiniteVolumesCellDescription::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::FiniteVolumesCellDescription::toString (std::ostream& out) const {
         out << "("; 
         out << "solverNumber:" << getSolverNumber();
         out << ",";
         out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "hasCompletedTimeStep:" << getHasCompletedTimeStep();
         out << ",";
         out << "timeStepSize:" << getTimeStepSize();
         out << ",";
         out << "timeStamp:" << getTimeStamp();
         out << ",";
         out << "previousTimeStepSize:" << getPreviousTimeStepSize();
         out << ",";
         out << "previousTimeStamp:" << getPreviousTimeStamp();
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
         out << "extrapolatedSolution:" << getExtrapolatedSolution();
         out << ",";
         out << "extrapolatedSolutionAverages:" << getExtrapolatedSolutionAverages();
         out << ",";
         out << "extrapolatedSolutionCompressed:" << getExtrapolatedSolutionCompressed();
         out << ",";
         out << "compressionState:" << toString(getCompressionState());
         out << ",";
         out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
         out << ",";
         out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
         out << ",";
         out << "bytesPerDoFInExtrapolatedSolution:" << getBytesPerDoFInExtrapolatedSolution();
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
         out << "type:" << toString(getType());
         out << ",";
         out << "parentIndex:" << getParentIndex();
         out << ",";
         out << "refinementEvent:" << toString(getRefinementEvent());
         out <<  ")";
      }
      
      
      exahype::records::FiniteVolumesCellDescription::PersistentRecords exahype::records::FiniteVolumesCellDescription::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::FiniteVolumesCellDescriptionPacked exahype::records::FiniteVolumesCellDescription::convert() const{
         return FiniteVolumesCellDescriptionPacked(
            getSolverNumber(),
            getNeighbourMergePerformed(),
            getHasCompletedTimeStep(),
            getTimeStepSize(),
            getTimeStamp(),
            getPreviousTimeStepSize(),
            getPreviousTimeStamp(),
            getSolution(),
            getSolutionAverages(),
            getSolutionCompressed(),
            getPreviousSolution(),
            getPreviousSolutionAverages(),
            getPreviousSolutionCompressed(),
            getExtrapolatedSolution(),
            getExtrapolatedSolutionAverages(),
            getExtrapolatedSolutionCompressed(),
            getCompressionState(),
            getBytesPerDoFInPreviousSolution(),
            getBytesPerDoFInSolution(),
            getBytesPerDoFInExtrapolatedSolution(),
            getLevel(),
            getOffset(),
            getSize(),
            getType(),
            getParentIndex(),
            getRefinementEvent()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::FiniteVolumesCellDescription::_log( "exahype::records::FiniteVolumesCellDescription" );
         
         MPI_Datatype exahype::records::FiniteVolumesCellDescription::Datatype = 0;
         MPI_Datatype exahype::records::FiniteVolumesCellDescription::FullDatatype = 0;
         
         
         void exahype::records::FiniteVolumesCellDescription::initDatatype() {
            {
               FiniteVolumesCellDescription dummyFiniteVolumesCellDescription[2];
               
               #ifdef MPI2
               const int Attributes = 25;
               #else
               const int Attributes = 26;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_CXX_BOOL		 //neighbourMergePerformed
                  , MPI_DOUBLE		 //timeStepSize
                  , MPI_DOUBLE		 //timeStamp
                  , MPI_DOUBLE		 //previousTimeStepSize
                  , MPI_DOUBLE		 //previousTimeStamp
                  , MPI_INT		 //solution
                  , MPI_INT		 //solutionAverages
                  , MPI_INT		 //solutionCompressed
                  , MPI_INT		 //previousSolution
                  , MPI_INT		 //previousSolutionAverages
                  , MPI_INT		 //previousSolutionCompressed
                  , MPI_INT		 //extrapolatedSolution
                  , MPI_INT		 //extrapolatedSolutionAverages
                  , MPI_INT		 //extrapolatedSolutionCompressed
                  , MPI_INT		 //compressionState
                  , MPI_INT		 //bytesPerDoFInPreviousSolution
                  , MPI_INT		 //bytesPerDoFInSolution
                  , MPI_INT		 //bytesPerDoFInExtrapolatedSolution
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_INT		 //type
                  , MPI_INT		 //parentIndex
                  , MPI_INT		 //refinementEvent
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
                  , 1		 //timeStepSize
                  , 1		 //timeStamp
                  , 1		 //previousTimeStepSize
                  , 1		 //previousTimeStamp
                  , 1		 //solution
                  , 1		 //solutionAverages
                  , 1		 //solutionCompressed
                  , 1		 //previousSolution
                  , 1		 //previousSolutionAverages
                  , 1		 //previousSolutionCompressed
                  , 1		 //extrapolatedSolution
                  , 1		 //extrapolatedSolutionAverages
                  , 1		 //extrapolatedSolutionCompressed
                  , 1		 //compressionState
                  , 1		 //bytesPerDoFInPreviousSolution
                  , 1		 //bytesPerDoFInSolution
                  , 1		 //bytesPerDoFInExtrapolatedSolution
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //type
                  , 1		 //parentIndex
                  , 1		 //refinementEvent
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[24] );
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
                  assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescription)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescription));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[1]))), 		&disp[25] );
               disp[25] -= base;
               disp[25] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescription::Datatype );
               MPI_Type_commit( &FiniteVolumesCellDescription::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescription::Datatype);
               MPI_Type_commit( &FiniteVolumesCellDescription::Datatype );
               #endif
               
            }
            {
               FiniteVolumesCellDescription dummyFiniteVolumesCellDescription[2];
               
               #ifdef MPI2
               const int Attributes = 26;
               #else
               const int Attributes = 27;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_CXX_BOOL		 //neighbourMergePerformed
                  , MPI_CXX_BOOL		 //hasCompletedTimeStep
                  , MPI_DOUBLE		 //timeStepSize
                  , MPI_DOUBLE		 //timeStamp
                  , MPI_DOUBLE		 //previousTimeStepSize
                  , MPI_DOUBLE		 //previousTimeStamp
                  , MPI_INT		 //solution
                  , MPI_INT		 //solutionAverages
                  , MPI_INT		 //solutionCompressed
                  , MPI_INT		 //previousSolution
                  , MPI_INT		 //previousSolutionAverages
                  , MPI_INT		 //previousSolutionCompressed
                  , MPI_INT		 //extrapolatedSolution
                  , MPI_INT		 //extrapolatedSolutionAverages
                  , MPI_INT		 //extrapolatedSolutionCompressed
                  , MPI_INT		 //compressionState
                  , MPI_INT		 //bytesPerDoFInPreviousSolution
                  , MPI_INT		 //bytesPerDoFInSolution
                  , MPI_INT		 //bytesPerDoFInExtrapolatedSolution
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_INT		 //type
                  , MPI_INT		 //parentIndex
                  , MPI_INT		 //refinementEvent
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
                  , 1		 //hasCompletedTimeStep
                  , 1		 //timeStepSize
                  , 1		 //timeStamp
                  , 1		 //previousTimeStepSize
                  , 1		 //previousTimeStamp
                  , 1		 //solution
                  , 1		 //solutionAverages
                  , 1		 //solutionCompressed
                  , 1		 //previousSolution
                  , 1		 //previousSolutionAverages
                  , 1		 //previousSolutionCompressed
                  , 1		 //extrapolatedSolution
                  , 1		 //extrapolatedSolutionAverages
                  , 1		 //extrapolatedSolutionCompressed
                  , 1		 //compressionState
                  , 1		 //bytesPerDoFInPreviousSolution
                  , 1		 //bytesPerDoFInSolution
                  , 1		 //bytesPerDoFInExtrapolatedSolution
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //type
                  , 1		 //parentIndex
                  , 1		 //refinementEvent
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._hasCompletedTimeStep))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._hasCompletedTimeStep))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStepSize))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousTimeStamp))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolution))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolution))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._compressionState))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedSolution))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._parentIndex))), 		&disp[24] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[25] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[25] );
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
                  assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescription)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescription));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[1]))), 		&disp[26] );
               disp[26] -= base;
               disp[26] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescription::FullDatatype );
               MPI_Type_commit( &FiniteVolumesCellDescription::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescription::FullDatatype);
               MPI_Type_commit( &FiniteVolumesCellDescription::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::FiniteVolumesCellDescription::shutdownDatatype() {
            MPI_Type_free( &FiniteVolumesCellDescription::Datatype );
            MPI_Type_free( &FiniteVolumesCellDescription::FullDatatype );
            
         }
         
         void exahype::records::FiniteVolumesCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
            // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescription " 
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
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescription "  
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
           msg << "testing for finished send task for exahype::records::FiniteVolumesCellDescription " 
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
             "exahype::records::FiniteVolumesCellDescription", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::FiniteVolumesCellDescription", 
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
         
         
         
         void exahype::records::FiniteVolumesCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node " 
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node " 
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
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: " 
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
        msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: " 
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
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescription", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node " 
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
         
         
         
         bool exahype::records::FiniteVolumesCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
      
      
      exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords::PersistentRecords() {
         if ((12 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((12 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
      _solverNumber(solverNumber),
      _neighbourMergePerformed(neighbourMergePerformed),
      _timeStepSize(timeStepSize),
      _timeStamp(timeStamp),
      _previousTimeStepSize(previousTimeStepSize),
      _previousTimeStamp(previousTimeStamp),
      _solution(solution),
      _solutionAverages(solutionAverages),
      _solutionCompressed(solutionCompressed),
      _previousSolution(previousSolution),
      _previousSolutionAverages(previousSolutionAverages),
      _previousSolutionCompressed(previousSolutionCompressed),
      _extrapolatedSolution(extrapolatedSolution),
      _extrapolatedSolutionAverages(extrapolatedSolutionAverages),
      _extrapolatedSolutionCompressed(extrapolatedSolutionCompressed),
      _level(level),
      _offset(offset),
      _size(size),
      _type(type),
      _parentIndex(parentIndex),
      _refinementEvent(refinementEvent) {
         setHasCompletedTimeStep(hasCompletedTimeStep);
         setCompressionState(compressionState);
         setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
         setBytesPerDoFInSolution(bytesPerDoFInSolution);
         setBytesPerDoFInExtrapolatedSolution(bytesPerDoFInExtrapolatedSolution);
         if ((12 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((12 < (8 * sizeof(int))));
         
      }
      
      exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked() {
         if ((12 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((12 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._solverNumber, persistentRecords._neighbourMergePerformed, persistentRecords.getHasCompletedTimeStep(), persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._previousTimeStepSize, persistentRecords._previousTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._extrapolatedSolution, persistentRecords._extrapolatedSolutionAverages, persistentRecords._extrapolatedSolutionCompressed, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInExtrapolatedSolution(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._type, persistentRecords._parentIndex, persistentRecords._refinementEvent) {
         if ((12 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((12 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& neighbourMergePerformed, const bool& hasCompletedTimeStep, const double& timeStepSize, const double& timeStamp, const double& previousTimeStepSize, const double& previousTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& extrapolatedSolution, const int& extrapolatedSolutionAverages, const int& extrapolatedSolutionCompressed, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInExtrapolatedSolution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type, const int& parentIndex, const RefinementEvent& refinementEvent):
      _persistentRecords(solverNumber, neighbourMergePerformed, hasCompletedTimeStep, timeStepSize, timeStamp, previousTimeStepSize, previousTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, extrapolatedSolution, extrapolatedSolutionAverages, extrapolatedSolutionCompressed, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInExtrapolatedSolution, level, offset, size, type, parentIndex, refinementEvent) {
         if ((12 >= (8 * sizeof(int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((12 < (8 * sizeof(int))));
         
      }
      
      
      exahype::records::FiniteVolumesCellDescriptionPacked::~FiniteVolumesCellDescriptionPacked() { }
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const CompressionState& param) {
         return exahype::records::FiniteVolumesCellDescription::toString(param);
      }
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::getCompressionStateMapping() {
         return exahype::records::FiniteVolumesCellDescription::getCompressionStateMapping();
      }
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const Type& param) {
         return exahype::records::FiniteVolumesCellDescription::toString(param);
      }
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::getTypeMapping() {
         return exahype::records::FiniteVolumesCellDescription::getTypeMapping();
      }
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const RefinementEvent& param) {
         return exahype::records::FiniteVolumesCellDescription::toString(param);
      }
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::getRefinementEventMapping() {
         return exahype::records::FiniteVolumesCellDescription::getRefinementEventMapping();
      }
      
      
      
      std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::FiniteVolumesCellDescriptionPacked::toString (std::ostream& out) const {
         out << "("; 
         out << "solverNumber:" << getSolverNumber();
         out << ",";
         out << "neighbourMergePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getNeighbourMergePerformed(i) << ",";
   }
   out << getNeighbourMergePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
         out << ",";
         out << "hasCompletedTimeStep:" << getHasCompletedTimeStep();
         out << ",";
         out << "timeStepSize:" << getTimeStepSize();
         out << ",";
         out << "timeStamp:" << getTimeStamp();
         out << ",";
         out << "previousTimeStepSize:" << getPreviousTimeStepSize();
         out << ",";
         out << "previousTimeStamp:" << getPreviousTimeStamp();
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
         out << "extrapolatedSolution:" << getExtrapolatedSolution();
         out << ",";
         out << "extrapolatedSolutionAverages:" << getExtrapolatedSolutionAverages();
         out << ",";
         out << "extrapolatedSolutionCompressed:" << getExtrapolatedSolutionCompressed();
         out << ",";
         out << "compressionState:" << toString(getCompressionState());
         out << ",";
         out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
         out << ",";
         out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
         out << ",";
         out << "bytesPerDoFInExtrapolatedSolution:" << getBytesPerDoFInExtrapolatedSolution();
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
         out << "type:" << toString(getType());
         out << ",";
         out << "parentIndex:" << getParentIndex();
         out << ",";
         out << "refinementEvent:" << toString(getRefinementEvent());
         out <<  ")";
      }
      
      
      exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords exahype::records::FiniteVolumesCellDescriptionPacked::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::FiniteVolumesCellDescription exahype::records::FiniteVolumesCellDescriptionPacked::convert() const{
         return FiniteVolumesCellDescription(
            getSolverNumber(),
            getNeighbourMergePerformed(),
            getHasCompletedTimeStep(),
            getTimeStepSize(),
            getTimeStamp(),
            getPreviousTimeStepSize(),
            getPreviousTimeStamp(),
            getSolution(),
            getSolutionAverages(),
            getSolutionCompressed(),
            getPreviousSolution(),
            getPreviousSolutionAverages(),
            getPreviousSolutionCompressed(),
            getExtrapolatedSolution(),
            getExtrapolatedSolutionAverages(),
            getExtrapolatedSolutionCompressed(),
            getCompressionState(),
            getBytesPerDoFInPreviousSolution(),
            getBytesPerDoFInSolution(),
            getBytesPerDoFInExtrapolatedSolution(),
            getLevel(),
            getOffset(),
            getSize(),
            getType(),
            getParentIndex(),
            getRefinementEvent()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::FiniteVolumesCellDescriptionPacked::_log( "exahype::records::FiniteVolumesCellDescriptionPacked" );
         
         MPI_Datatype exahype::records::FiniteVolumesCellDescriptionPacked::Datatype = 0;
         MPI_Datatype exahype::records::FiniteVolumesCellDescriptionPacked::FullDatatype = 0;
         
         
         void exahype::records::FiniteVolumesCellDescriptionPacked::initDatatype() {
            {
               FiniteVolumesCellDescriptionPacked dummyFiniteVolumesCellDescriptionPacked[2];
               
               #ifdef MPI2
               const int Attributes = 22;
               #else
               const int Attributes = 23;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_CXX_BOOL		 //neighbourMergePerformed
                  , MPI_DOUBLE		 //timeStepSize
                  , MPI_DOUBLE		 //timeStamp
                  , MPI_DOUBLE		 //previousTimeStepSize
                  , MPI_DOUBLE		 //previousTimeStamp
                  , MPI_INT		 //solution
                  , MPI_INT		 //solutionAverages
                  , MPI_INT		 //solutionCompressed
                  , MPI_INT		 //previousSolution
                  , MPI_INT		 //previousSolutionAverages
                  , MPI_INT		 //previousSolutionCompressed
                  , MPI_INT		 //extrapolatedSolution
                  , MPI_INT		 //extrapolatedSolutionAverages
                  , MPI_INT		 //extrapolatedSolutionCompressed
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_INT		 //type
                  , MPI_INT		 //parentIndex
                  , MPI_INT		 //refinementEvent
                  , MPI_INT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
                  , 1		 //timeStepSize
                  , 1		 //timeStamp
                  , 1		 //previousTimeStepSize
                  , 1		 //previousTimeStamp
                  , 1		 //solution
                  , 1		 //solutionAverages
                  , 1		 //solutionCompressed
                  , 1		 //previousSolution
                  , 1		 //previousSolutionAverages
                  , 1		 //previousSolutionCompressed
                  , 1		 //extrapolatedSolution
                  , 1		 //extrapolatedSolutionAverages
                  , 1		 //extrapolatedSolutionCompressed
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //type
                  , 1		 //parentIndex
                  , 1		 //refinementEvent
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[21] );
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
                  assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescriptionPacked)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescriptionPacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[1]))), 		&disp[22] );
               disp[22] -= base;
               disp[22] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescriptionPacked::Datatype );
               MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescriptionPacked::Datatype);
               MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::Datatype );
               #endif
               
            }
            {
               FiniteVolumesCellDescriptionPacked dummyFiniteVolumesCellDescriptionPacked[2];
               
               #ifdef MPI2
               const int Attributes = 22;
               #else
               const int Attributes = 23;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //solverNumber
                  , MPI_CXX_BOOL		 //neighbourMergePerformed
                  , MPI_DOUBLE		 //timeStepSize
                  , MPI_DOUBLE		 //timeStamp
                  , MPI_DOUBLE		 //previousTimeStepSize
                  , MPI_DOUBLE		 //previousTimeStamp
                  , MPI_INT		 //solution
                  , MPI_INT		 //solutionAverages
                  , MPI_INT		 //solutionCompressed
                  , MPI_INT		 //previousSolution
                  , MPI_INT		 //previousSolutionAverages
                  , MPI_INT		 //previousSolutionCompressed
                  , MPI_INT		 //extrapolatedSolution
                  , MPI_INT		 //extrapolatedSolutionAverages
                  , MPI_INT		 //extrapolatedSolutionCompressed
                  , MPI_INT		 //level
                  , MPI_DOUBLE		 //offset
                  , MPI_DOUBLE		 //size
                  , MPI_INT		 //type
                  , MPI_INT		 //parentIndex
                  , MPI_INT		 //refinementEvent
                  , MPI_INT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //solverNumber
                  , DIMENSIONS_TIMES_TWO		 //neighbourMergePerformed
                  , 1		 //timeStepSize
                  , 1		 //timeStamp
                  , 1		 //previousTimeStepSize
                  , 1		 //previousTimeStamp
                  , 1		 //solution
                  , 1		 //solutionAverages
                  , 1		 //solutionCompressed
                  , 1		 //previousSolution
                  , 1		 //previousSolutionAverages
                  , 1		 //previousSolutionCompressed
                  , 1		 //extrapolatedSolution
                  , 1		 //extrapolatedSolutionAverages
                  , 1		 //extrapolatedSolutionCompressed
                  , 1		 //level
                  , DIMENSIONS		 //offset
                  , DIMENSIONS		 //size
                  , 1		 //type
                  , 1		 //parentIndex
                  , 1		 //refinementEvent
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._neighbourMergePerformed))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStepSize))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousTimeStamp))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolution))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionAverages))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._extrapolatedSolutionCompressed))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._refinementEvent))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[21] );
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
                  assertion4(disp[i]<static_cast<int>(sizeof(FiniteVolumesCellDescriptionPacked)), i, disp[i], Attributes, sizeof(FiniteVolumesCellDescriptionPacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[1]))), 		&disp[22] );
               disp[22] -= base;
               disp[22] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &FiniteVolumesCellDescriptionPacked::FullDatatype );
               MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescriptionPacked::FullDatatype);
               MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::FiniteVolumesCellDescriptionPacked::shutdownDatatype() {
            MPI_Type_free( &FiniteVolumesCellDescriptionPacked::Datatype );
            MPI_Type_free( &FiniteVolumesCellDescriptionPacked::FullDatatype );
            
         }
         
         void exahype::records::FiniteVolumesCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
            // ============================= 
// start injected snippet/aspect 
// ============================= 
switch (mode) { 
  case ExchangeMode::Blocking: 
    {
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator()); 
       if  (result!=MPI_SUCCESS) { 
         std::ostringstream msg; 
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescriptionPacked " 
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
         msg << "was not able to send message exahype::records::FiniteVolumesCellDescriptionPacked "  
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
           msg << "testing for finished send task for exahype::records::FiniteVolumesCellDescriptionPacked " 
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
             "exahype::records::FiniteVolumesCellDescriptionPacked", 
             "send(int)", destination,tag,1 
           ); 
           triggeredTimeoutWarning = true; 
         } 
         if ( 
           tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
           (clock()>timeOutShutdown) 
         ) { 
           tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
             "exahype::records::FiniteVolumesCellDescriptionPacked", 
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
         
         
         
         void exahype::records::FiniteVolumesCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, ExchangeMode mode) {
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node " 
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
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node " 
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
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Test( sendRequestHandle, &flag, source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
        if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: " 
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
        msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: " 
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
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
          triggeredTimeoutWarning = true; 
        } 
        if ( 
          tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() && 
          (clock()>timeOutShutdown) 
        ) { 
          tarch::parallel::Node::getInstance().triggerDeadlockTimeOut( 
            "exahype::records::FiniteVolumesCellDescriptionPacked", 
            "receive(int)", source,tag,1 
          ); 
        } 
        tarch::parallel::Node::getInstance().receiveDanglingMessages(); 
        result = MPI_Iprobe(source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &flag, MPI_STATUS_IGNORE ); 
         if (result!=MPI_SUCCESS) { 
          std::ostringstream msg; 
          msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: " 
              << tarch::parallel::MPIReturnValueToString(result); 
          _log.error("receive(int)", msg.str() ); 
        } 
      } 
      result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), source==MPI_ANY_SOURCE ? &status : MPI_STATUS_IGNORE ); 
      if ( result != MPI_SUCCESS ) { 
        std::ostringstream msg; 
        msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node " 
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
         
         
         
         bool exahype::records::FiniteVolumesCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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


