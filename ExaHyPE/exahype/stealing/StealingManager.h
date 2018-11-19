#if !defined(_EXAHYPE_STEALING_STEALINGMANAGER_H_) && defined(SharedTBB)  && defined(Parallel)
#define _EXAHYPE_STEALING_STEALINGMANAGER_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tbb/concurrent_hash_map.h"
#include "tbb/concurrent_queue.h"

#include "exahype/solvers/Solver.h"

#include <mpi.h>
#include <atomic>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>

namespace exahype {
  namespace stealing {
    class StealingManager;
    /*
     * Different MPI request types for prioritizing message requests.
     */
    enum class RequestType {
      send = 0,
      sendBack = 1,
      receive = 2,
      receiveBack = 3
    };
  }
}

namespace exahype {
  namespace solvers {
    class Solver;
  }
}

/*
* This stealing manager manages the asynchronous MPI requests
* created e.g., when a task is stolen. In addition, the stealing
* manager serves as a connecting point to the ADERDGSolver
* for making decisions on whether a task should be stolen or not.
*/
class exahype::stealing::StealingManager {
  private:
    static tarch::logging::Log _log;

    /*
     *  MPI_Request are mapped to an internal id.
     */
    std::atomic<int> _nextRequestId;
    /*
     *  MPI_Request groups (e.g. all requests belonging to a task send) are mapped to
     *  an internal group id.
     */
    std::atomic<int> _nextGroupId;

    // queues for each message type
    tbb::concurrent_queue<MPI_Request> _requests[4];

    tbb::concurrent_hash_map<int, MPI_Request> _idToRequest[4];
    tbb::concurrent_hash_map<int, int> _requestToGroup[4];
    tbb::concurrent_hash_map<int, int> _outstandingReqsForGroup[4];
    tbb::concurrent_hash_map<int, int> _remoteRanksForGroup[4];
    tbb::concurrent_hash_map<int, int> _remoteTagsForGroup[4];
    tbb::concurrent_hash_map<int, std::function<void(exahype::solvers::Solver*, int , int)>> _handlers[4];
    tbb::concurrent_hash_map<int, exahype::solvers::Solver*> _solvers[4];
    tarch::multicore::BooleanSemaphore _semaphore;

    StealingManager();
    /*
     * This method makes progress on all current requests of the given request type.
     */
    void progressRequestsOfType(RequestType type);
    /*
     * Maps a request type to an integer that defines message queue.
     */
    int requestTypeToMap(RequestType requestType);
    /*
     * This method pop's all current requests of a given type from the request queue and
     * inserts them into an array on which MPI can make progress.
     */
    void createRequestArray(
        RequestType type,
		std::vector<MPI_Request> &requests,
		std::unordered_map<int, int> &vecIdToReqId);

    // for all stealing-related communication, a separate MPI communicator is used (needs to be created in runGlobalStep())
    MPI_Comm _stealingComm;
    MPI_Comm _stealingCommMapped;

    /*
     * The request handler job aims to distribute the work that is to be done
     * when a request group is finished evenly among the TBB worker threads.
     */
    class RequestHandlerJob
    {
      private:
        std::function<void(exahype::solvers::Solver*, int, int)> _handleRequest;
        exahype::solvers::Solver* _solver;
        int _tag;
        int _remoteRank;
      public:
        RequestHandlerJob(
            std::function<void(exahype::solvers::Solver*, int, int)> handleRequest,
            exahype::solvers::Solver* solver,
            int tag,
            int remoteRank);
        bool operator()();
    };

  public:
    /*
     * Submit a group of MPI requests with a given MPI message tag.
     * The handler call back function will be called when the MPI
     * request has been finished where tag and rank can be used to
     * keep track of the data that a finished MPI request belongs to
     * (e.g., for clean-up of allocated heap data).
     */
    void submitRequests(
    	MPI_Request *requests,
		int nRequests,
		int tag,
		int remoteRank,
        std::function<void(exahype::solvers::Solver*, int , int)> handleRequest,
        RequestType type,
		exahype::solvers::Solver *solver,
		bool block=false);
    void progressRequests();

    void createMPICommunicator();
    MPI_Comm getMPICommunicator();
    MPI_Comm getMPICommunicatorMapped();

    /*
     * Given the current load situation and global knowledge of the load
     * of the other ranks, this method selects a victim rank, i.e.,
     * a rank to which local work should be offloaded in order to
     * improve the load balance.
     */
    bool selectVictimRank(int& victim);

    static StealingManager& getInstance();
    virtual ~StealingManager();
};

#endif
