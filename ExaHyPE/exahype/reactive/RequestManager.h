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

#if !defined(_EXAHYPE_REQUESTMANAGER_H_)  && defined(Parallel) && defined(SharedTBB)
#define _EXAHYPE_REQUESTMANAGER_H_

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Jobs.h"

#include "tbb/concurrent_hash_map.h"
#include "tbb/concurrent_queue.h"

#include <mpi.h>
#include <atomic>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

#ifndef MPI_CHECK
#ifndef Asserts
#define MPI_CHECK(func, x) do { \
  ierr = (x); \
  if (ierr != MPI_SUCCESS) { \
    logError(#func, "Runtime error:"<<#x<<" returned "<<ierr<<" at " << __FILE__<< ":"<< __LINE__); \
  } \
} while (0)
#else
#define MPI_CHECK(func, x) do { \
  ierr = (x); \
  } while (0)
#endif
#endif

#if defined(UseMPIThreadSplit)
#ifndef MAX_THREADS
#define MAX_THREADS 48
#endif
#else
#ifndef MAX_THREADS
#define MAX_THREADS 1
#endif
#endif

namespace exahype {
  namespace reactive {
    class RequestManager;
    /**
     * Different MPI request types for prioritizing message requests.
     */
    enum class RequestType {
      Send = 0,
      SendBack = 1,
      Receive = 2,
      ReceiveBack = 3,
      SendOutcome = 4,
      ReceiveOutcome = 5
    };
  }
}

namespace exahype {
  namespace solvers {
    //forward declaration
    class Solver;
  }
}

/**
 * The request manager manages the asynchronous non-blocking MPI requests,
 * arising for instance when a task is offloaded or when a task outcome is shared.
 * It allows to progress any outstanding requests and invoke call back methods
 * whenever outstanding requests have been completed.
 */
class exahype::reactive::RequestManager {
  private:
    /**
     * The logging device.
     */
    static tarch::logging::Log _log;

    /**
     * Thread id this request manager belongs to.
     */
    int _threadId;

    /**
     * Semaphore to ensure that only one thread at a time
     * makes progress (in case UseMPIThreadSplit is not active).
     */
    tarch::multicore::BooleanSemaphore _progressSemaphore;

    /**
     *  Internal request id for next request.
     *  MPI_Request handles are mapped to an internal integer request id.
     */
    std::atomic<int> _nextRequestId;

    /**
     *  Internal group id for next request group.
     *  Groups of MPI_Requests (e.g. all requests belonging to a task send) are mapped to
     *  an internal integer group id. A request group is a collection of
     *  requests which logically belongs together.
     */
    std::atomic<int> _nextGroupId;

    /**
     * Request queues for each request type.
     */
    tbb::concurrent_queue<int> _outstandingRequests[6];

    /**
     * Maps internal integer request id to MPI_Request handle.
     */
    tbb::concurrent_hash_map<int, MPI_Request> _reqIdToReqHandle[6];

    /**
     * Maps internal request id to internal group id.
     */
    tbb::concurrent_hash_map<int, int>         _reqIdToGroup[6];

    /**
     * Map that tracks how many outstanding requests a request group still has.
     */
    tbb::concurrent_hash_map<int, int>         _outstandingReqsForGroup[6];

    /**
     * Maps a group id to the remote rank which its requests belong to.
     */
    tbb::concurrent_hash_map<int, int>         _groupIdToRank[6];

    /**
     * Maps a group id to the MPI tag which its requests belong to.
     */
    tbb::concurrent_hash_map<int, int>         _groupIdToTag[6];

    /**
     * Maps a group id to the handler function which is invoked when
     * the group's requests have been completed.
     */
    tbb::concurrent_hash_map<int, std::function<void(exahype::solvers::Solver*, int , int)>> _handlers[6];

    /**
     * Maps a group id to a pointer to the solver to which its requests
     * logically belong to.
     */
    tbb::concurrent_hash_map<int, exahype::solvers::Solver*> _solvers[6];

    /**
     * The vector of MPI_Request handles which the manager currently
     * makes progress on.
     *
     * @Note: Once this vector becomes empty, the code
     * tries to grab new outstanding requests from the
     * request queue and - if possible - creates a new
     * vector of outstanding requests.
     */
    std::vector<MPI_Request> _activeRequests[6];

    /**
     * This array maps the elements in _currentOutstandingRequests
     * back to the internal request ids.
     */
    std::unordered_map<int, int> _internalIdsOfActiveRequests[6];

    // some counters for debugging (RAII not easily possible since there is no copy constructor for atomics)
    std::atomic<int> *_postedSendsPerRank;
    std::atomic<int> *_postedReceivesPerRank;
    std::atomic<int> *_postedSendBacksPerRank;
    std::atomic<int> *_postedReceiveBacksPerRank;
    std::atomic<int> *_postedSendOutcomesPerRank;
    std::atomic<int> *_postedReceiveOutcomesPerRank;

    /**
     * @param Id of the thread this manager belongs to.
     */
    RequestManager(int threadId);
    virtual ~RequestManager();

    RequestManager(const RequestManager& other) = delete;
    RequestManager& operator=(const RequestManager& other) = delete;

    RequestManager(RequestManager&& other) = delete;
    RequestManager& operator=(RequestManager&& other) = delete;

    /**
     * This method makes progress on all current requests of the given request type.
     */
    bool progressRequestsOfType(RequestType type);

    void processFinishedRequests(int mapId, int numCompleted, std::vector<int>& indicesOfCompletedRequests);

    /**
     * Maps a request type to an integer that defines message queue.
     */
    static int requestTypeToMsgQueueIdx(RequestType requestType);

    /**
     * This method pops #limit current requests of a given type from the request queue and
     * inserts them into the active array of MPI requests on which we can make progress.
     */
    void createAndFillRequestArray(
        RequestType type,
        std::vector<MPI_Request> &requests,
        std::unordered_map<int, int> &vecIdToReqId,
        int limit = std::numeric_limits<int>::max());

    /**
     * @return Next internal request id.
     */
    int getNextRequestId() {
      // overflow is ignored for now, as wrap around would not cause any problems
      return _nextRequestId++;
    }

    /**
     * @return Next internal request group id.
     */
    int getNextGroupId() {
      // overflow is ignored for now, as wrap around would not cause any problems
      return _nextGroupId++;
    }

    /**
     * There may be multiple request managers with UseMPIThreadSplit
     */
    static RequestManager* StaticManagers[MAX_THREADS];

    /**
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

    /**
     * Special flag indicating that a request group contains requests for several MPI ranks (receivers, senders).
     */
    static constexpr int MULTIPLE_RANKS = -1;

    /**
     * Singleton (per thread).
     */
    static RequestManager& getInstance();

    /**
     * For debugging.
     */
    int getNumberOfOutstandingRequestsOfType(RequestType type) const;

    /**
     * For debugging.
     */
    void printPostedRequests() const;

    /**
     * For debugging.
     */
    void resetPostedRequests();

    /**
     * Submit a group of MPI requests with a given MPI message tag.
     * The handler call back function will be called when the MPI
     * request has been finished where tag and rank can be used to
     * keep track of the data that a finished MPI request belongs to
     * (e.g., for clean-up of allocated heap data).
     *
     * @param requests Pointer to array of nRequests non-blocking requests
     * @param nRequests Size of requests array
     * @param tag MPI tag of all requests
     * @param remoteRank MPI rank from/to which data is received/sent. Can be MULTIPLE_RANKS (@see MULTIPLE_RANKS).
     * @param reqHandler Function to be invoked when request group has finished
     * @param type Type of request (@see RequestType)
     * @param solver Pointer to solver from which MPI requests were sent
     * @param block If true, the call blocks until all requests have finished (MPI_Wait)
     */
    void submitRequests(
        MPI_Request *requests,
        int nRequests,
        int tag,
        int remoteRank,
        std::function<void(exahype::solvers::Solver*, int , int)> reqHandler,
        RequestType type,
        exahype::solvers::Solver *solver,
        bool block=false);

    /**
     * Makes progress on outstanding request groups and
     * calls completion handler if all requests in a group can be completed.
     */
    void progressRequests();

    /**
     * Progresses receive back requests only.
     * @return True if some requests could be completed.
     */
    bool progressReceiveBackRequests();

    /**
     * @param requestType RequestType to check for outstanding requests.
     * @return True if there are outstanding requests of the given type.
     */
    bool hasOutstandingRequestOfType(RequestType requestType) const;

#if defined (DirtyCleanUp)
    void cancelOutstandingRequests();
#endif

};

#endif
