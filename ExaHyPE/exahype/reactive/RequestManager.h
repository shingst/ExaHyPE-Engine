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
#define MAX_THREADS 48
#else
#define MAX_THREADS 1
#endif

namespace exahype {
  namespace reactive {
    class RequestManager;
    /**
     * Different MPI request types for prioritizing message requests.
     */
    enum class RequestType {
      send = 0,
      sendBack = 1,
      receive = 2,
      receiveBack = 3,
      sendOutcome = 4,
      receiveOutcome = 5
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
 * The request manager manages the created asynchronous non-blocking MPI requests,
 * arising for instance when a task is offloaded.
 * It allows to progress any outstanding requests and invoke call back methods
 * whenever outstanding requests have been completed.
 */
class exahype::reactive::RequestManager {
  private:
    /**
     * The logging device.
     */
    static tarch::logging::Log _log;

    int _threadId;

    /**
     * Semaphore to ensure that only one thread at a time
     * makes progress.
     */
    tarch::multicore::BooleanSemaphore _progressSemaphore;

    /**
     *  MPI_Request handles are mapped to an internal id.
     */
    std::atomic<int> _nextRequestId;

    /**
     *  MPI_Request groups (e.g. all requests belonging to a task send) are mapped to
     *  an internal group id. A request group is a logically belonging together collection of
     *  requests.
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

    // some counters for debugging
    std::atomic<int> *_postedSendsPerRank;
    std::atomic<int> *_postedReceivesPerRank;
    std::atomic<int> *_postedSendBacksPerRank;
    std::atomic<int> *_postedReceiveBacksPerRank;
    std::atomic<int> *_postedSendOutcomesPerRank;
    std::atomic<int> *_postedReceiveOutcomesPerRank;

    RequestManager(int threadId);

    /**
     * This method makes progress on all current requests of the given request type.
     */
    bool progressRequestsOfType(RequestType type);

    /**
     * Maps a request type to an integer that defines message queue.
     */
    static int requestTypeToMsgQueueIdx(RequestType requestType);

    /**
     * This method pop's #limit current requests of a given type from the request queue and
     * inserts them into the active array of MPI requests on which we can make progress.
     */
    void createRequestArray(
        RequestType type,
        std::vector<MPI_Request> &requests,
        std::unordered_map<int, int> &vecIdToReqId,
        int limit = std::numeric_limits<int>::max());

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


    inline int getNextRequestId() {
      // Todo: Deal with overflow
      return _nextRequestId++;
    }

    inline int getNextGroupId() {
      // Todo: Deal with overflow
      return _nextGroupId++;
    }

    static RequestManager* _static_managers[MAX_THREADS];

  public:
    static constexpr int MULTIPLE_SOURCES = -1;

    int getNumberOfOutstandingRequests(RequestType type);

    void printPostedRequests();
    void resetPostedRequests();

    /**
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
    void progressAnyRequests();
    bool progressReceiveBackRequests();
    bool hasOutstandingRequestOfType(RequestType requestType);

#if defined (DirtyCleanUp)
    void cancelOutstandingRequests();
#endif


#ifdef OffloadingUseProgressTask
    void resetHasNotifiedSendCompleted();
    void notifySendCompleted(int rank);
    void receiveCompleted(int rank, int rail=-1);
    void notifyAllVictimsSendCompletedIfNotNotified();
#endif

    static RequestManager& getInstance();
    virtual ~RequestManager();
};

#endif
