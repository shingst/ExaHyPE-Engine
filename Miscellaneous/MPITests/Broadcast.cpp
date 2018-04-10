#include "mpi.h"
#include "stdlib.h"
#include "time.h"

#include <iostream>
#include <chrono>

#if !defined(MessageSize)
#define MessageSize 100
#endif

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);

  int rank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  const int numberCount = MessageSize;
  int numbers[numberCount];

  int flag = 0;
  #if defined(Collectives)
  int numberOfProcessors = 0;
  MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcessors );

  auto t1 = std::chrono::high_resolution_clock::now();
  MPI_Bcast(numbers, numberCount, MPI_INTEGER, 0, MPI_COMM_WORLD);
  auto t2 = std::chrono::high_resolution_clock::now();

  if (rank == 0) {
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << rank << " broadcasted " << numberCount << " numbers to all " << numberOfProcessors <<" receivers in " << duration << " microseconds " <<  std::endl;
  }
  #else
  if (rank == 0) { // master
    int numberOfProcessors = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcessors );

    auto t1 = std::chrono::high_resolution_clock::now();
    #if defined(BlockingSend)
    for ( int workerRank = 1; workerRank < numberOfProcessors; workerRank++ ) {
      MPI_Send(numbers, MAX_NUMBERS, MPI_INTEGER, workerRank, 0, MPI_COMM_WORLD);
    }
    #else
    // post all sends
    MPI_Request* sendRequests = new MPI_Request[numberOfProcessors-1];
    for ( int workerRank = 1; workerRank < numberOfProcessors; workerRank++ ) {
      MPI_Isend(numbers, numberCount, MPI_INTEGER, workerRank, 0, MPI_COMM_WORLD, &sendRequests[workerRank-1]);
    }
    // wait till all sends are finished; use MPI_Test loop in order to allow timeout check
    bool allMessagesSent = false;
    #if defined(WaitOnSends)
    MPI_Waitall(numberOfProcessors-1,sendRequests,MPI_STATUSES_IGNORE);
    #else
    while (!allMessagesSent) {
      allMessagesSent = true;
      for ( int workerRank = 1; workerRank < numberOfProcessors; workerRank++ ) {
        MPI_Test(&sendRequests[workerRank-1],&flag,MPI_STATUS_IGNORE);
        allMessagesSent &= flag;
      }
    }
    delete[] sendRequests;
    #endif
    #endif

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << rank << " broadcasted " << numberCount << " numbers to all " << numberOfProcessors <<" receivers in " << duration << " microseconds " <<  std::endl;
  }
  // workers
  else if (rank > 0) {
    // Allocate a buffer to hold the incoming numbers
    int numberBuffer[numberCount];

    #if defined(BlockingReceive)
    MPI_Recv(numberBuffer, numberCount, MPI_INT, 0 /*fromRank*/, 0/*tag*/, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #else
    // post the receive
    flag = 0;
    MPI_Request receiveRequest;
    MPI_Irecv(numberBuffer, numberCount, MPI_INT, 0, 0, MPI_COMM_WORLD, &receiveRequest);
    // Now wait
    #if defined(WaitOnReceive)
    MPI_Wait(&receiveRequest,MPI_STATUS_IGNORE);
    #else
    while (!flag) {
      MPI_Test(&receiveRequest,&flag,MPI_STATUS_IGNORE);
    }
    #endif
    #endif
  }
  #endif

  MPI_Finalize();
}
