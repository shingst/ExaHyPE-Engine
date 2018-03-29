#include "mpi.h"
#include "stdlib.h"
#include "time.h"

#include <iostream>

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);

  int rank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  int flag = 0;
  int numberCount = 0;
  // master
  if (rank == 0) {
    const int MAX_NUMBERS = 100;
    int numbers[MAX_NUMBERS];
    // Pick a random amount of integers to send to process one
    srand(time(NULL));
    numberCount = (rand() / (float)RAND_MAX) * MAX_NUMBERS;

    int numberOfProcessors = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcessors );
    #if defined(BlockingSend)
    for ( int workerRank = 1; workerRank < numberOfProcessors; workerRank++ ) {
      MPI_Send(numbers, numberCount, MPI_INTEGER, workerRank, 0, MPI_COMM_WORLD);
    }
    #else
    // post all sends
    MPI_Request* sendRequests = new MPI_Request[numberOfProcessors-1];
    for ( int workerRank = 1; workerRank < numberOfProcessors; workerRank++ ) {
      MPI_Isend(numbers, numberCount, MPI_INTEGER, workerRank, 0, MPI_COMM_WORLD, &sendRequests[workerRank-1]);
    }

    // wait till all sends are finished; use MPI_Test loop in order to allow timeout check
    bool allMessagesSent = false;
    while (!allMessagesSent) {
      allMessagesSent = true;
      for ( int workerRank = 1; workerRank < numberOfProcessors; workerRank++ ) {
        MPI_Test(&sendRequests[workerRank-1],&flag,MPI_STATUS_IGNORE);
        if (flag) {
          std::cout << "0 sent " << numberCount << " numbers to " << workerRank << std::endl;
        }
        allMessagesSent &= flag;
      }
    }
    delete[] sendRequests;
    #endif
  }
  // workers
  else if (rank > 0) {
    // Probe for an incoming message from process zero

    MPI_Status status;
    #if defined(BlockingProbe)
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    #else
    flag = 0;
    while (!flag) {
      MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);
    }
    #endif

    // When probe returns, the status object has the size and other
    // attributes of the incoming message. Get the message size
    MPI_Get_count(&status, MPI_INT, &numberCount);

    // Allocate a buffer to hold the incoming numbers
    int* numberBuffer = new int[numberCount];

    #if defined(BlockingReceive)
    MPI_Recv(numberBuffer, numberCount, MPI_INT, 0 /*fromRank*/, 0/*tag*/, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #else
    // post the receive
    flag = 0;
    MPI_Request receiveRequest;
    MPI_Irecv(numberBuffer, numberCount, MPI_INT, 0, 0, MPI_COMM_WORLD, &receiveRequest);
    // Now wait
    while (!flag) {
      MPI_Test(&receiveRequest,&flag,MPI_STATUS_IGNORE);
    }
    #endif

    std::cout << rank << " received " << numberCount << " numbers from " << 0 << std::endl;
    delete[] numberBuffer;
  }

  MPI_Finalize();
}
