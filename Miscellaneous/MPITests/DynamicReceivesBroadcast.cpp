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

  int flag = 0;
  int numberCount = 0;
  // master
  if (rank == 0) {
    // Pick a random amount of integers to send to process one
    // const int MAX_NUMBERS = 100;
    // srand(time(NULL));
    // numberCount = (rand() / (float)RAND_MAX) * MAX_NUMBERS;
    numberCount = MessageSize;
    int numbers[numberCount];

    int numberOfProcessors = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcessors );

    auto t1 = std::chrono::high_resolution_clock::now();
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
        allMessagesSent &= flag;
      }
    }
    delete[] sendRequests;
    #endif

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << rank << " sent " << numberCount << " numbers to all " << numberOfProcessors <<" receivers in " << duration << " microseconds " <<  std::endl;
  }
  // workers
  else if (rank > 0) {
    // Probe for an incoming message from process zero

    flag = 0;
    MPI_Status status;
    #if defined(TobiasReceives)
    while (!flag) {
      MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);

      if (flag) {
        // When probe returns, the status object has the size and other
        // attributes of the incoming message. Get the message size
        MPI_Get_count(&status, MPI_INT, &numberCount);

        // Allocate a buffer to hold the incoming numbers
        int* numberBuffer = new int[numberCount];

        bool receiveComplete = false;
        int attempts = 0;
        while (!receiveComplete) {
          // post receive
          MPI_Request receiveRequest;
          MPI_Irecv(numberBuffer, numberCount, MPI_INT, 0, 0, MPI_COMM_WORLD, &receiveRequest);
          attempts++;

          // test receive
          int flag2 = 0;
          MPI_Test(&receiveRequest,&flag2,MPI_STATUS_IGNORE);
          receiveComplete = flag2;
        }

        delete[] numberBuffer;
        std::cout << rank << " received " << numberCount << " numbers from " << 0 << " in " << attempts << " attempts "<< std::endl;
      }
    }
    #else
    #if defined(BlockingProbe)
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    #else
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
    delete[] numberBuffer;
    #endif
    std::cout << rank << " received " << numberCount << " numbers from " << 0 << std::endl;
  }

  MPI_Finalize();
}
