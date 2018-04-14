#include "mpi.h"
#include "stdlib.h"
#include "time.h"

#include <iostream>
#include <sstream>

#include <map>
#include <list>
#include <vector>
#include <tuple>

#include <cmath>
#include <algorithm>

#include <chrono>
#include <iomanip>

#if defined(ReceiveDanglingMessagesBlocking) and !defined(ReceiveDanglingMessages)
#define ReceiveDanglingMessages
#endif

// kindly copied from: https://stackoverflow.com/a/365068
inline void pow2RoundUp (int& x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x +=1;
}

void parseArguments(int argc, char** argv, int& dimensions, int& maximumMessageSize, int& numberOfTests) {

  if (argc < 4) { // program name is first argument
    std::cerr << "Error: Please call the program the following way: ./CartesianExchange <dimensions> <maximumMessageSize> <numberOfTests>" <<
        " where dimensions, maximumMessageSize, and numberOfTests are all positive integers." << std::endl;
    std::terminate();
  }

  std::istringstream ss1(argv[1]);
  if (!(ss1 >> dimensions) || dimensions<=0 ) {
    std::cerr << "Error: First argument 'dimensions' must be positive integer but is: '" << argv[1] << '\n';
    std::terminate();
  }

  std::istringstream ss2(argv[2]);
  if (!(ss2 >> maximumMessageSize) || maximumMessageSize<=0 ) {
    std::cerr << "Error: Second argument 'maximumMessageSize' must be positive integer but is: '" << argv[2] << '\n';
    std::terminate();
  }
  pow2RoundUp(maximumMessageSize);

  std::istringstream ss3(argv[3]);
  if (!(ss3 >> numberOfTests) || numberOfTests<=0 ) {
    std::cerr << "Error: Third argument 'numberOfTest' must be positive integer but is: '" << argv[3] << '\n';
    std::terminate();
  }
}

typedef std::tuple<int,int,double,double,double> record; // messageSize,numberOfMessages,avg,min,max

void printRecords(std::vector<record>& records,std::string name) {
  std::cout << std::endl;
  std::cout << name << ":" << std::endl;
  std::cout << std::endl;
  std::cout << "| message size | #messages    | avg time/us  | min time/us  | max time/us  |" << std::endl;
  std::cout << "|--------------|--------------|--------------|--------------|--------------|" << std::endl;
  for (record r :  records) {
    std::cout <<  "| " << std::setw(12) << std::get<0>(r) << " ";
    std::cout <<  "| " << std::setw(12) << std::get<1>(r) << " ";

    std::stringstream t;
    t << std::fixed << std::setprecision(1) << std::get<2>(r);
    std::cout << "| " << std::setw(12) << t.str() << " ";

    t.str("");
    t << std::fixed << std::setprecision(0) << std::get<3>(r);
    std::cout << "| " << std::setw(12) << t.str() << " ";

    t.str("");
    t << std::fixed << std::setprecision(0) << std::get<4>(r);
    std::cout << "| " << std::setw(12) << t.str() << " ";

    std::cout << "|" << std::endl;
  }
}

void printIntro(const int dimensions, const int maximumMessageSize, const int numberOfTests) {
  std::cout << std::endl;

  std::cout << "Start experiment with parameters: "       << std::endl
      << std::endl
      << "dimensions         = " << dimensions            << std::endl
      << "maximumMessageSize = " << maximumMessageSize    << " (rounded to next power of 2)" << std::endl
      << "numberOfTests      = " << numberOfTests         << std::endl;

  
  std::cout << std::endl << "Compiler options: " << std::endl << std::endl;
  #if defined(UseVector)
  std::cout
      << "UseVector                       = yes" << std::endl;
  #endif
  #if defined(BlockPerRank)
  std::cout
      << "BlockPerRank                    = yes" << std::endl;
  #endif
  #if defined(ReceiveDanglingMessages)
  std::cout
      << "ReceiveDanglingMessages         = yes" << std::endl;
  #endif
  #if defined(ReceiveDanglingMessagesBlocking)
  std::cout
      << "ReceiveDanglingMessagesBlocking = yes" << std::endl;
  #endif
  #if defined(DynamicReceives)
  std::cout
      << "DynamicReceives                 = yes" << std::endl;
  #endif
}

#if defined(UseVector)
typedef std::vector<MPI_Request*> RequestContainer;
#else
typedef std::list<MPI_Request*> RequestContainer;
#endif

void receiveDanglingMessages(MPI_Comm& cartesianComm, const int myRank, double* receiveBuffer, std::map<int,RequestContainer>& receiveRequests) {
  // assumption: keys of receiveRequests hold the neighbour ranks
  int flag = 0;
  for (auto rankIt = receiveRequests.begin(); rankIt != receiveRequests.end(); rankIt++) {
    MPI_Status status;
    MPI_Iprobe(rankIt->first,0,cartesianComm,&flag,&status);
    if (flag) {
       int messageSize;
       MPI_Get_count(&status,MPI_DOUBLE,&messageSize);
       MPI_Request* receiveRequest = new MPI_Request();
       #if defined(ReceiveDanglingMessagesBlocking)
       MPI_Recv(receiveBuffer, messageSize, MPI_DOUBLE, rankIt->first, 0, cartesianComm,MPI_STATUS_IGNORE);
       #else
       MPI_Irecv(receiveBuffer, messageSize, MPI_DOUBLE, rankIt->first, 0, cartesianComm, receiveRequest);
       #endif
       receiveRequests[rankIt->first].push_back(receiveRequest); // receive request not necessary when blocking but we want to fill the receiveRequests list
    }
  }
}

#if defined(DynamicReceives) and  defined(ReceiveDanglingMessages)
#error DynamicReceives and ReceiveDanglingMessages must not be defined together!
#endif

#if defined(TestSendAndReceiveTogether) and defined(DynamicReceives)
#error TestSendAndReceiveTogether and DynamicReceives must not be defined together!
#endif

#if defined(TestSendAndReceiveTogether) and defined(ReceiveDanglingMessages) 
#error TestSendAndReceiveTogether and ReceiveDanglingMessages must not be defined together!
#endif

// Arguments:        dimensions, maximumMessageSize, numberOfTests
// Compiler Options: -DUseVector
int main(int argc, char** argv) {
  ////////////////////
  // Initialisation //
  ////////////////////

  MPI_Init(&argc,&argv);

  int dimensions          = -1;
  int maximumMessageSize  = -1;
  int numberOfTests       = -1;
  parseArguments(argc,argv,dimensions,maximumMessageSize,numberOfTests);

  // compute the dimensions
  int numberOfRanks;
  MPI_Comm_size( MPI_COMM_WORLD, &numberOfRanks);

  int ranksPerDimension[dimensions];
  ranksPerDimension[0] = std::round( std::pow( numberOfRanks, 1.0/dimensions) );
  std::fill_n(&ranksPerDimension[1],dimensions-1,ranksPerDimension[0]);
  int product = 1;
  for (int d=0; d<dimensions; d++) {
      product *= ranksPerDimension[d];
  }
  if ( product != numberOfRanks ) {
    std::cerr << "ERROR: Specified number of ranks cannot be distributed evenly among "<<dimensions<<" dimensions (number of ranks specified: "<<numberOfRanks,
    std::terminate();
  }

  // create the Cartesian topology
  MPI_Comm cartesianComm;
  int periodicBoundaries[dimensions];
  std::fill_n(periodicBoundaries, dimensions, 0);  

  MPI_Cart_create(
      MPI_COMM_WORLD, dimensions, ranksPerDimension, periodicBoundaries,
      1 /* reorder, i.e. redistribute the ranks */,
      &cartesianComm);

  // obtain my rank in Cartesian topology
  int myRank = -1;
  MPI_Comm_rank(cartesianComm, &myRank);
  if ( myRank==0 ) {
    printIntro(dimensions,maximumMessageSize,numberOfTests);
  }

  // create a container for neighbouring ranks
  #ifdef UseVector
  std::map<int,RequestContainer> sendRequests;
  std::map<int,RequestContainer> receiveRequests;
  #else
  std::map<int,RequestContainer> sendRequests;
  std::map<int,RequestContainer> receiveRequests;
  #endif

  /////////////////////////
  // Run the experiments //
  /////////////////////////

  // statistics
  std::vector<record> tPosting; // statistics for program phases
  std::vector<record> tClear;
  std::vector<record> tWait;
  tPosting.reserve(100);
  tClear.reserve(100);
  tWait.reserve(100);

  int l = 0;
  int messageSize      = std::round( std::pow( 2, l) );
  int numberOfMessages = maximumMessageSize/messageSize;
  while ( messageSize <= maximumMessageSize ) { // loop over message sizes
    double* sendBuffer    = new double[messageSize];
    double* receiveBuffer = new double[messageSize];

    tPosting.push_back( std::make_tuple(messageSize,numberOfMessages,0.0,99999999.0,0.0) );
    tClear.push_back( std::make_tuple(messageSize,numberOfMessages,0.0,99999999.0,0.0) );
    tWait.push_back( std::make_tuple(messageSize,numberOfMessages,0.0,99999999.0,0.0) );
    for (int test = 0; test < numberOfTests; ++test) {
      MPI_Barrier(cartesianComm);

      // send out the packets
      auto tPostingBegin = std::chrono::high_resolution_clock::now();
      for (int direction=0; direction<dimensions; direction++) {
        for (int orientation=0; orientation<2; orientation++) {
          const int displacement = -1 + 2*orientation;

          int sourceNeighbour      = -1; // this one wants to send to me       (blocking comm.)
          int destinationNeighbour = -1; // this one receives messages from me (blocking comm.)
          MPI_Cart_shift(cartesianComm, direction ,displacement,&sourceNeighbour,&destinationNeighbour);

          // We do nonblocking communication, so we do not care about the sourceNeighbour
          // We further always send/receive from/to the same buffer location.
          if ( destinationNeighbour!=MPI_PROC_NULL) {
            sendRequests.insert( std::pair<int,RequestContainer>( destinationNeighbour, RequestContainer(0) ) );
            receiveRequests.insert( std::pair<int,RequestContainer>( destinationNeighbour, RequestContainer(0) ) );
            #if defined(UseVector)
            sendRequests[destinationNeighbour].reserve(numberOfMessages);
            receiveRequests[destinationNeighbour].reserve(numberOfMessages);
            #endif

            for (int m = 0; m < numberOfMessages; ++m) {
              #if !defined(ReceiveDanglingMessages) && !defined(DynamicReceives)
              MPI_Request* receiveRequest = new MPI_Request();
              MPI_Irecv(receiveBuffer, messageSize, MPI_DOUBLE, destinationNeighbour, 0, cartesianComm, receiveRequest);
              receiveRequests[destinationNeighbour].push_back(receiveRequest);
              #endif

              MPI_Request* sendRequest = new MPI_Request();
              MPI_Isend(sendBuffer, messageSize, MPI_DOUBLE, destinationNeighbour, 0, cartesianComm, sendRequest);
              sendRequests[destinationNeighbour].push_back(sendRequest);
            }
          }
        }
      }
      auto tPostingEnd      = std::chrono::high_resolution_clock::now();
      double tPostingOfTest = static_cast<double>( std::chrono::duration_cast<std::chrono::microseconds>( tPostingEnd - tPostingBegin ).count() );
      std::get<2>(tPosting.back()) += tPostingOfTest;
      std::get<3>(tPosting.back())  = std::min<double>( std::get<3>(tPosting.back()), tPostingOfTest );
      std::get<4>(tPosting.back())  = std::max<double>( std::get<4>(tPosting.back()), tPostingOfTest );

      // wait for completion
      auto tWaitBegin = std::chrono::high_resolution_clock::now();
      int flag      = 0;
      #if defined(BlockPerRank)
      for (auto rankIt = sendRequests.begin(); rankIt != sendRequests.end(); rankIt++) {
        bool complete = false;
        while (!complete) {
          complete = receiveRequests[rankIt->first].size() == sendRequests[rankIt->first].size();
          int index = 0;
          for (auto requestIt = sendRequests[rankIt->first].begin(); requestIt != sendRequests[rankIt->first].end(); requestIt++) {
            #if defined(TestSendAndReceiveTogether)
            MPI_Request* receiveRequest = receiveRequests[rankIt->first][index];
            MPI_Test(receiveRequest,&flag,MPI_STATUS_IGNORE);
            complete &= flag;
            index++;
            #else
            MPI_Test(*requestIt,&flag,MPI_STATUS_IGNORE);
            complete &= flag;
            #endif

            #if defined(DynamicReceives)
            MPI_Status status;
            MPI_Iprobe(rankIt->first,0,cartesianComm,&flag,&status);
            if (flag) {
               int messageSize;
               MPI_Get_count(&status,MPI_DOUBLE,&messageSize);
               MPI_Recv(receiveBuffer, messageSize, MPI_DOUBLE, rankIt->first, 0, cartesianComm,MPI_STATUS_IGNORE);
               MPI_Request* receiveRequest = new MPI_Request(); // receive request not necessary when blocking but we want to fill the receiveRequests list
               receiveRequests[rankIt->first].push_back(receiveRequest);
            }
            #endif
          }

          #if defined(ReceiveDanglingMessages)
          receiveDanglingMessages(cartesianComm,myRank,receiveBuffer,receiveRequests);
          #endif

          #if !defined(DynamicReceives) and !defined(ReceiveDanglingMessagesBlocking) and !defined(TestSendAndReceiveTogether)
          for (auto requestIt = receiveRequests[rankIt->first].begin(); requestIt != receiveRequests[rankIt->first].end(); requestIt++) {
            MPI_Test(*requestIt,&flag,MPI_STATUS_IGNORE);
            complete &= flag;
          }
          #endif
        }
      }
      #else
      bool complete = false;
      while (!complete) {
        complete = true;
        for (auto rankIt = sendRequests.begin(); rankIt != sendRequests.end(); rankIt++) {
          complete &= receiveRequests[rankIt->first].size() == sendRequests[rankIt->first].size();

          int index = 0;
          for (auto requestIt = sendRequests[rankIt->first].begin(); requestIt != sendRequests[rankIt->first].end(); requestIt++) {
            #if defined(TestSendAndReceiveTogether)
            MPI_Request* receiveRequest = receiveRequests[rankIt->first][index];
            MPI_Test(receiveRequest,&flag,MPI_STATUS_IGNORE);
            complete &= flag;
            index++;
            #else
            MPI_Test(*requestIt,&flag,MPI_STATUS_IGNORE);
            complete &= flag;
            #endif
            
            #if defined(DynamicReceives)
            MPI_Status status;
            MPI_Iprobe(rankIt->first,0,cartesianComm,&flag,&status);
            if (flag) {
               int messageSize;
               MPI_Get_count(&status,MPI_DOUBLE,&messageSize);
               MPI_Recv(receiveBuffer, messageSize, MPI_DOUBLE, rankIt->first, 0, cartesianComm,MPI_STATUS_IGNORE);
               MPI_Request* receiveRequest = new MPI_Request(); // receive request not necessary when blocking but we want to fill the receiveRequests list
               receiveRequests[rankIt->first].push_back(receiveRequest);
            }
            #endif
          }

          #if defined(ReceiveDanglingMessages)
          receiveDanglingMessages(cartesianComm,myRank,receiveBuffer,receiveRequests);
          #endif

          #if !defined(DynamicReceives) and !defined(ReceiveDanglingMessagesBlocking) and !defined(TestSendAndReceiveTogether)
          for (auto requestIt = receiveRequests[rankIt->first].begin(); requestIt != receiveRequests[rankIt->first].end(); requestIt++) {
            MPI_Test(*requestIt,&flag,MPI_STATUS_IGNORE);
            complete &= flag;
          }
          #endif
        }
      }
      #endif
      auto tWaitEnd = std::chrono::high_resolution_clock::now();
      double tWaitOfTest = static_cast<double>( std::chrono::duration_cast<std::chrono::microseconds>( tWaitEnd - tWaitBegin ).count() );
      std::get<2>(tWait.back()) += tWaitOfTest;
      std::get<3>(tWait.back())  = std::min<double>( std::get<3>(tWait.back()), tWaitOfTest );
      std::get<4>(tWait.back())  = std::max<double>( std::get<4>(tWait.back()), tWaitOfTest );

      // free requests and clear the containers
      auto tClearBegin = std::chrono::high_resolution_clock::now();

      for (auto rankIt = sendRequests.begin(); rankIt != sendRequests.end(); rankIt++) {
        for (auto requestIt = sendRequests[rankIt->first].begin(); requestIt != sendRequests[rankIt->first].end(); requestIt++) {
          delete *requestIt;
        }
        sendRequests[rankIt->first].clear();
        for (auto requestIt = receiveRequests[rankIt->first].begin(); requestIt != receiveRequests[rankIt->first].end(); requestIt++) {
          delete *requestIt;
        }
        receiveRequests[rankIt->first].clear();
      }
      sendRequests.clear();
      receiveRequests.clear();

      auto tClearEnd = std::chrono::high_resolution_clock::now();
      double tClearOfTest = static_cast<double>( std::chrono::duration_cast<std::chrono::microseconds>( tClearEnd - tClearBegin ).count() );
      std::get<2>(tClear.back()) += tClearOfTest;
      std::get<3>(tClear.back())  = std::min<double>( std::get<3>(tClear.back()), tClearOfTest );
      std::get<4>(tClear.back())  = std::max<double>( std::get<4>(tClear.back()), tClearOfTest );
    }

    // compute the averages
    std::get<2>(tPosting.back()) /= numberOfTests;
    std::get<2>(tWait.back())     /= numberOfTests;
    std::get<2>(tClear.back())    /= numberOfTests;

    // clean up
    delete[] sendBuffer;
    delete[] receiveBuffer;

    // go into next iteration
    l++;
    messageSize      = std::round( std::pow( 2, l) );
    numberOfMessages = maximumMessageSize/messageSize;
  }

  // outro
  bool isCloseToCenter = true;
  int coordinates[dimensions];
  MPI_Cart_coords(cartesianComm,myRank,dimensions,coordinates);
  for (int d = 0; d < dimensions; ++d) {
    isCloseToCenter &= (coordinates[d] == ranksPerDimension[0]/2);
  }
  if (isCloseToCenter) {
    std::cout << std::endl;

    printRecords(tPosting,"Posting Messages");
    printRecords(tWait,"Waiting on Completion");
    printRecords(tClear,"Clearing Buffers");

    std::cout <<
        std::endl <<
        "info: measurements obtained for rank "<< myRank <<
        " located at x=(";
    for (int direction=0; direction<dimensions; direction++) { 
       std::cout <<ranksPerDimension[0]/2;
       if (direction<dimensions-1) {
         std::cout << ",";
       } else {
         std::cout << ")" << " in domain=[0,"<<ranksPerDimension[0]-1<<"]^"<<dimensions<<""<< std::endl;
       }
    }
 
    std::cout <<
        "info: rank "<<myRank<<" is neighbour of the ranks [";
    for (int direction=0; direction<dimensions; direction++) {
      for (int orientation=0; orientation<2; orientation++) {
        const int displacement = -1 + 2*orientation;

        int sourceNeighbour      = -1; // this one wants to send to me       (blocking comm.)
        int destinationNeighbour = -1; // this one receives messages from me (blocking comm.)
        MPI_Cart_shift(cartesianComm, direction ,displacement,&sourceNeighbour,&destinationNeighbour);
        
        const int faceIndex = 2*direction+orientation;
        if (faceIndex < 2*dimensions-1) {
           std::cout << destinationNeighbour << ",";
        } else {
           std::cout << destinationNeighbour << "] (-1: no neighbour)" << std::endl;
        }
      }
    }
  }

  MPI_Finalize();
}
