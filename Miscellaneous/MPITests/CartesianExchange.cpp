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

void parseArguments(int argc, char** argv, int* dimensions, int* maximumMessageSize, int* numberOfTests) {

  if (argc < 4) { // program name is first argument
    std::cerr << "Error: Please call the program the following way: ./CartesianExchange <dimensions> <maximumMessageSize> <numberOfTests>" <<
        " where dimensions, maximumMessageSize, and numberOfTests are all positive integers." << std::endl;
    std::terminate();
  }

  std::istringstream ss1(argv[1]);
  if (!(ss1 >> *dimensions) || *dimensions<=0 ) {
    std::cerr << "Error: First argument 'dimensions' must be positive integer but is: '" << argv[1] << '\n';
    std::terminate();
  }

  std::istringstream ss2(argv[2]);
  if (!(ss2 >> *maximumMessageSize) || *maximumMessageSize<=0 ) {
    std::cerr << "Error: Second argument 'maximumMessageSize' must be positive integer but is: '" << argv[2] << '\n';
    std::terminate();
  }

  std::istringstream ss3(argv[3]);
  if (!(ss3 >> *numberOfTests) || *numberOfTests<=0 ) {
    std::cerr << "Error: Thrid argument 'numberOfTest' must be positive integer but is: '" << argv[3] << '\n';
    std::terminate();
  }
}

typedef std::tuple<int,int,double,double,double> record; // packetSize,numberOfPackets,avg,min,max

void printRecords(std::vector<record>& records,std::string name) {
  std::cout << std::endl;
  std::cout << name << ":" << std::endl;
  std::cout << std::endl;
  std::cout << "| packet size | #packets    | avg time/us | min time/us | max time/us |" << std::endl;
  std::cout << "|-------------|-------------|-------------|-------------|-------------|" << std::endl;
  record& rPrevious = records.front();
  for (record r :  records) {
    if ( (std::get<0>(r) * std::get<1>(r)) !=  (std::get<0>(rPrevious) * std::get<1>(rPrevious)) ) {
      std::cout << "|-------------|-------------|-------------|-------------|-------------|" << std::endl;
    }

    std::cout <<  "| " << std::setw(11) << std::get<0>(r) << " ";
    std::cout <<  "| " << std::setw(11) << std::get<1>(r) << " ";

    std::stringstream t;
    t << std::fixed << std::setprecision(2) << std::get<2>(r);
    std::cout << "| " << std::setw(11) << t.str() << " ";

    t.str("");
    t << std::fixed << std::setprecision(2) << std::get<3>(r);
    std::cout << "| " << std::setw(11) << t.str() << " ";

    t.str("");
    t << std::fixed << std::setprecision(2) << std::get<4>(r);
    std::cout << "| " << std::setw(11) << t.str() << " ";

    std::cout << "|" << std::endl;

    rPrevious = r;
  }
}


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
  parseArguments(argc,argv,&dimensions,&maximumMessageSize,&numberOfTests);

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
    std::cout << std::endl;

    std::cout << "Start experiment with parameters: "    << std::endl
        << std::endl
        << "dimensions         = " << dimensions         << std::endl
        << "maximumMessageSize = " << maximumMessageSize << std::endl
        << "dimensions         = " << numberOfTests      << std::endl;
  }

  // create a container for neighbouring ranks
  #ifdef UseVector
  std::map<int,std::vector<MPI_Request*>> sendRequests;
  std::map<int,std::vector<MPI_Request*>> receiveRequests;
  #else
  std::map<int,std::list<MPI_Request*>> sendRequests;
  std::map<int,std::list<MPI_Request*>> receiveRequests;
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
  int messageSize = std::round( std::pow( 2, l) );
  while ( messageSize <= maximumMessageSize ) { // loop over message sizes
    double* sendBuffer    = new double[messageSize];
    double* receiveBuffer = new double[messageSize];

    for (int p = 0; p < l+1; ++p) {            // loop over packet sizes
      int packetSize       = std::round( std::pow( 2, p) );
      int numberOfPackets = messageSize/packetSize;

      tPosting.push_back( std::make_tuple(packetSize,numberOfPackets,0.0,99999999.0,0.0) );
      tClear.push_back( std::make_tuple(packetSize,numberOfPackets,0.0,99999999.0,0.0) );
      tWait.push_back( std::make_tuple(packetSize,numberOfPackets,0.0,99999999.0,0.0) );
      for (int test = 0; test < numberOfTests; ++test) {

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
              #ifdef UseVector
              sendRequests.insert( std::pair<int,std::vector<MPI_Request*>>( destinationNeighbour, std::vector<MPI_Request*>(0) ) );
              receiveRequests.insert( std::pair<int,std::vector<MPI_Request*>>( destinationNeighbour, std::vector<MPI_Request*>(0) ) );
              sendRequests[destinationNeighbour].reserve(numberOfPackets);
              receiveRequests[destinationNeighbour].reserve(numberOfPackets);
              #else
              sendRequests.insert( std::pair<int,std::list<MPI_Request*>>( destinationNeighbour, std::list<MPI_Request*>(0) ) );
              receiveRequests.insert( std::pair<int,std::list<MPI_Request*>>( destinationNeighbour, std::list<MPI_Request*>(0) ) );
              #endif

              for (int m = 0; m < numberOfPackets; ++m) {
                MPI_Request* sendRequest    = new MPI_Request();
                MPI_Request* receiveRequest = new MPI_Request();

                // TODO(Dominic): Post sends and receives at the same time vs. post sends first and then receives
                MPI_Isend(sendBuffer,    packetSize, MPI_DOUBLE, destinationNeighbour, 0, cartesianComm, sendRequest);
                MPI_Irecv(receiveBuffer, packetSize, MPI_DOUBLE, destinationNeighbour, 0, cartesianComm, receiveRequest);

                sendRequests[destinationNeighbour].push_back(sendRequest);
                receiveRequests[destinationNeighbour].push_back(receiveRequest);
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
        bool complete = false;
        int flag      = 0;
        while (!complete) {
          complete = true;
          for (auto rankIt = sendRequests.begin(); rankIt != sendRequests.end(); rankIt++) {
            for (auto requestIt = sendRequests[rankIt->first].begin(); requestIt != sendRequests[rankIt->first].end(); requestIt++) {
              MPI_Test(*requestIt,&flag,MPI_STATUS_IGNORE);
              complete &= flag;
            }
            for (auto requestIt = receiveRequests[rankIt->first].begin(); requestIt != receiveRequests[rankIt->first].end(); requestIt++) {
              MPI_Test(*requestIt,&flag,MPI_STATUS_IGNORE);
              complete &= flag;
            }
          }
        }
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
    }

    // clean up
    delete[] sendBuffer;
    delete[] receiveBuffer;

    // go into next iteration
    l++;
    messageSize = std::round( std::pow( 2, l) );
  }

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
        "Info: Measurements obtained for rank "<< myRank <<
        " located at x_i="<<ranksPerDimension[0]/2<<", i=0,..,"<<dimensions-1<< "." << std::endl;
  }

  MPI_Finalize();
}
