#include "SHMPinningObserver.h"
#include "SHMMacros.h"
#include <sched.h>
#include <sys/resource.h>
#include <iostream>
#include <thread>

namespace {
  //const int MaxNumberOfSupportedCPUs = 16*1024;
  const int MaxNumberOfSupportedCPUs = sizeof(long int)*8;
}



shminvade::SHMPinningObserver::SHMPinningObserver():
  _mask( nullptr ) {

  for ( _ncpus = sizeof(cpu_set_t)/CHAR_BIT; _ncpus < MaxNumberOfSupportedCPUs; _ncpus <<= 1 ) {
    _mask = CPU_ALLOC( _ncpus );
    if ( !_mask ) break;
    const size_t size = CPU_ALLOC_SIZE( _ncpus );
    CPU_ZERO_S( size, _mask );
    const int err = sched_getaffinity( 0, size, _mask );
    if ( !err ) break;

    CPU_FREE( _mask );
    _mask = NULL;
    if ( errno != EINVAL )  break;
  }
  if ( _mask ) {
    #if SHM_INVADE_DEBUG>=1
    CPUSetBitfield bitfield = cpuSetMaskToBitfield( *_mask );
    std::cout << SHM_DEBUG_PREFIX <<  "Process mask is " << bitfield
    		  << ", i.e. system has " << bitfield.count() << " logical cores (second half usually hypercores)"
    		  << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;

    if ( bitfield.count() != std::thread::hardware_concurrency()) {
      std::cerr << SHM_DEBUG_PREFIX <<  "Process affinity mask has only " << bitfield.count() << " entries although hardware concurrency is " << std::thread::hardware_concurrency() << std::endl;
    }
    #endif

    // thread 0 will be pinned twice, but I want to be sure it is in our data
    // base right from the start
    _numThreads++;
    pinCurrentThread();
    _numThreads--;
  }
  else {
    std::cerr << SHM_DEBUG_PREFIX <<  "Failed to obtain process affinity mask" << std::endl;
  }
}


shminvade::SHMPinningObserver::~SHMPinningObserver() {
  if ( _mask != nullptr ) {
    CPU_FREE( _mask );
  }
}


void shminvade::SHMPinningObserver::pinCurrentThread() {
  const size_t size = CPU_ALLOC_SIZE( _ncpus );
  const int num_cpus = CPU_COUNT_S( size, _mask );
  int thr_idx =  tbb::task_arena::current_thread_index();
  thr_idx %= num_cpus; // To limit unique number in [0; num_cpus-1] range
  // Place threads with specified step
  int cpu_idx = 0;
  for ( int i = 0, offset = 0; i<thr_idx; ++i ) {
    cpu_idx ++;
    if ( cpu_idx >= num_cpus )
      cpu_idx = ++offset;
  }


  // Find index of 'cpu_idx'-th bit equal to 1
  int mapped_idx = -1;
  while ( cpu_idx >= 0 ) {
    if ( CPU_ISSET_S( ++mapped_idx, size, _mask ) )
      --cpu_idx;
  }

  cpu_set_t *target_mask = CPU_ALLOC( _ncpus );
  CPU_ZERO_S( size, target_mask );
  CPU_SET_S( mapped_idx, size, target_mask );
  const int err = sched_setaffinity( 0, size, target_mask );

  if ( err ) {
    std::cerr << SHM_DEBUG_PREFIX <<  "Failed to set thread affinity!" << std::endl;
    exit( EXIT_FAILURE );
  }

  CPU_FREE( target_mask );
}


void shminvade::SHMPinningObserver::on_scheduler_entry( bool ) {
  ++_numThreads;

  if ( _mask ) {
    pinCurrentThread();
  }
}


void shminvade::SHMPinningObserver::on_scheduler_exit( bool ) {
  --_numThreads;
}


int shminvade::SHMPinningObserver::getNumberOfRegisteredThreads() const {
  return _numThreads;
}


shminvade::SHMPinningObserver::CPUSetBitfield shminvade::SHMPinningObserver::cpuSetMaskToBitfield( cpu_set_t   _mask) {
  CPUSetBitfield result = 0;

  for (long i = 0; i < MaxNumberOfSupportedCPUs; i++) {
    if (CPU_ISSET(i, &_mask)) {
      result[i] = true;
    }
  }

  return result;
}

