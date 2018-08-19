#include "SHMPinningObserver.h"
#include "SHMMacros.h"


#include <iostream>
#include <cassert>
#include <sys/resource.h>
#include <sched.h>
#include <stdlib.h>
#include <unistd.h>


shminvade::SHMPinningObserver::SHMPinningObserver():
  _availableCores(0),
  _mutex() {
}


void shminvade::SHMPinningObserver::observe(bool toggle) {
  assert( _availableCores.count()==0 );

  // reconstruct Unix mask
  cpu_set_t* mask;
  for ( _numberOfCPUs = sizeof(cpu_set_t)/CHAR_BIT; _numberOfCPUs < MaxCores; _numberOfCPUs <<= 1 ) {
    mask = CPU_ALLOC( _numberOfCPUs );
    if ( !mask ) break;
    const size_t size = CPU_ALLOC_SIZE( _numberOfCPUs );
    CPU_ZERO_S( size, mask );
    const int err = sched_getaffinity( 0, size, mask );
    if ( !err ) break;

    CPU_FREE( mask );
    mask = NULL;
    if ( errno != EINVAL )  break;
  }

  // Convert into bitset to make it more C++ish
  if ( mask ) {
    _availableCores = 0;
    for (int i=0; i<MaxCores; i++) {
      if (CPU_ISSET(i, mask)) {
        _availableCores[i] = true;
      }
    }

    #if SHM_INVADE_DEBUG>=2
    std::cout << SHM_DEBUG_PREFIX << "identified available cores: " << _availableCores << std::endl;
    #endif

    tbb::task_scheduler_observer::observe(toggle);
  }
  else {
    std::cerr << SHM_DEBUG_PREFIX << SHM_DEBUG_SEPARATOR << "failed to obtain process affinity mask" << std::endl;
  }
}


shminvade::SHMPinningObserver::~SHMPinningObserver() {
}


void shminvade::SHMPinningObserver::on_scheduler_entry( bool ) {
  _mutex.lock();

  if (_availableCores.count()==0) {
    std::cerr << SHM_DEBUG_PREFIX << SHM_DEBUG_SEPARATOR << "too many threads, i.e. no idle core available anymore" << std::endl;
  }
  else {
    int targetCore = 0;
    for (int i=0; i<MaxCores; i++) {
      if (_availableCores[i]) {
    	targetCore = i;
        _availableCores[i] = false;
        i = MaxCores;
      }
    }

	cpu_set_t*   target_mask = CPU_ALLOC( _numberOfCPUs );
    const size_t size        = CPU_ALLOC_SIZE( _numberOfCPUs );
    CPU_ZERO_S( size, target_mask );
    CPU_SET_S( targetCore, size, target_mask );
	const int err = sched_setaffinity( 0, size, target_mask );

    if ( err ) {
      std::cerr << SHM_DEBUG_PREFIX << SHM_DEBUG_SEPARATOR << "pinning new thread to core " << targetCore << " failed with error code " << err << std::endl;
    }
    else {
      const int   currentCore     = sched_getcpu();
      #if SHM_INVADE_DEBUG>=2
      std::cout << SHM_DEBUG_PREFIX << SHM_DEBUG_SEPARATOR << "pinned new thread currently running on core " << currentCore << " to core " << targetCore << std::endl;
      #endif
    }

    CPU_FREE( target_mask );
  }

  _mutex.unlock();
}


void shminvade::SHMPinningObserver::on_scheduler_exit( bool ) {
  _mutex.lock();

  const int   currentCore     = sched_getcpu();
  if (_availableCores[currentCore]) {
    #if SHM_INVADE_DEBUG>=2
    std::cerr << SHM_DEBUG_PREFIX << SHM_DEBUG_SEPARATOR << "thread running on core " << currentCore << " exits though core is already marked as idle" << std::endl;
    #endif
  }
  else {
    _availableCores[currentCore] = true;
    #if SHM_INVADE_DEBUG>=2
    std::cout << SHM_DEBUG_PREFIX << SHM_DEBUG_SEPARATOR << "thread exits and frees core " << currentCore << std::endl;
    #endif
  }

  _mutex.unlock();
}

