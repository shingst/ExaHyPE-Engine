#include "SHMMacros.h"

#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <signal.h>


std::string shminvade::getSHMDebugPrefix() {
  static const pid_t process = (pid_t) syscall (__NR_getpid);
  return std::string(SHM_DEBUG_PREFIX) + std::string(SHM_DEBUG_SEPARATOR) + std::to_string(process) + std::string(SHM_DEBUG_SEPARATOR);
}



