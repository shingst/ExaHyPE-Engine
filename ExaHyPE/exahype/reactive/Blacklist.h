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

#if !defined(_EXAHYPE_BLACKLIST_H_)  && defined(Parallel)
#define _EXAHYPE_BLACKLIST_H_

#include "tarch/logging/Log.h"

namespace exahype {
namespace reactive{

/**
 * Singleton for the blacklist.
 * Ranks from which offloaded tasks arrive to late are put on the blacklist.
 * The blacklist is passed on to the PerformanceMonitor which distributes
 * the blacklist information to other ranks.
 */
class Blacklist {
  private:
  /**
   * The logging device.
   */
  static tarch::logging::Log _log;

  //singleton
  Blacklist();
  virtual ~Blacklist();

  /**
   * Local instance of blacklist.
   */
  double *_localBlacklist;

  public:

  static Blacklist& getInstance();

  /**
   * @param MPI rank
   * @return True if rank is blacklisted
   */
  bool isBlacklisted(int rank) const;

  /**
   * Notifies the blacklist that a task has come too late from a victim rank and that the victim rank
   * needs to be blacklisted.
   * Should be called whenever an emergency arises.
   * @param Rank which needs to be put on the blacklist.
   */
  void triggerEmergencyAndBlacklistRank(int rank);

  /**
   * If called, decreases blacklist values such that blacklisted ranks eventually recover.
   */
  void recoverBlacklistedRanks();

  /**
   * Prints blacklist.
   */
  void printBlacklist();
};

}
}


#endif
