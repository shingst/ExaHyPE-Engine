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

#ifndef EXAHYPE_EXAHYPE_SOLVERS_OUTCOMEDATABASE_H_
#define EXAHYPE_EXAHYPE_SOLVERS_OUTCOMEDATABASE_H_

#if defined(SharedTBB)
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_queue.h>

#include <cstdio>

namespace exahype {
namespace solvers {

enum class DeliveryStatus{Received, Transit};

/*class OutcomeKey {
public:
  OutcomeKey();
  virtual ~OutcomeKey() = 0;
  virtual bool operator==(const OutcomeKey &other) const = 0;
  virtual operator std::size_t() const = 0;
};*/

template <typename T>
class OutcomeEntry {

private:
  T* _data;
  DeliveryStatus _status;
public:
  OutcomeEntry(T* data, DeliveryStatus status);

  T* getData() const;
  DeliveryStatus getStatus() const;

  virtual ~OutcomeEntry();
};

template <typename K, typename T>
class OutcomeDatabase {

private:
  tbb::concurrent_hash_map<K, OutcomeEntry<T>> _database;
public:
  OutcomeDatabase();
  virtual ~OutcomeDatabase();

  void insertOutcome(const K& key, T* data, const DeliveryStatus status);
  bool tryFindAndExtractOutcome(K key, T **data, DeliveryStatus& status);

  int size();
};

} /* namespace solvers */
} /* namespace exahype */

#endif /* SharedTBB*/
#endif /* EXAHYPE_EXAHYPE_SOLVERS_OUTCOMEDATABASE_H_ */
