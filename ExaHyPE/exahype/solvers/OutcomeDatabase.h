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
#endif

#include <cstdio>

namespace exahype {
namespace solvers {

/**
 * Status flag that indicates whether an outcome is still in transit or has already been received.
 */
enum class DeliveryStatus{Received, Transit};

/**
 * Template class for storing an outcome entry of type T.
 * Each outcome has a flag that indicates whether it has fully been received.
 */
template <typename T>
class OutcomeEntry {

private:
  /**
   * The actual data.
   */
  T* _data;

  /**
   * The delivery status.
   */
  DeliveryStatus _status;
public:
  /**
   * Creates an outcome entry.
   */
  OutcomeEntry(T* data, DeliveryStatus status);

  /**
   * Returns a pointer to the data.
   */
  T* getData() const;

  /**
   * Returns the delivery status of this entry.
   */
  DeliveryStatus getStatus() const;

  virtual ~OutcomeEntry();
};

/**
 * Outcome database.
 * The template argument K is a key used to search for an outcome.
 * The template argument T is the type of the outcome.
 */
template <typename K, typename T>
class OutcomeDatabase {

private:
  /**
   * Outcomes are currently stored in a concurrent hash map.
   * We need TBB.
   */
#if defined(SharedTBB)
  tbb::concurrent_hash_map<K, OutcomeEntry<T>> _database;
#endif
public:
  OutcomeDatabase();
  ~OutcomeDatabase();

  /**
   * Inserts an outcome with key K into the database.
   * @param key A unique identifier that is hashable.
   * @param data The data to be inserted.
   * @param status The status of the outcome entry.
   *
   */
  void insertOutcome(const K& key, T* data, const DeliveryStatus status);

  /**
   * Searches for an outcome and extracts it if possible (calling function may want to re-insert).
   * @return Returns true if the outcome was found.
   */
  bool tryFindAndExtractOutcome(K key, T **data, DeliveryStatus& status);

  /**
   * Returns the unsafe size of the job outcome database.
   */
  int size() const;
};

} /* namespace solvers */
} /* namespace exahype */

#endif /* EXAHYPE_EXAHYPE_SOLVERS_OUTCOMEDATABASE_H_ */
