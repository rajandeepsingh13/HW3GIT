#pragma once

#include <upcxx/upcxx.hpp>
#include "kmer_t.hpp"

/*

bool HashMap::insert(const kmer_pair &kmer) {
  uint64_t hash = kmer.hash();



  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    success = request_slot(slot);
    if (success) {
      write_slot(slot, kmer);
    }
  } while (!success && probe < size());
  return success;
}








bool HashMap::request_slot(uint64_t slot) {
  if (used[slot] != 0) {
    return false;
  } else {
    used[slot] = 1;
    return true;
  }
}

size_t HashMap::size() const noexcept {
  return my_size;
}
*/