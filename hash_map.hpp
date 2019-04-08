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


bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer,  int currentRank) {

  uint64_t hash = key_kmer.hash() % my_size;

  int sizePerProc = (my_size+rank_n-1)/rank_n;
  int sizePerProcLast = my_size - sizePerProc*(rank_n - 1);
  int procBasedOnHash;
  int localSlotID;


  uint64_t probeRank = 0;

  bool success = false;
  int localSlotCount;

  do {
    hash = (hash + probeRank) % my_size;
    procBasedOnHash = hash / sizePerProc;
    localSlotID = hash % sizePerProc;

    if (rget(globalUsed[procBasedOnHash] + localSlotID).wait() != 0){
      val_kmer = rget(globalData[procBasedOnHash] + localSlotID).wait();
      if (val_kmer.kmer == key_kmer) {
          success = true;
        }
        else {
          if(procBasedOnHash == rank_n-1){
        localSlotCount = my_size%rank_n;
      } else{
        localSlotCount = my_size/rank_n;
      }

      if (localSlotID<localSlotCount)
        localSlotID++;
      else{
        localSlotID=0;
        procBasedOnHash++;
        probeRank++;
      }

        }
  } else {

    if(procBasedOnHash == rank_n-1){
      localSlotCount = my_size%rank_n;
    } else{
      localSlotCount = my_size/rank_n;
    }

    if (localSlotID<localSlotCount)
      localSlotID++;
    else{
      localSlotID=0;
      procBasedOnHash++;
      probeRank++;
    }

  } 
  } while (!success && probeRank < rank_n);
  return success;


}
*/