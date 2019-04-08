#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <set>
#include <numeric>
#include <cstddef>
#include <chrono>
#include <upcxx/upcxx.hpp>

#include "kmer_t.hpp"
#include "read_kmers.hpp"

#include "butil.hpp"

std::vector<upcxx::global_ptr<kmer_pair>> globalData;
std::vector<upcxx::global_ptr<int>> globalUsed;

//		COPIED FROM HASH_MAP.HPP

struct HashMap {

  size_t my_size;
  int rank_n;

  size_t size() const noexcept;

  
  

  HashMap(size_t size, int nprocs);

  void initialize_localData();

  // Most important functions: insert and retrieve
  // k-mers from the hash table.
  bool insert(const kmer_pair &kmer);
  bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

  // Helper functions

  // Write and read to a logical data slot in the table.
  void write_slot(uint64_t slot, const kmer_pair &kmer);
  kmer_pair read_slot(uint64_t slot);

  // Request a slot or check if it's already used.
  bool request_slot(uint64_t slot);
  bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size, int nprocs) {
  my_size = size;
  rank_n = nprocs;
  globalData = std::vector<upcxx::global_ptr<kmer_pair>>(rank_n);      // This can fuck me up as data is declared later
  globalUsed = std::vector<upcxx::global_ptr<int>>(rank_n, 0);   // This can fuck me up as used is declared later
}

void HashMap::initialize_localData(){
  //std::vector <kmer_pair> data;
  upcxx::global_ptr<kmer_pair> data_local = new_array<kmer_pair>(my_size/rank_n);
  for (int i = 0; i < rank_n; i++){
    globalHashMap[i] = broadcast(data_local, i).wait();
  }

  //std::vector <int> used;
  upcxx::global_ptr<int> used_local = new_array<int>(my_size/rank_n);
  for (int i = 0; i < rank_n; i++){
    globalHashMap[i] = broadcast(used_local, i).wait();
  }
}

bool HashMap::insert(const kmer_pair &kmer) {
  uint64_t hash = kmer.hash();

  int sizePerProc = my_size/rank_n;
  int sizePerProcLast = my_size%rank_n;
  int procBasedOnHash = hash / (my_size/rank_n);
  int localSlotID = hash - procBasedOnHash * (my_size/rank_n);


  uint64_t probeRank = 0;

  bool success = false;
  int localSlotCount;
  do {
  	if (rget(globalUsed[procBasedOnHash] + localSlotID).wait() == 0){
  		rput(1, globalUsed[procBasedOnHash] + localSlotID).wait();
	    rput(kmer_pair(kmer), globalData[procBasedOnHash] + localSlotID).wait();
	    success = true;
	    break;
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

bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
  uint64_t hash = key_kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    if (slot_used(slot)) {
      val_kmer = read_slot(slot);
      if (val_kmer.kmer == key_kmer) {
        success = true;
      }
    }
  } while (!success && probe < size());
  return success;
}

bool HashMap::slot_used(uint64_t slot) {
  return used[slot] != 0;
}

void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
  data[slot] = kmer;
}

kmer_pair HashMap::read_slot(uint64_t slot) {
  return data[slot];

int main(int argc, char **argv) {
  upcxx::init();

  // TODO: remove this, when you start writing
  // parallel implementation.
  //if (upcxx::rank_n() > 1) {
  //  throw std::runtime_error("Error: parallel implementation not started yet!"
  //    " (remove this when you start working.)");
  //}

  if (argc < 2) {
    BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test]\n");
    upcxx::finalize();
    exit(1);
  }

  std::string kmer_fname = std::string(argv[1]);
  std::string run_type = "";

  if (argc >= 3) {
    run_type = std::string(argv[2]);
  }

  int ks = kmer_size(kmer_fname);

  if (ks != KMER_LEN) {
    throw std::runtime_error("Error: " + kmer_fname + " contains " +
      std::to_string(ks) + "-mers, while this binary is compiled for " +
      std::to_string(KMER_LEN) + "-mers.  Modify packing.hpp and recompile.");
  }

  size_t n_kmers = line_count(kmer_fname);

  // Load factor of 0.5
  size_t hash_table_size = n_kmers * (1.0 / 0.5);
  HashMap hashmap(hash_table_size, upcxx::rank_n());

  hashmap.initialize_localData()
  upcxx::barrier();

  if (run_type == "verbose") {
    BUtil::print("Initializing hash table of size %d for %d kmers.\n",
      hash_table_size, n_kmers);
  }

  std::vector <kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

  if (run_type == "verbose") {
    BUtil::print("Finished reading kmers.\n");
  }

  auto start = std::chrono::high_resolution_clock::now();

  std::vector <kmer_pair> start_nodes;

  if (rank_me() == MASTER){

	  for (auto &kmer : kmers) {
	    bool success = hashmap.insert(kmer);
	    if (!success) {
	      throw std::runtime_error("Error: HashMap is full!");
	    }

	    if (kmer.backwardExt() == 'F') {
	      start_nodes.push_back(kmer);
	    }
	  }

  }

  auto end_insert = std::chrono::high_resolution_clock::now();
  upcxx::barrier();

  double insert_time = std::chrono::duration <double> (end_insert - start).count();
  if (run_type != "test") {
    BUtil::print("Finished inserting in %lf\n", insert_time);
  }
  upcxx::barrier();

  auto start_read = std::chrono::high_resolution_clock::now();

  std::list <std::list <kmer_pair>> contigs;
  for (const auto &start_kmer : start_nodes) {
    std::list <kmer_pair> contig;
    contig.push_back(start_kmer);
    while (contig.back().forwardExt() != 'F') {
      kmer_pair kmer;
      bool success = hashmap.find(contig.back().next_kmer(), kmer);
      if (!success) {
        throw std::runtime_error("Error: k-mer not found in hashmap.");
      }
      contig.push_back(kmer);
    }
    contigs.push_back(contig);
  }

  auto end_read = std::chrono::high_resolution_clock::now();
  upcxx::barrier();
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration <double> read = end_read - start_read;
  std::chrono::duration <double> insert = end_insert - start;
  std::chrono::duration <double> total = end - start;

  int numKmers = std::accumulate(contigs.begin(), contigs.end(), 0,
    [] (int sum, const std::list <kmer_pair> &contig) {
      return sum + contig.size();
    });

  if (run_type != "test") {
    BUtil::print("Assembled in %lf total\n", total.count());
  }

  if (run_type == "verbose") {
    printf("Rank %d reconstructed %d contigs with %d nodes from %d start nodes."
      " (%lf read, %lf insert, %lf total)\n", upcxx::rank_me(), contigs.size(),
      numKmers, start_nodes.size(), read.count(), insert.count(), total.count());
  }

  if (run_type == "test") {
    std::ofstream fout("test_" + std::to_string(upcxx::rank_me()) + ".dat");
    for (const auto &contig : contigs) {
      fout << extract_contig(contig) << std::endl;
    }
    fout.close();
  }

  upcxx::finalize();
  return 0;
}