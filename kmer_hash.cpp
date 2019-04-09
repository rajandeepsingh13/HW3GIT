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
//std::list<upcxx::global_ptr<kmer_pair> contigsMaster;




//		COPIED FROM HASH_MAP.HPP



struct HashMap {

  size_t my_size;
  int rank_n;

  size_t size() const noexcept;

  upcxx::global_ptr<kmer_pair> data_local;
  upcxx::global_ptr<int> used_local;

  HashMap(size_t size, int nprocs);

  void initialize_localData();

  // Most important functions: insert and retrieve
  // k-mers from the hash table.
  bool insert(const kmer_pair &kmer);
  bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer, int currentRank);

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
  data_local = upcxx::new_array<kmer_pair>(my_size/rank_n);
  for (int i = 0; i < rank_n; i++){
    globalData[i] = broadcast(data_local, i).wait();
  }

  //std::vector <int> used;
  used_local = upcxx::new_array<int>(my_size/rank_n);
  for (int i = 0; i < rank_n; i++){
    globalUsed[i] = broadcast(used_local, i).wait();
  }



}

bool HashMap::insert(const kmer_pair &kmer) {
  uint64_t hash = kmer.hash() % my_size;

  int sizePerProc = (my_size+rank_n-1)/rank_n;
  int sizePerProcLast = my_size - sizePerProc*(rank_n - 1);
  int procBasedOnHash;
  int localSlotID;


  uint64_t probeRank = 0;

  bool success = false;

  do {
  	hash = (hash + probeRank) % my_size;
  	procBasedOnHash = hash / sizePerProc;
  	localSlotID = hash % sizePerProc;

  	if (upcxx::rget(globalUsed[procBasedOnHash] + localSlotID).wait() == 0){
  		upcxx::rput(1, globalUsed[procBasedOnHash] + localSlotID).wait();
	    upcxx::rput(kmer_pair(kmer), globalData[procBasedOnHash] + localSlotID).wait();
	    success = true;
	    break;
	} else {
		probeRank++;

	}	
  } while (!success && probeRank < my_size);
  return success;
}

bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer,  int currentRank) {

  uint64_t hash = key_kmer.hash() % my_size;

  int sizePerProc = (my_size+rank_n-1)/rank_n;
  int sizePerProcLast = my_size - sizePerProc*(rank_n - 1);
  int procBasedOnHash;
  int localSlotID;


  int *ptr_used_local = used_local.local();
  kmer_pair *ptr_data_local = data_local.local();

  uint64_t probeRank = 0;

  bool success = false;
  int localSlotCount;
  localSlotCount = currentRank == rank_n-1 ? sizePerProcLast : sizePerProc;


  do {
  	hash = (hash + probeRank) % my_size;
  	procBasedOnHash = hash / sizePerProc;
  	localSlotID = hash % sizePerProc;

  	if (upcxx::rget(globalUsed[procBasedOnHash] + localSlotID).wait() != 0){
	    val_kmer = upcxx::rget(globalData[procBasedOnHash] + localSlotID).wait();
	    if (val_kmer.kmer == key_kmer) {
        	success = true;
      	}
      	else {
      		probeRank++;
      	}
	} else {
		probeRank++;

	}	
  } while (!success && probeRank < my_size);
  return success;

/*
  do {
  	hash = (hash + probeRank) % my_size;
  	localSlotID = hash % sizePerProc;

  	if (ptr_used_local[localSlotID] != 0){
  		val_kmer = ptr_data_local[localSlotID];
  		if (val_kmer.kmer == key_kmer) {
        	success = true;
      	}
      	else {
      		probeRank++;
      	}
	} else {
		probeRank++;
	}	
  } while (!success && probeRank < localSlotCount);
  return success;
*/

}


int main(int argc, char **argv) {
  upcxx::init();
  //std::cout<<" "<<upcxx::rank_n()<<"\n";
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

  hashmap.initialize_localData();
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
  int sizeSplit = (kmers.size()+upcxx::rank_n()-1)/upcxx::rank_n();
  int startPoint = sizeSplit * upcxx::rank_me();
	  for (int i = startPoint; i<startPoint+sizeSplit||i<kmers.size(); i++) {
	  	kmer_pair kmer = kmers[i];
	    bool success = hashmap.insert(kmer);
	    if (!success) {
	      throw std::runtime_error("Error: HashMap is full!");
	    }

	    if (kmer.backwardExt() == 'F') {
	      start_nodes.push_back(kmer);
	    }
	  }

//std::cout<<" "<<start_nodes.size()<<"\n";

  auto end_insert = std::chrono::high_resolution_clock::now();
  upcxx::barrier();

  double insert_time = std::chrono::duration <double> (end_insert - start).count();
  if (run_type != "test") {
    BUtil::print("Finished inserting in %lf\n", insert_time);
  }
  upcxx::barrier();

  auto start_read = std::chrono::high_resolution_clock::now();


  std::cout<<"\n this is it:"<<upcxx::rget(globalUsed[0] + 2).wait();
  upcxx::barrier();


std::list <std::list <kmer_pair>> contigs;
	  
	  for (const auto &start_kmer : start_nodes) {
	    std::list <kmer_pair> contig;
	    contig.push_back(start_kmer);
	    while (contig.back().forwardExt() != 'F') {
	      kmer_pair kmer;
	      bool success = hashmap.find(contig.back().next_kmer(), kmer, upcxx::rank_me()); //fix this
	      if (!success) {
	        throw std::runtime_error("Error: k-mer not found in hashmap.");
	        //int lol = 1;
	      }//else{
	      	//std::cout<<"Success\n";
	      	//contig.push_back(kmer);
	     // }
	      contig.push_back(kmer);
	      //upcxx::barrier();
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