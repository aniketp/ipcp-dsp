// Code analyzed for the Third Data Prefetching Championship
//
// Author: Alberto Ros, University of Murcia
//
// Paper #13: Berti: A Per-Page Best-Request-Time Delta Prefetcher

#include "cache.h"

void CACHE::llc_prefetcher_initialize() 
{
    cout << "LLC Next Line or No Prefetcher" << endl;
}

uint32_t CACHE::llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{
  if (NUM_CPUS == 1) {
    uint64_t pf_addr = ((addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;
    prefetch_line(ip, addr, pf_addr, FILL_LLC, 0);
  }
  
  return metadata_in;
}

uint32_t CACHE::llc_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
  return metadata_in;
}

void CACHE::llc_prefetcher_final_stats()
{
  cout << "LLC Next Line or No Prefetcher Final Stats: none" << endl;
}
