// Submission For DPC-3 Team 12
// L1D Prefetcher is throttled next-N-lines in single-core and no prefetcher in multi-core

#include "cache.h"

#define BAD_SCORE 0
#define MAX_SCORE 30
#define SCORE_ENTRIES 1024 
// Keeping a score of 1024 pages. Score is a measure of being NextLine prefetch friendly

#define PAGES_TRACKED 6 // Tracking distinct pages in last 6 accesses

int score_table[SCORE_ENTRIES];
int history[SCORE_ENTRIES];

uint64_t page_tracker[PAGES_TRACKED];

void CACHE::l1d_prefetcher_initialize() 
{
    cout << "CPU " << cpu << " L1D next line prefetcher" << endl;
	for(int a = 0; a < SCORE_ENTRIES; a++){
		score_table[a] = 0;
		history[a] = -1;
	}
}

void CACHE::l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type)
{

// Single Core configuration
if (NUM_CPUS == 1) {

	uint64_t page = addr >> LOG2_PAGE_SIZE; 
	for (int i = PAGES_TRACKED-1; i>0; i--) { // N down to 1
		page_tracker[i] = page_tracker[i-1];
	}
	page_tracker[0] = page; // Keeping history of last 8 accesses

	int distinct_pages = 0;
	for (int i=0; i < PAGES_TRACKED; i++) {
		int j;
		for (j=0; j<i; j++) {
			if (page_tracker[i] == page_tracker[j])
				break;
		}
		if (i==j)
			distinct_pages++;
	} // Distinct pages being accessed in last-6 accesses


	int index = (int)((addr >> 12) % SCORE_ENTRIES);
	int delta = 0; 

	if(history[index] == -1){
		history[index] = (int)((addr >> LOG2_BLOCK_SIZE) & 0x3F);
		return;
	}else{
		delta = (((addr >> LOG2_BLOCK_SIZE) & 0x3F) - history[index]);	
		
		if(delta == 1 && abs(score_table[index]) != MAX_SCORE)
				score_table[index]++;
		else if(abs(score_table[index]) > BAD_SCORE) 
			score_table[index]--;
	 
		history[index] = (int)((addr >> LOG2_BLOCK_SIZE) & 0x3F); 

		if(score_table[index] > BAD_SCORE){ 
			//Page is NextLine prefetch friendly
		
			uint64_t pf_addr = addr;
			for (int i = 0; i < ceil(PQ.SIZE / distinct_pages)-1; i++) {
				//
				pf_addr = ((pf_addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;
    			prefetch_line(ip, addr, pf_addr, FILL_L1, 3);
			}
			pf_addr = ((pf_addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;
    		prefetch_line(ip, addr, pf_addr, FILL_L1, 6);
			// Last Prefetch with a different metadata
			// LLC Prefetcher doesn't ignore prefetch requests with this metadata value
			// In a way, continuing where L1D left off
		}
	}

    DP ( if (warmup_complete[cpu]) {
    cout << "[" << NAME << "] " << __func__ << hex << " base_cl: " << (addr>>LOG2_BLOCK_SIZE);
    cout << " pf_cl: " << (pf_addr>>LOG2_BLOCK_SIZE) << " ip: " << ip << " cache_hit: " << +cache_hit << " type: " << +type << endl; });
}

// No L1D Prefetcher in multi-core
if (NUM_CPUS == 4) {
}
}

void CACHE::l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{

}

void CACHE::l1d_prefetcher_final_stats()
{
}

