#include "cache.h"
#define MIN(a,b) (((a)<(b))?(a):(b))

#define word_offset 3 // (observing byte addresses with a 64-bit alignment)

// Delta Cache 
#define DC_range_l1d (128*8)	// 1024 (-1) possible deltas
#define DC_ways_l1d 16			// Delta cache assoc.
#define DC_LFUmax_l1d 128		// Maximum LFU counter value

int DC_deltanext_l1d[DC_range_l1d][DC_ways_l1d]; // Next deltas
int DC_LFUbits_l1d[DC_range_l1d][DC_ways_l1d];	 // Frequency counters


// Page Cache
#define PC_sets_l1d 256		// 256 sets
#define PC_ways_l1d 12		// 12 ways
#define PC_tag_bits_l1d 10 	// 10-bit page tags

int PC_ldelta_l1d[PC_sets_l1d][PC_ways_l1d]; 	// Last deltas
int PC_loffset_l1d[PC_sets_l1d][PC_ways_l1d];	// Last offsets
int PC_ptag_l1d[PC_sets_l1d][PC_ways_l1d];		// Page tag
int PC_NRUbit_l1d[PC_sets_l1d][PC_ways_l1d];	// Not-Recently Used (NRU) bits


// Supplementary function to caclulate the page tag
int get_page_tag_l1d(uint64_t page){
	return (page/PC_sets_l1d)&((1<<PC_tag_bits_l1d)-1);
}


// Prefetcher initialisation
void CACHE::l1d_prefetcher_initialize() 
{
	printf("Ultra pref. initializing...\n"); fflush(stdout);

	// Initialise the Delta Cache 
	for (int i=0; i<DC_range_l1d; i++){ 
		for (int j=0; j<DC_ways_l1d; j++){ 
			DC_deltanext_l1d[i][j]= 1 + DC_range_l1d/2; // (Optional (currently not used): fallback to delta=1)
			DC_LFUbits_l1d[i][j] = 0;
		}
	}
		
	// Note: the Page Cache initialisation is less important, since we are looking for page tag hits of 10-bits
	printf("Ultra pref. initialized\n"); fflush(stdout);
}


// This function updates the Delta Cache with a new delta transition (delta_from -> delta_next)
void update_DC_l1d (int delta_from, int delta_to){

	// Look for hits
	int dhit = 0;
	for (int i=0; i<DC_ways_l1d; i++){
		
		// If there is a hit, increment the respective counter
		if (DC_deltanext_l1d[delta_from][i]==delta_to){
			DC_LFUbits_l1d[delta_from][i]++;
			
			// If there is an overflow, 
			if (DC_LFUbits_l1d[delta_from][i]==DC_LFUmax_l1d){
				for (int j=0; j<DC_ways_l1d; j++){
					
					// Decrement all counters in the set (frequency proportions shall remain about the same)
					DC_LFUbits_l1d[delta_from][j]/=2;
				}
			}		
			dhit=1;			
			break;
		}
	}
	
	// If the delta transition is not in Delta Cache,
	if (dhit==0){
	
		// Evict the least-frequent delta transition in the set
		int min_freq=DC_LFUbits_l1d[delta_from][0];
		int min_freq_way=0;		
		for (int i=1; i<DC_ways_l1d; i++){
			if (DC_LFUbits_l1d[delta_from][i] < min_freq){
				min_freq = DC_LFUbits_l1d[delta_from][i];
				min_freq_way = i;
			}
		}
		
		// And replace with the current one		
		DC_deltanext_l1d[delta_from][min_freq_way] = delta_to;
		DC_LFUbits_l1d [delta_from][min_freq_way] = 1; 
	}
}


// This function returns the most probably immidiately next delta based on the current delta
int get_next_best_transition_l1d (int delta){

	// Caclulate the sum of the LFU counters for the current set
	int probs_sum = 0;
	for (int i=0; i<DC_ways_l1d; i++){
		probs_sum += DC_LFUbits_l1d[delta][i];
	}
	
	// Find the maximum LFU value
	int max_freq=DC_LFUbits_l1d[delta][0];
	int max_freq_way=0;	
	for (int i=1; i<DC_ways_l1d; i++){
		if (DC_LFUbits_l1d[delta][i] > max_freq){
			max_freq = DC_LFUbits_l1d[delta][i];
			max_freq_way = i;
		}
	}
	
	// Discard, if it represents a probability lower than 1/3
	if ((float)DC_LFUbits_l1d[delta][max_freq_way]/probs_sum<1/3.0) 
		return -1;
		
	return DC_deltanext_l1d[delta][max_freq_way];	
}


void CACHE::l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type)
{

	unsigned long long int ref_addr = addr>>word_offset;
	unsigned long long int block 	  = addr>>LOG2_BLOCK_SIZE;
	unsigned long long int page = addr >>(6+LOG2_BLOCK_SIZE);
	int page_offset = ref_addr & ( (1<<(6+(LOG2_BLOCK_SIZE-word_offset))) -1); // (Last 9 bits, for a 64-bit alignment)
	
	// Page Cache lookup: Looking for previous entries of the same page in the Page Cache
	int way=-1;
	for (int i=0; i<PC_ways_l1d; i++){
		if (PC_ptag_l1d[page%PC_sets_l1d][i]==get_page_tag_l1d(page)){
			way=i;
			break;
		}
	}

	int cur_delta = 1 + DC_range_l1d/2;  // (fallback to delta=1, when there is no page match)
	int matched = 0;	
	
	// If there was a Page Cache hit, 
	if (way!=-1){
		int ldelta_l1d=PC_ldelta_l1d[page%PC_sets_l1d][way];
		int loff_l1d=PC_loffset_l1d[page%PC_sets_l1d][way];		
		
		// Calculate current delta,
		cur_delta = page_offset-loff_l1d + DC_range_l1d/2;
		matched=1;
	
		// And update the Delta cache with the new delta transition
		update_DC_l1d(ldelta_l1d, cur_delta);		
	} 

	int next_delta=cur_delta;	
	uint64_t addr_n=addr;
	int count=0;
	
	// Note: Since we are requesting byte addresses, the degree here is not the actual prefetch degree (PQ length is 8), (The line limit is enforced by the simulator)
	int degree = 36; // (PQ.SIZE-PQ.occupancy)*2/3; uint64_t line_requests[PQ.SIZE*8]; // Uncomment for counting lines instead
	
	// Decrease the degree for the multi-core configuration
	if (NUM_CPUS>1) degree/=4;
	
	for (int i_=0; (i_<128) && (count<degree); i_++){				
			
		// Find the most probable next delta from the current delta (1 generation) 
		int best_delta = get_next_best_transition_l1d(next_delta);
		
		// Abort if probability less than 1/3
		if (best_delta==-1) break;
					
		{													
			// Aggregate all counter values in the set
			int sum=0;
			for (int j=0; j<DC_ways_l1d; j++){
				sum+=DC_LFUbits_l1d[next_delta][j];
			}
			
			// Looking for the top 2 child deltas (no more than 2 having a probability > 1/3)
			int used[DC_ways_l1d] = {0};			
			for (int i=0; i<2/*DC_ways_l1d*/; i++){ 
				int max_way = -1;
				int max_value =  -1;
				for (int j=0; j<DC_ways_l1d; j++){
					if((DC_LFUbits_l1d[next_delta][j]>max_value) && (!used[j])){
						max_way = j;
						max_value = DC_LFUbits_l1d[next_delta][j];													
					}
				}											
				if (max_way==-1) break; 
				
				// If the probability is greater than 1/3,
				if((float)DC_LFUbits_l1d[next_delta][max_way]/sum > 1/3.0){
					used[max_way]=1;
					
					uint64_t pf_addr = ((addr_n>>word_offset)+(DC_deltanext_l1d[next_delta][max_way]-DC_range_l1d/2)) << word_offset;
					uint64_t pf_block = pf_addr >> LOG2_BLOCK_SIZE;
					unsigned long long int pf_page = pf_addr>>12;

					// And it falls in the same page, but in a diferent cache line
					if ((page==pf_page) && (block!=pf_block)) {						
					
						//// Uncomment for counting cache lines instead
						//int already_prefetched=0;
						//for (int k=0; k<count; k++){
						//	if (line_requests[k]==pf_block){
						//		already_prefetched=1; break;
						//	}
						//}	
						//if (!already_prefetched){	
						//	line_requests[count]=pf_block;
						
							// Prefetch block
							prefetch_line(ip, addr, pf_addr /*pf_block<<LOG2_BLOCK_SIZE*/, FILL_L1, 0);
							count++;
							
							// Stop if the prefetch degree is reached
							if (count==degree) break;
						//}
					}						
				} 			
			}
		}
		
		// Update values for moving to the next delta generation based on the top next delta	
		next_delta = best_delta;
		uint64_t pf_addr = ((addr_n>>word_offset)+(best_delta-DC_range_l1d/2)) << word_offset;
		addr_n = pf_addr;		
	}

	// If there was a Page cache miss, evict the Not-Recently used
	if (way==-1) {		
	
		// Look for NRU bit equal to 0 
		for (int i=0; i<PC_ways_l1d; i++){
			if (PC_NRUbit_l1d[page%PC_sets_l1d][i]==0){
				way=i;
				break;
			}				
		}
		
		// If all are equal to 1, flip them 
		if (way==-1){
			way=0;
			for (int i=0; i<PC_ways_l1d; i++)
				PC_NRUbit_l1d[page%PC_sets_l1d][i]=0;
		}
	}


	// Update the respective Page Cache entry
	if (matched)
		PC_ldelta_l1d[page%PC_sets_l1d][way]=cur_delta; 
	else	
		// If we have this entry for the first time, the delta value is invalid 
		// (0 represents delta=-512, which falls in a different page)		
		PC_ldelta_l1d[page%PC_sets_l1d][way]=0;
		
	PC_loffset_l1d[page%PC_sets_l1d][way]=page_offset;
	PC_ptag_l1d[page%PC_sets_l1d][way]=get_page_tag_l1d(page);
	PC_NRUbit_l1d[page%PC_sets_l1d][way]=1;
}


// Not used
void CACHE::l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
	unsigned long long int ref_addr = addr>>6;
	unsigned long long int page = ref_addr>>6;
	int page_offset = ref_addr&63;
}


void CACHE::l1d_prefetcher_final_stats()
{
	// Print the information stored in Delta Cache
	printf("\n");
	for (int i=0; i<DC_range_l1d; i++){
		printf ("L1D %d>\t",i-DC_range_l1d/2);
		for (int j=0; j<DC_ways_l1d; j++){
			printf("(%d,%d)\t", DC_deltanext_l1d[i][j]-DC_range_l1d/2,DC_LFUbits_l1d[i][j]);
		}		
		printf("\n");	
	}	
}
