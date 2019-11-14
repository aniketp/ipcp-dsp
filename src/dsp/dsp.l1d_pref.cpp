#include "cache.h"
#include "ooo_cpu.h"
#include <algorithm> 

#define STREAM_TABLE_SIZE 8
#define MAX_L1_PREFETCH_DEGREE 8
#define MAX_L1_PREFETCH_DIST 0
int PREF_DEGREE = 4;
int PREF_DIST = 0;
int STREAM_DIR = 0;
int num_ins = 0;
int num_hits = 0;
float upbnd = 0.75;
float lwbnd = 0.25;
int monitor_time = 1000;
std::vector<uint64_t> stream_table;
std::vector<uint64_t> st_lru_counter;
int prefetch_acc[8] = {0};
uint64_t prefetched_addr[8] = {0};
uint64_t last_blk_addr = -1;
uint64_t previous_core_cycles = 0;
uint64_t previous_num_retired = 0;

void update_degree(int strength, int inc_acc){
	// cout << strength <<":"<< PREF_DEGREE<< endl;
	num_hits = 0;
	int inc_ss =0;
	if (strength >= 0.5*STREAM_TABLE_SIZE)
	{
		inc_ss = 2;
	}
	else if( strength >= 0.25*STREAM_TABLE_SIZE && strength < 0.5*STREAM_TABLE_SIZE)
		inc_ss = 1;
	else if( strength >= 0.125*STREAM_TABLE_SIZE && strength < 0.25*STREAM_TABLE_SIZE)
		inc_ss = 0;
	else
		inc_ss = -1;

	// cout<< inc_ss << ":" << inc_acc << endl;
	if(inc_acc == 1 && inc_ss == 2)
		PREF_DEGREE += 2;
	else if(inc_acc == 0 && inc_ss == 2)
		PREF_DEGREE += 1;
	else if(inc_acc == -1 && inc_ss == 2)
		PREF_DEGREE += 0;
	else if(inc_acc == 1 && inc_ss == 1)
		PREF_DEGREE += 1;
	else if(inc_acc == 0 && inc_ss == 1)
		PREF_DEGREE += 1;
	else if(inc_acc == -1 && inc_ss == 1)
		PREF_DEGREE += 0;
	else if(inc_acc == 1 && inc_ss == 0)
		PREF_DEGREE += 1;
	else if(inc_acc == 0 && inc_ss == 0)
		PREF_DEGREE += 0;
	else if(inc_acc == -1 && inc_ss == 0)
		PREF_DEGREE += 0;
	else if(inc_acc == 1 && inc_ss == -1)
		PREF_DEGREE += 0;
	else if(inc_acc == 0 && inc_ss == -1)
		PREF_DEGREE += -1;
	else if(inc_acc == -1 && inc_ss == -1)
		PREF_DEGREE += -1;
		
	// cout << "IPC-B:" << (float) (ooo_cpu[0].num_retired - previous_num_retired)/(cycle_window) << endl;

}


void CACHE::l1d_prefetcher_initialize() 
{
    for (int i = 0; i < STREAM_TABLE_SIZE; i++) {
        stream_table.push_back(-1);
        st_lru_counter.push_back(0);
    }
    //initialize number of cycles and number of instructions retired
    if(previous_core_cycles == 0)
	    previous_core_cycles = current_core_cycle[0];
	if(previous_num_retired == 0)
		previous_num_retired = ooo_cpu[0].num_retired;
}

//find the table's current size
int var_table_size(){
    int size = 0;
    for (int i = 0; i < STREAM_TABLE_SIZE; i++) {
    	if(stream_table[i] != -1)
    		size++;
    }
    return size;
}

// find lru element's index
int find_min()
{
    int mic = st_lru_counter[0];
    int mindex = 0;
    for(int i=0;i<STREAM_TABLE_SIZE;i++)
    {
        if(st_lru_counter[i]<mic){
        	mic=st_lru_counter[i];
        	mindex = i;
        }
    }
    return mindex;
}

//find block in the stream table
int find_block(uint64_t blk_addr){
	for (int i = 0; i < STREAM_TABLE_SIZE; ++i)
	{
		if(blk_addr == stream_table[i])
			return i;
	}
	return -1;
}


void CACHE::l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type)
{
    uint64_t curr_page = addr >> LOG2_PAGE_SIZE;
    uint64_t curr_blk_addr = addr >> LOG2_BLOCK_SIZE;
    uint64_t pf_addr;
	int inc_degree = 0;
	int pos = 0;
	int neg =0;
	int stream_strength =0;

	// intervention at every 1000 cycle interval
	if(current_core_cycle[0] - previous_core_cycles >= monitor_time){
		previous_core_cycles = current_core_cycle[0];
		previous_num_retired = ooo_cpu[0].num_retired;
		num_ins=0;
	}
	if (last_blk_addr == -1)
		last_blk_addr = curr_blk_addr;
	else{
		if(last_blk_addr == curr_blk_addr)
			return; // if last block address same as current block address, do nothing
		last_blk_addr = curr_blk_addr; // update last block address
	}

	//update stream table
	int it = find_block(curr_blk_addr);
	if( (std::find(stream_table.begin(),stream_table.end(),curr_blk_addr)) == stream_table.end())
	{
		int min_index = find_min();
		//evict the lru element
		stream_table.erase(stream_table.begin() + min_index); 
		st_lru_counter.erase(st_lru_counter.begin() + min_index);
		//push back the new element 
		stream_table.push_back(curr_blk_addr);
		st_lru_counter.push_back(0);
	}
	else{
		st_lru_counter[it]++;
	}
	// count memory instructions hit rate for every 1000 cycle window
	if (cache_hit == 1)
    {
		num_hits++;
    }
	if (var_table_size() == STREAM_TABLE_SIZE){
		//check for positive or negative stream
		while(1)
		{
			curr_blk_addr--;
			if((curr_blk_addr << LOG2_BLOCK_SIZE)>>LOG2_PAGE_SIZE != curr_page) // check for all elements till the page boundary
				break;
			if((std::find(stream_table.begin(),stream_table.end(),curr_blk_addr)) != stream_table.end())
				pos++;
		}
		curr_blk_addr = addr >> LOG2_BLOCK_SIZE;//reset current block address
		while(1)
		{
			curr_blk_addr++;
			if((curr_blk_addr << LOG2_BLOCK_SIZE)>>LOG2_PAGE_SIZE != curr_page)// check for all elements till the page boundary
				break;
			if((std::find(stream_table.begin(),stream_table.end(),curr_blk_addr)) != stream_table.end())
				neg++;
		}
		curr_blk_addr = addr >> LOG2_BLOCK_SIZE;//reset current block address
		if (neg > pos){
			STREAM_DIR = -1;
			stream_strength = neg;
		}
		else{
			STREAM_DIR = 1;
			stream_strength = pos;
		}

		//degree and distance increment based on prefetch accuracy
		for (int i = 0; i < PREF_DEGREE; ++i)
		{
			if (curr_blk_addr == prefetched_addr[i] && cache_hit == 1)
			{
				prefetch_acc[i] = 1;
			}
		}
		if (num_ins == 0)
		{
			int sum=0, inc_acc =0;
			for (int i = 0; i < MAX_L1_PREFETCH_DEGREE; ++i)
		    {
		    	sum+=prefetch_acc[i];
		    	prefetch_acc[i] = 0;//flush all accuracy data from the last window
		    	prefetched_addr[i] = 0;//flush all addresses prefetched in the last window
		    }
		    if (sum > upbnd*PREF_DEGREE)
		    {
		    	inc_acc++;
		    }
		    else if(sum < lwbnd*PREF_DEGREE){
		    	inc_acc--;
		    }
		    else
		    	inc_acc = 0;
		    update_degree(stream_strength, inc_acc);

		    if (PREF_DEGREE > MAX_L1_PREFETCH_DEGREE)
		    {
		    	PREF_DEGREE = MAX_L1_PREFETCH_DEGREE;
		    	upbnd = 0.875;
		    	lwbnd = 0.25;
		    }
		    if (PREF_DIST > MAX_L1_PREFETCH_DIST)
		    	PREF_DIST = MAX_L1_PREFETCH_DIST;

		}
		if(PREF_DEGREE < 1)
			PREF_DEGREE = 1;//make it next line
	}	
	//prefetch
	if (num_ins == 0){
	    for (int i = 1; i <= PREF_DEGREE; ++i)
	    {
	    	pf_addr = (curr_blk_addr + STREAM_DIR*(PREF_DIST + i))<<LOG2_BLOCK_SIZE;
	    	if (curr_page != (pf_addr>>LOG2_PAGE_SIZE)){
	    		break;
	    	}
			prefetched_addr[i-1] = (pf_addr) >> LOG2_BLOCK_SIZE;
	    	prefetch_line(ip, addr, pf_addr, FILL_L1, 0);
	    }
	}
	else{
	    pf_addr = ((addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;
	    prefetch_line(ip, addr, pf_addr, FILL_L1, 0);

	}
	num_ins++;
	return;
}

void CACHE::l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
	// cout<<"FILL" << "," << (addr>>LOG2_BLOCK_SIZE) << "," << set << "," << way<< "," << metadata_in << endl;

}

void CACHE::l1d_prefetcher_final_stats()
{

}
