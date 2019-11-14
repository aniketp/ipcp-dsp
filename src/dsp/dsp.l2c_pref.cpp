#include "cache.h"
#include "ooo_cpu.h"
#include <algorithm> 

#define L2_STREAM_TABLE_SIZE 16
#define MAX_L2_PREFETCH_DEGREE 12
#define MAX_L2_PREFETCH_DIST 0
int L2_PREF_DEGREE = 4;
int L2_PREF_DIST = 2;
int L2_STREAM_DIR = 0;
int L2_num_ins = 0;
int L2_num_hits = 0;
float L2_upbnd = 0.75;
float L2_lwbnd = 0.25;
int L2_monitor_time = 1000;
std::vector<uint64_t> L2_stream_table;
std::vector<uint64_t> L2_st_lru_counter;
int L2_prefetch_acc[12] = {0};
uint64_t L2_prefetched_addr[12] = {0};
uint64_t L2_last_blk_addr = -1;
uint64_t L2_previous_core_cycles = 0;
uint64_t L2_previous_num_retired = 0;

void L2_update_degree(int strength, int inc_acc){
	// cout << strength <<":"<< L2_num_hits <<":" << L2_PREF_DEGREE<< endl;
	L2_num_hits = 0;
	int inc_ss =0;
	if (strength >= 0.5*L2_STREAM_TABLE_SIZE)
	{
		inc_ss = 2;
	}
	else if( strength >= 0.25*L2_STREAM_TABLE_SIZE && strength < 0.5*L2_STREAM_TABLE_SIZE)
		inc_ss = 1;
	else if( strength >= 0.125*L2_STREAM_TABLE_SIZE && strength < 0.25*L2_STREAM_TABLE_SIZE)
		inc_ss = 0;
	else
		inc_ss = -1;

	// cout<< inc_ss << ":" << inc_acc << endl;
	if(inc_acc == 1 && inc_ss == 2)
		L2_PREF_DEGREE += 2;
	else if(inc_acc == 0 && inc_ss == 2)
		L2_PREF_DEGREE += 1;
	else if(inc_acc == -1 && inc_ss == 2)
		L2_PREF_DEGREE += 0;
	else if(inc_acc == 1 && inc_ss == 1)
		L2_PREF_DEGREE += 1;
	else if(inc_acc == 0 && inc_ss == 1)
		L2_PREF_DEGREE += 1;
	else if(inc_acc == -1 && inc_ss == 1)
		L2_PREF_DEGREE += 0;
	else if(inc_acc == 1 && inc_ss == 0)
		L2_PREF_DEGREE += 1;
	else if(inc_acc == 0 && inc_ss == 0)
		L2_PREF_DEGREE += 0;
	else if(inc_acc == -1 && inc_ss == 0)
		L2_PREF_DEGREE += 0;
	else if(inc_acc == 1 && inc_ss == -1)
		L2_PREF_DEGREE += 0;
	else if(inc_acc == 0 && inc_ss == -1)
		L2_PREF_DEGREE += -1;
	else if(inc_acc == -1 && inc_ss == -1)
		L2_PREF_DEGREE += -1;

	// cout << "IPC-B:" << (float) (ooo_cpu[0].num_retired - L2_previous_num_retired)/(cycle_window) << endl;

}
int L2_var_table_size(){
    int size = 0;
    for (int i = 0; i < L2_STREAM_TABLE_SIZE; i++) {
    	if(L2_stream_table[i] != -1)
    		size++;
    }
    return size;
}

// find lru element's index
int L2_find_min()
{
    int mic = L2_st_lru_counter[0];
    int mindex = 0;
    for(int i=0;i<L2_STREAM_TABLE_SIZE;i++)
    {
        if(L2_st_lru_counter[i]<mic){
        	mic=L2_st_lru_counter[i];
        	mindex = i;
        }
    }
    return mindex;
}

//find block in the stream table
int L2_find_block(uint64_t blk_addr){
	for (int i = 0; i < L2_STREAM_TABLE_SIZE; ++i)
	{
		if(blk_addr == L2_stream_table[i])
			return i;
	}
	return -1;
}


void CACHE::l2c_prefetcher_initialize() 
{
    for (int i = 0; i < L2_STREAM_TABLE_SIZE; i++) {
        L2_stream_table.push_back(-1);
        L2_st_lru_counter.push_back(0);
    }
    if(L2_previous_core_cycles == 0)
	    L2_previous_core_cycles = current_core_cycle[0];
	if(L2_previous_num_retired == 0)
		L2_previous_num_retired = ooo_cpu[0].num_retired;
}

uint32_t CACHE::l2c_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{
    uint64_t curr_page = addr >> LOG2_PAGE_SIZE;
    uint64_t curr_blk_addr = addr >> LOG2_BLOCK_SIZE;
    uint64_t pf_addr;
	int inc_degree = 0;
	int pos = 0;
	int neg =0;
	int stream_strength =0;

	if(current_core_cycle[0] - L2_previous_core_cycles >= L2_monitor_time){
		L2_previous_core_cycles = current_core_cycle[0];
		L2_previous_num_retired = ooo_cpu[0].num_retired;
		L2_num_ins=0;
	}
	if (L2_last_blk_addr == -1)
		L2_last_blk_addr = curr_blk_addr;
	else{
		if(L2_last_blk_addr == curr_blk_addr)
			return metadata_in; // if last block address same as current block address, do nothing
		L2_last_blk_addr = curr_blk_addr; // update last block address
	}

	//update stream table
	int it = L2_find_block(curr_blk_addr);
	if( (std::find(L2_stream_table.begin(),L2_stream_table.end(),curr_blk_addr)) == L2_stream_table.end())
	{
		int min_index = L2_find_min();
		//evict the lru element
		L2_stream_table.erase(L2_stream_table.begin() + min_index); 
		L2_st_lru_counter.erase(L2_st_lru_counter.begin() + min_index);
		//push back the new element 
		L2_stream_table.push_back(curr_blk_addr);
		L2_st_lru_counter.push_back(0);
	}
	else{
		L2_st_lru_counter[it]++;
	}
	if (cache_hit == 1)
    {
		L2_num_hits++;
    }
	if (L2_var_table_size() == L2_STREAM_TABLE_SIZE){
		//check for positive or negative stream
		while(1)
		{
			curr_blk_addr--;
			if((curr_blk_addr << LOG2_BLOCK_SIZE)>>LOG2_PAGE_SIZE != curr_page) // check for all elements till the page boundary
				break;
			if((std::find(L2_stream_table.begin(),L2_stream_table.end(),curr_blk_addr)) != L2_stream_table.end())
				pos++;
		}
		curr_blk_addr = addr >> LOG2_BLOCK_SIZE;//reset current block address
		while(1)
		{
			curr_blk_addr++;
			if((curr_blk_addr << LOG2_BLOCK_SIZE)>>LOG2_PAGE_SIZE != curr_page)// check for all elements till the page boundary
				break;
			if((std::find(L2_stream_table.begin(),L2_stream_table.end(),curr_blk_addr)) != L2_stream_table.end())
				neg++;
		}
		curr_blk_addr = addr >> LOG2_BLOCK_SIZE;//reset current block address
		if (neg > pos){
			L2_STREAM_DIR = -1;
			stream_strength = neg;
		}
		else{
			L2_STREAM_DIR = 1;
			stream_strength = pos;
		}

		//degree and distance increment based on prefetch accuracy
		for (int i = 0; i < L2_PREF_DEGREE; ++i)
		{
			if (curr_blk_addr == L2_prefetched_addr[i] && cache_hit == 1)
			{
				L2_prefetch_acc[i] = 1;
			}
		}
		if (L2_num_ins == 0)
		{
			int sum=0, inc_acc =0;
			for (int i = 0; i < MAX_L2_PREFETCH_DEGREE; ++i)
		    {
		    	sum+=L2_prefetch_acc[i];
		    	L2_prefetch_acc[i] = 0;//flush all accuracy data from the last window
		    	L2_prefetched_addr[i] = 0;//flush all addresses prefetched in the last window
		    }
		    if (sum > L2_upbnd*L2_PREF_DEGREE)
		    {
		    	inc_acc++;
		    }
		    else if(sum < L2_lwbnd*L2_PREF_DEGREE){
		    	inc_acc--;
		    }
		    else
		    	inc_acc = 0;
		    L2_update_degree(stream_strength, inc_acc);

		    if (L2_PREF_DEGREE > MAX_L2_PREFETCH_DEGREE)
		    {
		    	L2_PREF_DEGREE = MAX_L2_PREFETCH_DEGREE;
		    	L2_upbnd = 0.875;
		    	L2_lwbnd = 0.25;
		    }
		    if (L2_PREF_DIST > MAX_L2_PREFETCH_DIST)
		    	L2_PREF_DIST = MAX_L2_PREFETCH_DIST;

		}
		if(L2_PREF_DEGREE < 1)
			L2_PREF_DEGREE = 1;//make it next line
	}	
	//prefetch
	if (L2_num_ins == 0){
	    for (int i = 1; i <= L2_PREF_DEGREE; ++i)
	    {
	    	pf_addr = (curr_blk_addr + L2_STREAM_DIR*(L2_PREF_DIST + i))<<LOG2_BLOCK_SIZE;
	    	if (curr_page != (pf_addr>>LOG2_PAGE_SIZE)){
	    		break;
	    	}
			L2_prefetched_addr[i-1] = (pf_addr) >> LOG2_BLOCK_SIZE;
	    	prefetch_line(ip, addr, pf_addr, FILL_L2, 0);
	    }
	}
	else{
	    pf_addr = ((addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;
	    prefetch_line(ip, addr, pf_addr, FILL_L2, 0);

	}
	L2_num_ins++;
    return metadata_in;
}

uint32_t CACHE::l2c_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
  return metadata_in;
}

void CACHE::l2c_prefetcher_final_stats()
{
    cout << "CPU " << cpu << " L2C next line prefetcher final stats" << endl;
}
