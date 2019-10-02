// Submission For DPC-3 Team 12
// L2C Prefetcher is PPF in both the configurations: single-core or multi-core

#ifndef SPP_H
#define SPP_H

#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <cmath>
#include "cache.h"

using namespace std;

// SPP functional knobs
#define LOOKAHEAD_ON
#define FILTER_ON
#define GHR_ON
#define SPP_SANITY_CHECK

//#define SPP_DEBUG_PRINT
#ifdef SPP_DEBUG_PRINT
#define SPP_DP(x) x
#else
#define SPP_DP(x)
#endif

//#define SPP_PERC_WGHT
#ifdef SPP_PERC_WGHT
#define SPP_PW(x) x
#else 
#define SPP_PW(x)
#endif

// Signature table parameters
#define ST_SET 1
#define ST_WAY 256
#define ST_TAG_BIT 16
#define ST_TAG_MASK ((1 << ST_TAG_BIT) - 1)
#define SIG_SHIFT 3
#define SIG_BIT 12
#define SIG_MASK ((1 << SIG_BIT) - 1)
#define SIG_DELTA_BIT 7

// Pattern table parameters
#define PT_SET 2048
#define PT_WAY 4
#define C_SIG_BIT 4
#define C_DELTA_BIT 4
#define C_SIG_MAX ((1 << C_SIG_BIT) - 1)
#define C_DELTA_MAX ((1 << C_DELTA_BIT) - 1)

// Prefetch filter parameters
#define QUOTIENT_BIT  10
#define REMAINDER_BIT 6
#define HASH_BIT (QUOTIENT_BIT + REMAINDER_BIT + 1)
#define FILTER_SET (1 << QUOTIENT_BIT)

#define QUOTIENT_BIT_REJ  10
#define REMAINDER_BIT_REJ 8
#define HASH_BIT_REJ (QUOTIENT_BIT_REJ + REMAINDER_BIT_REJ + 1)
#define FILTER_SET_REJ (1 << QUOTIENT_BIT_REJ)

// Global register parameters
#define GLOBAL_COUNTER_BIT 10
#define GLOBAL_COUNTER_MAX ((1 << GLOBAL_COUNTER_BIT) - 1) 
#define MAX_GHR_ENTRY 8
#define PAGES_TRACKED 6

// Perceptron paramaters
#define PERC_ENTRIES 4096 //Upto 12-bit addressing in hashed perceptron
#define PERC_FEATURES 9 //Keep increasing based on new features
#define PERC_COUNTER_MAX 15 //-16 to +15: 5 bits counter 
#define PERC_THRESHOLD_HI  -5
#define PERC_THRESHOLD_LO  -15
#define POS_UPDT_THRESHOLD  90
#define NEG_UPDT_THRESHOLD -80

enum FILTER_REQUEST {SPP_L2C_PREFETCH, SPP_LLC_PREFETCH, L2C_DEMAND, L2C_EVICT, SPP_PERC_REJECT}; // Request type for prefetch filter

uint64_t get_hash(uint64_t key);

class SIGNATURE_TABLE {
  public:
    bool     valid[ST_SET][ST_WAY];
    uint32_t tag[ST_SET][ST_WAY],
             last_offset[ST_SET][ST_WAY],
             sig[ST_SET][ST_WAY],
             lru[ST_SET][ST_WAY];

    SIGNATURE_TABLE() {
        cout << "Initialize SIGNATURE TABLE" << endl;
        cout << "ST_SET: " << ST_SET << endl;
        cout << "ST_WAY: " << ST_WAY << endl;
        cout << "ST_TAG_BIT: " << ST_TAG_BIT << endl;
        cout << "ST_TAG_MASK: " << hex << ST_TAG_MASK << dec << endl;

        for (uint32_t set = 0; set < ST_SET; set++)
            for (uint32_t way = 0; way < ST_WAY; way++) {
                valid[set][way] = 0;
                tag[set][way] = 0;
                last_offset[set][way] = 0;
                sig[set][way] = 0;
                lru[set][way] = way;
            }
    };

    void read_and_update_sig(uint64_t page, uint32_t page_offset, uint32_t &last_sig, uint32_t &curr_sig, int32_t &delta);
};

class PATTERN_TABLE {
  public:
    int      delta[PT_SET][PT_WAY];
    uint32_t c_delta[PT_SET][PT_WAY],
             c_sig[PT_SET];

    PATTERN_TABLE() {
        cout << endl << "Initialize PATTERN TABLE" << endl;
        cout << "PT_SET: " << PT_SET << endl;
        cout << "PT_WAY: " << PT_WAY << endl;
        cout << "SIG_DELTA_BIT: " << SIG_DELTA_BIT << endl;
        cout << "C_SIG_BIT: " << C_SIG_BIT << endl;
        cout << "C_DELTA_BIT: " << C_DELTA_BIT << endl;

        for (uint32_t set = 0; set < PT_SET; set++) {
            for (uint32_t way = 0; way < PT_WAY; way++) {
                delta[set][way] = 0;
                c_delta[set][way] = 0;
            }
            c_sig[set] = 0;
        }
    }

    void update_pattern(uint32_t last_sig, int curr_delta),
         read_pattern(uint32_t curr_sig, int *prefetch_delta, uint32_t *confidence_q, int32_t *perc_sum_q, uint32_t &lookahead_way, uint32_t &lookahead_conf, uint32_t &pf_q_tail, uint32_t &depth, uint64_t addr, uint64_t base_addr, uint64_t train_addr, uint64_t curr_ip, int32_t train_delta, uint32_t last_sig, uint32_t pq_occupancy, uint32_t pq_SIZE, uint32_t mshr_occupancy, uint32_t mshr_SIZE);
};

class PREFETCH_FILTER {
  public:
    uint64_t remainder_tag[FILTER_SET],
			 pc[FILTER_SET],
			 pc_1[FILTER_SET],
			 pc_2[FILTER_SET],
			 pc_3[FILTER_SET],
			 address[FILTER_SET];
    bool     valid[FILTER_SET],  // Consider this as "prefetched"
             useful[FILTER_SET]; // Consider this as "used"
	int32_t	 delta[FILTER_SET],
			 perc_sum[FILTER_SET];
	uint32_t last_signature[FILTER_SET],
			 confidence[FILTER_SET],
			 cur_signature[FILTER_SET],
			 la_depth[FILTER_SET];

	uint64_t remainder_tag_reject[FILTER_SET_REJ],
			 pc_reject[FILTER_SET_REJ],
			 pc_1_reject[FILTER_SET_REJ],
			 pc_2_reject[FILTER_SET_REJ],
			 pc_3_reject[FILTER_SET_REJ],
			 address_reject[FILTER_SET_REJ];
	bool 	 valid_reject[FILTER_SET_REJ]; // Entries which the perceptron rejected
	int32_t	 delta_reject[FILTER_SET_REJ],
			 perc_sum_reject[FILTER_SET_REJ];
	uint32_t last_signature_reject[FILTER_SET_REJ],
			 confidence_reject[FILTER_SET_REJ],
			 cur_signature_reject[FILTER_SET_REJ],
		     la_depth_reject[FILTER_SET_REJ];

    PREFETCH_FILTER() {
        cout << endl << "Initialize PREFETCH FILTER" << endl;
        cout << "FILTER_SET: " << FILTER_SET << endl;

        for (uint32_t set = 0; set < FILTER_SET; set++) {
        	remainder_tag[set] = 0;
        	valid[set] = 0;
        	useful[set] = 0;
		}
        for (uint32_t set = 0; set < FILTER_SET_REJ; set++) {
        	valid_reject[set] = 0;
        	remainder_tag_reject[set] = 0;
		}
    }

    bool     check(uint64_t pf_addr, uint64_t base_addr, uint64_t ip, FILTER_REQUEST filter_request, int32_t cur_delta, uint32_t last_sign, uint32_t cur_sign, uint32_t confidence, int32_t sum, uint32_t depth);
};

class PERCEPTRON {
  public:
	// Perc Weights
	int32_t perc_weights[PERC_ENTRIES][PERC_FEATURES];

	// CONST depths for different features
	int32_t PERC_DEPTH[PERC_FEATURES];

	PERCEPTRON() {
		cout << "\nInitialize PERCEPTRON" << endl;
		cout << "PERC_ENTRIES: " << PERC_ENTRIES << endl;
		cout << "PERC_FEATURES: " << PERC_FEATURES << endl;

		PERC_DEPTH[0] = 2048;   //base_addr;
		PERC_DEPTH[1] = 4096;   //cache_line;
		PERC_DEPTH[2] = 4096;  	//page_addr;
		PERC_DEPTH[3] = 4096;   //confidence ^ page_addr;
		PERC_DEPTH[4] = 1024;	//curr_sig ^ sig_delta;
		PERC_DEPTH[5] = 4096; 	//ip_1 ^ ip_2 ^ ip_3;		
		PERC_DEPTH[6] = 1024; 	//ip ^ depth;
		PERC_DEPTH[7] = 2048;   //ip ^ sig_delta;
		PERC_DEPTH[8] = 128;   	//confidence;

		for (int i = 0; i < PERC_ENTRIES; i++) {
			for (int j = 0;j < PERC_FEATURES; j++) {
				perc_weights[i][j] = 0;
			}
		}
	}

	void	 perc_update(uint64_t check_addr, uint64_t ip, uint64_t ip_1, uint64_t ip_2, uint64_t ip_3, int32_t cur_delta, uint32_t last_sig, uint32_t curr_sig, uint32_t confidence, uint32_t depth, bool direction, int32_t perc_sum);
	int32_t	perc_predict(uint64_t check_addr, uint64_t ip, uint64_t ip_1, uint64_t ip_2, uint64_t ip_3, int32_t cur_delta, uint32_t last_sig, uint32_t curr_sig, uint32_t confidence, uint32_t depth);
};


class GLOBAL_REGISTER {
  public:
    // Global counters to calculate global prefetching accuracy
    uint64_t pf_useful,
             pf_issued,
             global_accuracy; // Alpha value in Section III. Equation 3

    // Global History Register (GHR) entries
    uint8_t  valid[MAX_GHR_ENTRY];
    uint32_t sig[MAX_GHR_ENTRY],
             confidence[MAX_GHR_ENTRY],
             offset[MAX_GHR_ENTRY];
    int      delta[MAX_GHR_ENTRY];

	uint64_t ip_0,
			 ip_1,
			 ip_2,
			 ip_3;

	uint64_t page_tracker[PAGES_TRACKED];
	
    GLOBAL_REGISTER() {
        pf_useful = 0;
        pf_issued = 0;
        global_accuracy = 0;
		ip_0 = 0;
		ip_1 = 0;
		ip_2 = 0;
		ip_3 = 0;

        for (uint32_t i = 0; i < MAX_GHR_ENTRY; i++) {
            valid[i] = 0;
            sig[i] = 0;
            confidence[i] = 0;
            offset[i] = 0;
            delta[i] = 0;
        }
    }

    void update_entry(uint32_t pf_sig, uint32_t pf_confidence, uint32_t pf_offset, int pf_delta);
    uint32_t check_entry(uint32_t page_offset);
};

#endif



SIGNATURE_TABLE ST;
PATTERN_TABLE   PT;
PREFETCH_FILTER FILTER;
GLOBAL_REGISTER GHR;
PERCEPTRON PERC;

void CACHE::l2c_prefetcher_initialize() 
{

}

uint32_t CACHE::l2c_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{

    uint64_t page = addr >> LOG2_PAGE_SIZE;
    uint32_t page_offset = (addr >> LOG2_BLOCK_SIZE) & (PAGE_SIZE / BLOCK_SIZE - 1),
             last_sig = 0,
             curr_sig = 0,
             confidence_q[100*L2C_MSHR_SIZE],
             depth = 0;

    int32_t  delta = 0,
             delta_q[100*L2C_MSHR_SIZE],
			 perc_sum_q[100*L2C_MSHR_SIZE];

    for (uint32_t i = 0; i < 100*L2C_MSHR_SIZE; i++){
        confidence_q[i] = 0;
        delta_q[i] = 0;
		perc_sum_q[i] = 0;
    }
    confidence_q[0] = 100;
    GHR.global_accuracy = GHR.pf_issued ? ((100 * GHR.pf_useful) / GHR.pf_issued)  : 0;
   	
	for (int i = PAGES_TRACKED-1; i>0; i--) { // N down to 1
		GHR.page_tracker[i] = GHR.page_tracker[i-1];
	}
	GHR.page_tracker[0] = page;

	int distinct_pages = 0;
	uint8_t num_pf = 0;
	for (int i=0; i < PAGES_TRACKED; i++) {
		int j;
		for (j=0; j<i; j++) {
			if (GHR.page_tracker[i] == GHR.page_tracker[j])
				break;
		}
		if (i==j)
			distinct_pages++;
	}
	//cout << "Distinct Pages: " << distinct_pages << endl;
	
    SPP_DP (
        cout << endl << "[ChampSim] " << __func__ << " addr: " << hex << addr << " cache_line: " << (addr >> LOG2_BLOCK_SIZE);
        cout << " page: " << page << " page_offset: " << dec << page_offset << endl;
    );

    // Stage 1: Read and update a sig stored in ST
    // last_sig and delta are used to update (sig, delta) correlation in PT
    // curr_sig is used to read prefetch candidates in PT 
    ST.read_and_update_sig(page, page_offset, last_sig, curr_sig, delta);
    
	// Also check the prefetch filter in parallel to update global accuracy counters 
    FILTER.check(addr, 0, 0, L2C_DEMAND, 0, 0, 0, 0, 0, 0); 

    // Stage 2: Update delta patterns stored in PT
    if (last_sig) PT.update_pattern(last_sig, delta);

    // Stage 3: Start prefetching
    uint64_t base_addr = addr;
	uint64_t curr_ip = ip;
    uint32_t lookahead_conf = 100,
             pf_q_head = 0, 
             pf_q_tail = 0;
    uint8_t  do_lookahead = 0;
	int32_t  prev_delta = 0;

	uint64_t train_addr  = addr;
	int32_t  train_delta = 0;

	GHR.ip_3 = GHR.ip_2;
	GHR.ip_2 = GHR.ip_1;
	GHR.ip_1 = GHR.ip_0;
	GHR.ip_0 = ip;

#ifdef LOOKAHEAD_ON
    do {
#endif
        uint32_t lookahead_way = PT_WAY;

		train_addr  = addr; train_delta = prev_delta;
		// Remembering the original addr here and accumulating the deltas in lookahead stages
		
		// Read the PT. Also passing info required for perceptron inferencing as PT calls perc_predict()
        PT.read_pattern(curr_sig, delta_q, confidence_q, perc_sum_q, lookahead_way, lookahead_conf, pf_q_tail, depth, addr, base_addr, train_addr, curr_ip, train_delta, last_sig, PQ.occupancy, PQ.SIZE, MSHR.occupancy, MSHR.SIZE);

        do_lookahead = 0;
        for (uint32_t i = pf_q_head; i < pf_q_tail; i++) {

			uint64_t pf_addr = (base_addr & ~(BLOCK_SIZE - 1)) + (delta_q[i] << LOG2_BLOCK_SIZE);
			int32_t perc_sum   = perc_sum_q[i];

			SPP_DP(
				cout << "[ChampSim] State of features: \nTrain addr: " << train_addr << "\tCurr IP: " << curr_ip << "\tIP_1: " << GHR.ip_1 << "\tIP_2: " << GHR.ip_2 << "\tIP_3: " << GHR.ip_3 << "\tDelta: " << train_delta + delta_q[i] << "\t:LastSig " << last_sig << "\t:CurrSig " << curr_sig << "\t:Conf " << confidence_q[i] << "\t:Depth " << depth << "\tSUM: "<< perc_sum  << endl;
			);
			FILTER_REQUEST fill_level = (perc_sum >= PERC_THRESHOLD_HI) ? SPP_L2C_PREFETCH : SPP_LLC_PREFETCH;
			
            if ((addr & ~(PAGE_SIZE - 1)) == (pf_addr & ~(PAGE_SIZE - 1))) { // Prefetch request is in the same physical page
				
				// Filter checks for redundancy and returns FALSE if redundant
				// Else it returns TRUE and logs the features for future retrieval 
				if ( num_pf < ceil(((PQ.SIZE)/distinct_pages)) ) {					
                	if (FILTER.check(pf_addr, train_addr, curr_ip, fill_level, train_delta + delta_q[i], last_sig, curr_sig, confidence_q[i], perc_sum, (depth-1))) {

						//[DO NOT TOUCH]:	
						// Use addr (not base_addr) to obey the same physical page boundary
						prefetch_line(ip, addr, pf_addr, ((fill_level == SPP_L2C_PREFETCH) ? FILL_L2 : FILL_LLC),5); 
						num_pf++;
						
						//FILTER.valid_reject[quotient] = 0;
                    	if (fill_level == SPP_L2C_PREFETCH) {
                    	    GHR.pf_issued++;
                    	    if (GHR.pf_issued > GLOBAL_COUNTER_MAX) {
                    	        GHR.pf_issued >>= 1;
                    	        GHR.pf_useful >>= 1;
                    	    }
                    	    SPP_DP (cout << "[ChampSim] SPP L2 prefetch issued GHR.pf_issued: " << GHR.pf_issued << " GHR.pf_useful: " << GHR.pf_useful << endl;);
                    	}

                    	SPP_DP (
                    	    cout << "[ChampSim] " << __func__ << " base_addr: " << hex << base_addr << " pf_addr: " << pf_addr;
                    	    cout << " pf_cache_line: " << (pf_addr >> LOG2_BLOCK_SIZE);
                    	    cout << " prefetch_delta: " << dec << delta_q[i] << " confidence: " << confidence_q[i];
                    	    cout << " depth: " << i << " fill_level: " << ((fill_level == SPP_L2C_PREFETCH) ? FILL_L2 : FILL_LLC) << endl;
                    	);
                	}
				}	
            } else { // Prefetch request is crossing the physical page boundary
#ifdef GHR_ON
                    // Store this prefetch request in GHR to bootstrap SPP learning when we see a ST miss (i.e., accessing a new page)
                    GHR.update_entry(curr_sig, confidence_q[i], (pf_addr >> LOG2_BLOCK_SIZE) & 0x3F, delta_q[i]); 
#endif
            }
            do_lookahead = 1;
            pf_q_head++;
        }

		// Update base_addr and curr_sig
        if (lookahead_way < PT_WAY) {
            uint32_t set = get_hash(curr_sig) % PT_SET;
            base_addr += (PT.delta[set][lookahead_way] << LOG2_BLOCK_SIZE);
			prev_delta += PT.delta[set][lookahead_way]; 

            // PT.delta uses a 7-bit sign magnitude representation to generate sig_delta
            //int sig_delta = (PT.delta[set][lookahead_way] < 0) ? ((((-1) * PT.delta[set][lookahead_way]) & 0x3F) + 0x40) : PT.delta[set][lookahead_way];
            int sig_delta = (PT.delta[set][lookahead_way] < 0) ? (((-1) * PT.delta[set][lookahead_way]) + (1 << (SIG_DELTA_BIT - 1))) : PT.delta[set][lookahead_way];
            curr_sig = ((curr_sig << SIG_SHIFT) ^ sig_delta) & SIG_MASK;
        }

        SPP_DP (
            cout << "Looping curr_sig: " << hex << curr_sig << " base_addr: " << base_addr << dec;
            cout << " pf_q_head: " << pf_q_head << " pf_q_tail: " << pf_q_tail << " depth: " << depth << endl;
        );
#ifdef LOOKAHEAD_ON
    } while (do_lookahead);
#endif

	return metadata_in;
}

uint32_t CACHE::l2c_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t match, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{

#ifdef FILTER_ON
    SPP_DP (cout << endl;);
    FILTER.check(evicted_addr, 0, 0, L2C_EVICT, 0, 0, 0, 0, 0, 0);
#endif

	return metadata_in;
}

void CACHE::l2c_prefetcher_final_stats()
{

}

// TODO: Find a good 64-bit hash function
uint64_t get_hash(uint64_t key)
{
    // Robert Jenkins' 32 bit mix function
    key += (key << 12);
    key ^= (key >> 22);
    key += (key << 4);
    key ^= (key >> 9);
    key += (key << 10);
    key ^= (key >> 2);
    key += (key << 7);
    key ^= (key >> 12);

    // Knuth's multiplicative method
    key = (key >> 3) * 2654435761;

    return key;
}

void SIGNATURE_TABLE::read_and_update_sig(uint64_t page, uint32_t page_offset, uint32_t &last_sig, uint32_t &curr_sig, int32_t &delta)
{

    uint32_t set = get_hash(page) % ST_SET,
             match = ST_WAY,
             partial_page = page & ST_TAG_MASK;
    uint8_t  ST_hit = 0;
    int      sig_delta = 0;

    SPP_DP (cout << "[ST] " << __func__ << " page: " << hex << page << " partial_page: " << partial_page << dec << endl;);

    // Case 1: Hit
    for (match = 0; match < ST_WAY; match++) {
        if (valid[set][match] && (tag[set][match] == partial_page)) {
            last_sig = sig[set][match];
            delta = page_offset - last_offset[set][match];

            if (delta) {
                // Build a new sig based on 7-bit sign magnitude representation of delta
                //sig_delta = (delta < 0) ? ((((-1) * delta) & 0x3F) + 0x40) : delta;
                sig_delta = (delta < 0) ? (((-1) * delta) + (1 << (SIG_DELTA_BIT - 1))) : delta;
                sig[set][match] = ((last_sig << SIG_SHIFT) ^ sig_delta) & SIG_MASK;
                curr_sig = sig[set][match];
                last_offset[set][match] = page_offset;

                SPP_DP (
                    cout << "[ST] " << __func__ << " hit set: " << set << " way: " << match;
                    cout << " valid: " << valid[set][match] << " tag: " << hex << tag[set][match];
                    cout << " last_sig: " << last_sig << " curr_sig: " << curr_sig;
                    cout << " delta: " << dec << delta << " last_offset: " << page_offset << endl;
                );
            } else last_sig = 0; // Hitting the same cache line, delta is zero

            ST_hit = 1;
            break;
        }
    }

    // Case 2: Invalid
    if (match == ST_WAY) {
        for (match = 0; match < ST_WAY; match++) {
            if (valid[set][match] == 0) {
                valid[set][match] = 1;
                tag[set][match] = partial_page;
                sig[set][match] = 0;
                curr_sig = sig[set][match];
                last_offset[set][match] = page_offset;

                SPP_DP (
                    cout << "[ST] " << __func__ << " invalid set: " << set << " way: " << match;
                    cout << " valid: " << valid[set][match] << " tag: " << hex << partial_page;
                    cout << " sig: " << sig[set][match] << " last_offset: " << dec << page_offset << endl;
                );

                break;
            }
        }
    }

    // Case 3: Miss
    if (match == ST_WAY) {
        for (match = 0; match < ST_WAY; match++) {
            if (lru[set][match] == ST_WAY - 1) { // Find replacement victim
                tag[set][match] = partial_page;
                sig[set][match] = 0;
                curr_sig = sig[set][match];
                last_offset[set][match] = page_offset;

                SPP_DP (
                    cout << "[ST] " << __func__ << " miss set: " << set << " way: " << match;
                    cout << " valid: " << valid[set][match] << " victim tag: " << hex << tag[set][match] << " new tag: " << partial_page;
                    cout << " sig: " << sig[set][match] << " last_offset: " << dec << page_offset << endl;
                );

                break;
            }
        }

        #ifdef SPP_SANITY_CHECK
        // Assertion
        if (match == ST_WAY) {
            cout << "[ST] Cannot find a replacement victim!" << endl;
            assert(0);
        }
        #endif
    }

#ifdef GHR_ON
    if (ST_hit == 0) {
        uint32_t GHR_found = GHR.check_entry(page_offset);
        if (GHR_found < MAX_GHR_ENTRY) {
            sig_delta = (GHR.delta[GHR_found] < 0) ? (((-1) * GHR.delta[GHR_found]) + (1 << (SIG_DELTA_BIT - 1))) : GHR.delta[GHR_found];
            sig[set][match] = ((GHR.sig[GHR_found] << SIG_SHIFT) ^ sig_delta) & SIG_MASK;
            curr_sig = sig[set][match];
        }
    }
#endif

    // Update LRU
    for (uint32_t way = 0; way < ST_WAY; way++) {
        if (lru[set][way] < lru[set][match]) {
            lru[set][way]++;

            #ifdef SPP_SANITY_CHECK
            // Assertion
            if (lru[set][way] >= ST_WAY) {
                cout << "[ST] LRU value is wrong! set: " << set << " way: " << way << " lru: " << lru[set][way] << endl;
                assert(0);
            }
            #endif
        }
    }
    lru[set][match] = 0; // Promote to the MRU position
}

void PATTERN_TABLE::update_pattern(uint32_t last_sig, int curr_delta)
{
    // Update (sig, delta) correlation
    uint32_t set = get_hash(last_sig) % PT_SET,
             match = 0;

    // Case 1: Hit
    for (match = 0; match < PT_WAY; match++) {
        if (delta[set][match] == curr_delta) {
            c_delta[set][match]++;
            c_sig[set]++;
            if (c_sig[set] > C_SIG_MAX) {
                for (uint32_t way = 0; way < PT_WAY; way++)
                    c_delta[set][way] >>= 1;
                c_sig[set] >>= 1;
            }

            SPP_DP (
                cout << "[PT] " << __func__ << " hit sig: " << hex << last_sig << dec << " set: " << set << " way: " << match;
                cout << " delta: " << delta[set][match] << " c_delta: " << c_delta[set][match] << " c_sig: " << c_sig[set] << endl;
            );

            break;
        }
    }

    // Case 2: Miss
    if (match == PT_WAY) {
        uint32_t victim_way = PT_WAY,
                 min_counter = C_SIG_MAX;

        for (match = 0; match < PT_WAY; match++) {
            if (c_delta[set][match] < min_counter) { // Select an entry with the minimum c_delta
                victim_way = match;
                min_counter = c_delta[set][match];
            }
        }

        delta[set][victim_way] = curr_delta;
        c_delta[set][victim_way] = 0;
        c_sig[set]++;
        if (c_sig[set] > C_SIG_MAX) {
            for (uint32_t way = 0; way < PT_WAY; way++)
                c_delta[set][way] >>= 1;
            c_sig[set] >>= 1;
        }

        SPP_DP (
            cout << "[PT] " << __func__ << " miss sig: " << hex << last_sig << dec << " set: " << set << " way: " << victim_way;
            cout << " delta: " << delta[set][victim_way] << " c_delta: " << c_delta[set][victim_way] << " c_sig: " << c_sig[set] << endl;
        );

        #ifdef SPP_SANITY_CHECK
        // Assertion
        if (victim_way == PT_WAY) {
            cout << "[PT] Cannot find a replacement victim!" << endl;
            assert(0);
        }
        #endif
    }
}

void PATTERN_TABLE::read_pattern(uint32_t curr_sig, int *delta_q, uint32_t *confidence_q, int32_t *perc_sum_q, uint32_t &lookahead_way, uint32_t &lookahead_conf, uint32_t &pf_q_tail, uint32_t &depth, uint64_t addr, uint64_t base_addr, uint64_t train_addr, uint64_t curr_ip, int32_t train_delta, uint32_t last_sig, uint32_t pq_occupancy, uint32_t pq_SIZE, uint32_t mshr_occupancy, uint32_t mshr_SIZE)
{
    // Update (sig, delta) correlation
    uint32_t set = get_hash(curr_sig) % PT_SET,
             local_conf = 0,
             pf_conf = 0,
             max_conf = 0;
	
	bool found_candidate = false;

    if (c_sig[set]) {
        for (uint32_t way = 0; way < PT_WAY; way++) {
            local_conf = (100 * c_delta[set][way]) / c_sig[set];
            pf_conf = depth ? (GHR.global_accuracy * c_delta[set][way] / c_sig[set] * lookahead_conf / 100) : local_conf;

			int32_t perc_sum = PERC.perc_predict(train_addr, curr_ip, GHR.ip_1, GHR.ip_2, GHR.ip_3, train_delta + delta[set][way], last_sig, curr_sig, pf_conf, depth);
			bool do_pf = (perc_sum >= PERC_THRESHOLD_LO) ? 1 : 0;
			bool fill_l2 = (perc_sum >= PERC_THRESHOLD_HI) ? 1 : 0;

			if (fill_l2 && (mshr_occupancy >= mshr_SIZE || pq_occupancy >= pq_SIZE) )
				continue;
			// Now checking against the L2C_MSHR_SIZE
			// Saving some slots in the internal PF queue by checking against do_pf
            if (pf_conf && do_pf && pf_q_tail < 100 ) {

				confidence_q[pf_q_tail] = pf_conf;
            	delta_q[pf_q_tail] = delta[set][way];
				perc_sum_q[pf_q_tail] = perc_sum;
				//cout << "WAY:  "<< way << "\tPF_CONF: " << pf_conf <<  "\tIndex: " << pf_q_tail << endl;
				SPP_DP (
					cout << "[PT] State of Features: \nTrain addr: " << train_addr << "\tCurr IP: " << curr_ip << "\tIP_1: " << GHR.ip_1 << "\tIP_2: " << GHR.ip_2 << "\tIP_3: " << GHR.ip_3 << "\tDelta: " << train_delta + delta[set][way] << "\tLastSig: " << last_sig << "\tCurrSig: " << curr_sig << "\tConf: " << pf_conf << "\tDepth: " << depth << "\tSUM: "<< perc_sum  << endl;
				);
            	// Lookahead path follows the most confident entry
            	if (pf_conf > max_conf) {
            	    lookahead_way = way;
            	    max_conf = pf_conf;
            	}
            	pf_q_tail++;
				found_candidate = true;
                
                SPP_DP (
                    cout << "[PT] " << __func__ << " HIGH CONF: " << pf_conf << " sig: " << hex << curr_sig << dec << " set: " << set << " way: " << way;
                    cout << " delta: " << delta[set][way] << " c_delta: " << c_delta[set][way] << " c_sig: " << c_sig[set];
                    cout << " conf: " << local_conf << " pf_q_tail: " << (pf_q_tail-1) << " depth: " << depth << endl;
                );
            } else {
                SPP_DP (
                    cout << "[PT] " << __func__ << "  LOW CONF: " << pf_conf << " sig: " << hex << curr_sig << dec << " set: " << set << " way: " << way;
                    cout << " delta: " << delta[set][way] << " c_delta: " << c_delta[set][way] << " c_sig: " << c_sig[set];
                    cout << " conf: " << local_conf << " pf_q_tail: " << (pf_q_tail) << " depth: " << depth << endl;
                );
            }
			// Recording Perc negatives
            if (pf_conf && pf_q_tail < L2C_MSHR_SIZE && (perc_sum < PERC_THRESHOLD_HI) ) {
				// Note: Using PERC_THRESHOLD_HI as the decising factor for negative case
				// Because 'trueness' of a prefetch is decisded based on the feedback from L2C
				// So even though LLC prefetches go through, they are treated as false wrt L2C in this case
				uint64_t pf_addr = (base_addr & ~(BLOCK_SIZE - 1)) + (delta[set][way] << LOG2_BLOCK_SIZE);
    			
                if ((addr & ~(PAGE_SIZE - 1)) == (pf_addr & ~(PAGE_SIZE - 1))) { // Prefetch request is in the same physical page
                	FILTER.check(pf_addr, train_addr, curr_ip, SPP_PERC_REJECT, train_delta + delta[set][way], last_sig, curr_sig, pf_conf, perc_sum, depth);
				}
			}
        }
        lookahead_conf = max_conf;
        if (found_candidate) depth++;

        SPP_DP (cout << "global_accuracy: " << GHR.global_accuracy << " lookahead_conf: " << lookahead_conf << endl;);
    } else confidence_q[pf_q_tail] = 0;
}

bool PREFETCH_FILTER::check(uint64_t check_addr, uint64_t base_addr, uint64_t ip, FILTER_REQUEST filter_request, int cur_delta, uint32_t last_sig, uint32_t curr_sig, uint32_t conf, int32_t sum, uint32_t depth)
{
    uint64_t cache_line = check_addr >> LOG2_BLOCK_SIZE,
             hash = get_hash(cache_line);
	
	//MAIN FILTER
	uint64_t quotient = (hash >> REMAINDER_BIT) & ((1 << QUOTIENT_BIT) - 1),
             remainder = hash % (1 << REMAINDER_BIT);
	
	//REJECT FILTER
	uint64_t quotient_reject = (hash >> REMAINDER_BIT_REJ) & ((1 << QUOTIENT_BIT_REJ) - 1),
             remainder_reject = hash % (1 << REMAINDER_BIT_REJ);

    SPP_DP (
        cout << "[FILTER] check_addr: " << hex << check_addr << " check_cache_line: " << (check_addr >> LOG2_BLOCK_SIZE);
		cout << " request type: " << filter_request;
        cout << " hash: " << hash << dec << " quotient: " << quotient << " remainder: " << remainder << endl;
    );

    switch (filter_request) {
		
		case SPP_PERC_REJECT: // To see what would have been the prediction given perceptron has rejected the PF
            if ((valid[quotient] || useful[quotient]) && remainder_tag[quotient] == remainder) { 
				// We want to check if the prefetch would have gone through had perc not rejected
				// So even in perc reject case, I'm checking in the accept filter for redundancy
                SPP_DP (
                    cout << "[FILTER] " << __func__ << " line is already in the filter check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient] << " useful: " << useful[quotient] << endl; 
                );
                return false; // False return indicates "Do not prefetch"
            } else {
				valid_reject[quotient_reject] = 1;
				remainder_tag_reject[quotient_reject] = remainder_reject;

				// Logging perc features
				address_reject[quotient_reject] = base_addr;
				pc_reject[quotient_reject] = ip;
				pc_1_reject[quotient_reject] = GHR.ip_1;
				pc_2_reject[quotient_reject] = GHR.ip_2;
				pc_3_reject[quotient_reject] = GHR.ip_3;
				delta_reject[quotient_reject] = cur_delta;
				perc_sum_reject[quotient_reject] = sum;
				last_signature_reject[quotient_reject] = last_sig;
				cur_signature_reject[quotient_reject] = curr_sig;
				confidence_reject[quotient_reject] = conf;
				la_depth_reject[quotient_reject] = depth;

				SPP_DP (
                    cout << "[FILTER] " << __func__ << " PF rejected by perceptron. Set valid_reject for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " remainder_tag: " << remainder_tag_reject[quotient_reject] << endl; 
					cout << " More Recorded Metadata: Addr: " << hex << address_reject[quotient_reject] << dec << " PC: " << pc_reject[quotient_reject] << " Delta: " << delta_reject[quotient_reject] << " Last Signature: " << last_signature_reject[quotient_reject] << " Current Signature: " << cur_signature_reject[quotient_reject] << " Confidence: " << confidence_reject[quotient_reject] << endl;
                );
			}
			break;
		
		case SPP_L2C_PREFETCH:
            if ((valid[quotient] || useful[quotient]) && remainder_tag[quotient] == remainder) { 
                SPP_DP (
                    cout << "[FILTER] " << __func__ << " line is already in the filter check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient] << " useful: " << useful[quotient] << endl; 
                );

                return false; // False return indicates "Do not prefetch"
            } else {

                valid[quotient] = 1;  // Mark as prefetched
                useful[quotient] = 0; // Reset useful bit
                remainder_tag[quotient] = remainder;

				// Logging perc features
				delta[quotient] = cur_delta;
				pc[quotient] = ip;
				pc_1[quotient] = GHR.ip_1;
				pc_2[quotient] = GHR.ip_2;
				pc_3[quotient] = GHR.ip_3;
				last_signature[quotient] = last_sig; 
				cur_signature[quotient] = curr_sig;
				confidence[quotient] = conf;
				address[quotient] = base_addr; 
				perc_sum[quotient] = sum;
				la_depth[quotient] = depth;
				
				SPP_DP (
                    cout << "[FILTER] " << __func__ << " set valid for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " remainder_tag: " << remainder_tag[quotient] << " valid: " << valid[quotient] << " useful: " << useful[quotient] << endl; 
					cout << " More Recorded Metadata: Addr:" << hex << address[quotient] << dec << " PC: " << pc[quotient] << " Delta: " << delta[quotient] << " Last Signature: " << last_signature[quotient] << " Current Signature: " << cur_signature[quotient] << " Confidence: " << confidence[quotient] << endl;
                );
            }
            break;

        case SPP_LLC_PREFETCH:
            if ((valid[quotient] || useful[quotient]) && remainder_tag[quotient] == remainder) { 
                SPP_DP (
                    cout << "[FILTER] " << __func__ << " line is already in the filter check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient] << " useful: " << useful[quotient] << endl; 
                );

                return false; // False return indicates "Do not prefetch"
            } else {
                // NOTE: SPP_LLC_PREFETCH has relatively low confidence 
                // Therefore, it is safe to prefetch this cache line in the large LLC and save precious L2C capacity
                // If this prefetch request becomes more confident and SPP eventually issues SPP_L2C_PREFETCH,
                // we can get this cache line immediately from the LLC (not from DRAM)
                // To allow this fast prefetch from LLC, SPP does not set the valid bit for SPP_LLC_PREFETCH
				
				SPP_DP (
                    cout << "[FILTER] " << __func__ << " don't set valid for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient] << " useful: " << useful[quotient] << endl; 
                );
            }
            break;

        case L2C_DEMAND:
            if ((remainder_tag[quotient] == remainder) && (useful[quotient] == 0)) {
                useful[quotient] = 1;
                if (valid[quotient]) {
					GHR.pf_useful++; // This cache line was prefetched by SPP and actually used in the program
				}

                SPP_DP (
                    cout << "[FILTER] " << __func__ << " set useful for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient] << " useful: " << useful[quotient];
                    cout << " GHR.pf_issued: " << GHR.pf_issued << " GHR.pf_useful: " << GHR.pf_useful << endl; 
					if (valid[quotient])
						cout << " Calling Perceptron Update (INC) as L2C_DEMAND was useful" << endl;
                );

                if (valid[quotient]) {
					// Prefetch leads to a demand hit
					PERC.perc_update(address[quotient], pc[quotient], pc_1[quotient], pc_2[quotient], pc_3[quotient], delta[quotient], last_signature[quotient], cur_signature[quotient], confidence[quotient], la_depth[quotient], 1, perc_sum[quotient]);
				}
            }
			//If NOT Prefetched
			if (!(valid[quotient] && remainder_tag[quotient] == remainder)) {
				// AND If Rejected by Perc
				if (valid_reject[quotient_reject] && remainder_tag_reject[quotient_reject] == remainder_reject) {
               		 SPP_DP (
               		     cout << "[FILTER] " << __func__ << " not doing anything for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
               		     cout << " quotient: " << quotient << " valid_reject:" << valid_reject[quotient_reject];
               		     cout << " GHR.pf_issued: " << GHR.pf_issued << " GHR.pf_useful: " << GHR.pf_useful << endl; 
			   		 	 cout << " Calling Perceptron Update (DEC) as a useful L2C_DEMAND was rejected and reseting valid_reject" << endl;
               		 );
					// Not prefetched but could have been a good idea to prefetch
					PERC.perc_update(address_reject[quotient_reject], pc_reject[quotient_reject], pc_1_reject[quotient_reject], pc_2_reject[quotient_reject], pc_3_reject[quotient_reject], delta_reject[quotient_reject], last_signature_reject[quotient_reject], cur_signature_reject[quotient_reject], confidence_reject[quotient_reject], la_depth_reject[quotient_reject], 0, perc_sum_reject[quotient_reject]);
					valid_reject[quotient_reject] = 0;
					remainder_tag_reject[quotient_reject] = 0;
				}
			}
            break;

        case L2C_EVICT:
            // Decrease global pf_useful counter when there is a useless prefetch (prefetched but not used)
            if (valid[quotient] && !useful[quotient]) {
				if (GHR.pf_useful) 
					GHR.pf_useful--;
				
				SPP_DP (
               		cout << "[FILTER] " << __func__ << " eviction for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
               		cout << " quotient: " << quotient << " valid: " << valid[quotient] << " useful: " << useful[quotient] << endl; 
					cout << " Calling Perceptron Update (DEC) as L2C_DEMAND was not useful" << endl;
					cout << " Reseting valid_reject" << endl;
            	);

				// Prefetch leads to eviction
				PERC.perc_update(address[quotient], pc[quotient], pc_1[quotient], pc_2[quotient], pc_3[quotient], delta[quotient], last_signature[quotient], cur_signature[quotient], confidence[quotient], la_depth[quotient], 0, perc_sum[quotient]);
			}
            // Reset filter entry
            valid[quotient] = 0;
            useful[quotient] = 0;
            remainder_tag[quotient] = 0;

			// Reset reject filter too
			valid_reject[quotient_reject] = 0;
			remainder_tag_reject[quotient_reject] = 0;

            break;

        default:
            // Assertion
            cout << "[FILTER] Invalid filter request type: " << filter_request << endl;
            assert(0);
    }

    return true;
}

void GLOBAL_REGISTER::update_entry(uint32_t pf_sig, uint32_t pf_confidence, uint32_t pf_offset, int pf_delta) 
{
    // NOTE: GHR implementation is slightly different from the original paper
    // Instead of matching (last_offset + delta), GHR simply stores and matches the pf_offset
    uint32_t min_conf = 100,
             victim_way = MAX_GHR_ENTRY;

    SPP_DP (
        cout << "[GHR] Crossing the page boundary pf_sig: " << hex << pf_sig << dec;
        cout << " confidence: " << pf_confidence << " pf_offset: " << pf_offset << " pf_delta: " << pf_delta << endl;
    );

    for (uint32_t i = 0; i < MAX_GHR_ENTRY; i++) {
        //if (sig[i] == pf_sig) { // TODO: Which one is better and consistent?
            // If GHR already holds the same pf_sig, update the GHR entry with the latest info
        if (valid[i] && (offset[i] == pf_offset)) {
            // If GHR already holds the same pf_offset, update the GHR entry with the latest info
            sig[i] = pf_sig;
            confidence[i] = pf_confidence;
            //offset[i] = pf_offset;
            delta[i] = pf_delta;

            SPP_DP (cout << "[GHR] Found a matching index: " << i << endl;);

            return;
        }

        // GHR replacement policy is based on the stored confidence value
        // An entry with the lowest confidence is selected as a victim
        if (confidence[i] < min_conf) {
            min_conf = confidence[i];
            victim_way = i;
        }
    }

    // Assertion
    if (victim_way >= MAX_GHR_ENTRY) {
        cout << "[GHR] Cannot find a replacement victim!" << endl;
        assert(0);
    }

    SPP_DP (
        cout << "[GHR] Replace index: " << victim_way << " pf_sig: " << hex << sig[victim_way] << dec;
        cout << " confidence: " << confidence[victim_way] << " pf_offset: " << offset[victim_way] << " pf_delta: " << delta[victim_way] << endl;
    );

    valid[victim_way] = 1;
    sig[victim_way] = pf_sig;
    confidence[victim_way] = pf_confidence;
    offset[victim_way] = pf_offset;
    delta[victim_way] = pf_delta;
}

uint32_t GLOBAL_REGISTER::check_entry(uint32_t page_offset)
{
    uint32_t max_conf = 0,
             max_conf_way = MAX_GHR_ENTRY;

    for (uint32_t i = 0; i < MAX_GHR_ENTRY; i++) {
        if ((offset[i] == page_offset) && (max_conf < confidence[i])) {
            max_conf = confidence[i];
            max_conf_way = i;
        }
    }

    return max_conf_way;
}

void get_perc_index(uint64_t base_addr, uint64_t ip, uint64_t ip_1, uint64_t ip_2, uint64_t ip_3, int32_t cur_delta, uint32_t last_sig, uint32_t curr_sig, uint32_t confidence, uint32_t depth, uint64_t perc_set[PERC_FEATURES])
{
	// Returns the imdexes for the perceptron tables
    uint64_t cache_line = base_addr >> LOG2_BLOCK_SIZE,
			 page_addr  = base_addr >> LOG2_PAGE_SIZE;

	int sig_delta = (cur_delta < 0) ? (((-1) * cur_delta) + (1 << (SIG_DELTA_BIT - 1))) : cur_delta;
	uint64_t  pre_hash[PERC_FEATURES];

	pre_hash[0] = base_addr;
	pre_hash[1] = cache_line;
	pre_hash[2] = page_addr;
	pre_hash[3] = confidence ^ page_addr;
	pre_hash[4] = curr_sig ^ sig_delta;
	pre_hash[5] = ip_1 ^ (ip_2>>1) ^ (ip_3>>2);
	pre_hash[6] = ip ^ depth;
	pre_hash[7] = ip ^ sig_delta;
	pre_hash[8] = confidence;

	for (int i = 0; i < PERC_FEATURES; i++) {
		perc_set[i] = (pre_hash[i]) % PERC.PERC_DEPTH[i]; // Variable depths
		SPP_DP (
			cout << "  Perceptron Set Index#: " << i << " = " <<  perc_set[i];
		);
	}
	SPP_DP (
		cout << endl;
	);		
}

int32_t	PERCEPTRON::perc_predict(uint64_t base_addr, uint64_t ip, uint64_t ip_1, uint64_t ip_2, uint64_t ip_3, int32_t cur_delta, uint32_t last_sig, uint32_t curr_sig, uint32_t confidence, uint32_t depth)
{
	SPP_DP (
		int sig_delta = (cur_delta < 0) ? (((-1) * cur_delta) + (1 << (SIG_DELTA_BIT - 1))) : cur_delta;
		cout << "[PERC_PRED] Current IP: " << ip << "  and  Memory Adress: " << hex << base_addr << endl;
		cout << " Last Sig: " << last_sig << " Curr Sig: " << curr_sig << dec << endl;
		cout << " Cur Delta: " << cur_delta << " Sign Delta: " << sig_delta << " Confidence: " << confidence<< endl;
		cout << " ";
	);

	uint64_t perc_set[PERC_FEATURES];
	// Get the indexes in perc_set[]
	get_perc_index(base_addr, ip, ip_1, ip_2, ip_3, cur_delta, last_sig, curr_sig, confidence, depth, perc_set);
	
	int32_t sum = 0;
	for (int i = 0; i < PERC_FEATURES; i++) {
		sum += perc_weights[perc_set[i]][i];	
		// Calculate Sum
	}
	SPP_DP (
		cout << " Sum of perceptrons: " << sum << " Prediction made: " << ((sum >= PERC_THRESHOLD_LO) ?  ((sum >= PERC_THRESHOLD_HI) ? FILL_L2 : FILL_LLC) : 0)  << endl;
	);
	// Return the sum
	return sum;
}

void 	PERCEPTRON::perc_update(uint64_t base_addr, uint64_t ip, uint64_t ip_1, uint64_t ip_2, uint64_t ip_3, int32_t cur_delta, uint32_t last_sig, uint32_t curr_sig, uint32_t confidence, uint32_t depth, bool direction, int32_t perc_sum)
{
	SPP_DP (
		int sig_delta = (cur_delta < 0) ? (((-1) * cur_delta) + (1 << (SIG_DELTA_BIT - 1))) : cur_delta;
		cout << "[PERC_UPD] (Recorded) IP: " << ip << "  and  Memory Adress: " << hex << base_addr << endl;
		cout << " Last Sig: " << last_sig << " Curr Sig: " << curr_sig << dec << endl;
		cout << " Cur Delta: " << cur_delta << " Sign Delta: " << sig_delta << " Confidence: "<< confidence << " Update Direction: " << direction << endl;
		cout << " ";
	);

	uint64_t perc_set[PERC_FEATURES];
	// Get the perceptron indexes
	get_perc_index(base_addr, ip, ip_1, ip_2, ip_3, cur_delta, last_sig, curr_sig, confidence, depth, perc_set);
	
	int32_t sum = 0;
	// Restore the sum that led to the prediction
	sum = perc_sum;
	
	if (!direction) { // direction = 1 means the sum was in the correct direction, 0 means it was in the wrong direction
		// Prediction wrong
		for (int i = 0; i < PERC_FEATURES; i++) {
			if (sum >= PERC_THRESHOLD_HI) {
				// Prediction was to prefectch -- so decrement counters
				if (perc_weights[perc_set[i]][i] > -1*(PERC_COUNTER_MAX+1) )
					perc_weights[perc_set[i]][i]--;
			}
			if (sum < PERC_THRESHOLD_HI) {
				// Prediction was to not prefetch -- so increment counters
				if (perc_weights[perc_set[i]][i] < PERC_COUNTER_MAX)
					perc_weights[perc_set[i]][i]++;
			}
		}
		SPP_DP (
			int differential = (sum >= PERC_THRESHOLD_HI) ? -1 : 1;
			cout << " Direction is: " << direction << " and sum is:" << sum;
			cout << " Overall Differential: " << differential << endl;
		);
	}
	if (direction && sum > NEG_UPDT_THRESHOLD && sum < POS_UPDT_THRESHOLD) {
		// Prediction correct but sum not 'saturated' enough
		for (int i = 0; i < PERC_FEATURES; i++) {
			if (sum >= PERC_THRESHOLD_HI) {
				// Prediction was to prefetch -- so increment counters
				if (perc_weights[perc_set[i]][i] < PERC_COUNTER_MAX)
					perc_weights[perc_set[i]][i]++;
			}
			if (sum < PERC_THRESHOLD_HI) {
				// Prediction was to not prefetch -- so decrement counters
				if (perc_weights[perc_set[i]][i] > -1*(PERC_COUNTER_MAX+1) )
					perc_weights[perc_set[i]][i]--;
			}
		}
		SPP_DP (
			int differential = 0;
			if (sum >= PERC_THRESHOLD_HI) differential =  1;
			if (sum  < PERC_THRESHOLD_HI) differential = -1;
			cout << " Direction is: " << direction << " and sum is:" << sum;
			cout << " Overall Differential: " << differential << endl;
		);
	}
}
