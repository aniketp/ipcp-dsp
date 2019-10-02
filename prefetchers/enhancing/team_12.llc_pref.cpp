// Submission For DPC-3 Team 12
// LLC Prefetcher is next-line in single-core and SPP in multi-core

#include "cache.h"
#include <algorithm>
#include <cmath>

/*
 * Header for LLC SPP
 */
#ifndef SPP_H
#define SPP_H

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
#define PF_THRESHOLD_LLC 75

// Global register parameters
#define GLOBAL_COUNTER_BIT 10
#define GLOBAL_COUNTER_MAX ((1 << GLOBAL_COUNTER_BIT) - 1) 
#define MAX_GHR_ENTRY 8
#define PAGES_TRACKED 8

enum FILTER_REQUEST {SPP_LLC_PREFETCH, LLC_DEMAND, LLC_EVICT}; // Request type for prefetch filter
uint64_t get_hash_l(uint64_t key);

class SIGNATURE_TABLE_L {
  public:
    bool     valid[ST_SET][ST_WAY][NUM_CPUS];
    uint32_t tag[ST_SET][ST_WAY][NUM_CPUS],
             last_offset[ST_SET][ST_WAY][NUM_CPUS],
             sig[ST_SET][ST_WAY][NUM_CPUS],
             lru[ST_SET][ST_WAY][NUM_CPUS];

    SIGNATURE_TABLE_L() {
        cout << "Initialize SIGNATURE TABLE" << endl;
        cout << "ST_SET: " << ST_SET << endl;
        cout << "ST_WAY: " << ST_WAY << endl;
        cout << "ST_TAG_BIT: " << ST_TAG_BIT << endl;
        cout << "ST_TAG_MASK: " << hex << ST_TAG_MASK << dec << endl;

        	for (uint32_t set = 0; set < ST_SET; set++){
            	for (uint32_t way = 0; way < ST_WAY; way++) {
        			for (uint32_t i = 0; i < NUM_CPUS; i++) {
            	    	valid[set][way][i] = 0;
                		tag[set][way][i] = 0;
                		last_offset[set][way][i] = 0;
                		sig[set][way][i] = 0;
                		lru[set][way][i] = way;
					}
            	}
			}
    };

    void read_and_update_sig(uint64_t page, uint32_t page_offset, uint32_t &last_sig, uint32_t &curr_sig, int32_t &delta, int cpu);
};

class PATTERN_TABLE_L {
  public:
    int      delta[PT_SET][PT_WAY][NUM_CPUS];
    uint32_t c_delta[PT_SET][PT_WAY][NUM_CPUS],
             c_sig[PT_SET][NUM_CPUS];

    PATTERN_TABLE_L() {
        cout << endl << "Initialize PATTERN TABLE" << endl;
        cout << "PT_SET: " << PT_SET << endl;
        cout << "PT_WAY: " << PT_WAY << endl;
        cout << "SIG_DELTA_BIT: " << SIG_DELTA_BIT << endl;
        cout << "C_SIG_BIT: " << C_SIG_BIT << endl;
        cout << "C_DELTA_BIT: " << C_DELTA_BIT << endl;

        for (uint32_t i = 0; i < NUM_CPUS; i++) {
        	for (uint32_t set = 0; set < PT_SET; set++) {
            	for (uint32_t way = 0; way < PT_WAY; way++) {
            	    delta[set][way][i] = 0;
            	    c_delta[set][way][i] = 0;
            	}
            	c_sig[set][i] = 0;
        	}
		}
    }

    void update_pattern(uint32_t last_sig, int curr_delta, int cpu),
         read_pattern(uint32_t curr_sig, int *prefetch_delta, uint32_t *confidence_q, uint32_t &lookahead_way, uint32_t &lookahead_conf, uint32_t &pf_q_tail, uint32_t &depth, uint32_t pq_occupancy, uint32_t pq_SIZE, uint32_t mshr_occupancy, uint32_t mshr_SIZE, int cpu);
};

class PREFETCH_FILTER_L {
  public:
    uint64_t remainder_tag[FILTER_SET][NUM_CPUS];
    bool     valid[FILTER_SET][NUM_CPUS],  // Consider this as "prefetched"
             useful[FILTER_SET][NUM_CPUS]; // Consider this as "used"

    PREFETCH_FILTER_L() {
        cout << endl << "Initialize PREFETCH FILTER" << endl;
        cout << "FILTER_SET: " << FILTER_SET << endl;

        for (uint32_t i = 0; i < NUM_CPUS; i++) {
        	for (uint32_t set = 0; set < FILTER_SET; set++) {
        	    remainder_tag[set][i] = 0;
        	    valid[set][i] = 0;
        	    useful[set][i] = 0;
        	}
		}

    }

    bool     check(uint64_t pf_addr, FILTER_REQUEST filter_request, int cpu);
};

class GLOBAL_REGISTER_L {
  public:
    // Global counters to calculate global prefetching accuracy
	
	uint64_t page_tracker[PAGES_TRACKED][NUM_CPUS];

    uint64_t pf_useful[NUM_CPUS],
             pf_issued[NUM_CPUS],
             global_accuracy[NUM_CPUS]; // Alpha value in Section III. Equation 3

	// Global History Register (GHR) entries
    uint8_t  valid[MAX_GHR_ENTRY][NUM_CPUS];
    uint32_t sig[MAX_GHR_ENTRY][NUM_CPUS],
             confidence[MAX_GHR_ENTRY][NUM_CPUS],
             offset[MAX_GHR_ENTRY][NUM_CPUS];
    int      delta[MAX_GHR_ENTRY][NUM_CPUS];

    GLOBAL_REGISTER_L() {
		for (int j = 0; j < NUM_CPUS; j++) {
        	pf_useful[j] = 0;
        	pf_issued[j] = 0;
        	global_accuracy[j] = 0;

			for (uint32_t i = 0; i < MAX_GHR_ENTRY; i++) {
        	    		valid[i][j] = 0;
        	              sig[i][j] = 0;
        	       confidence[i][j] = 0;
        	    	   offset[i][j] = 0;
        	    		delta[i][j] = 0;
        	}
		}
    }

    void update_entry(uint32_t pf_sig, uint32_t pf_confidence, uint32_t pf_offset, int pf_delta, int cpu);
    uint32_t check_entry(uint32_t page_offset, int cpu);
};

#endif
/*
 * Header for LLC SPP ends
 */


SIGNATURE_TABLE_L ST_L;
PATTERN_TABLE_L   PT_L;
PREFETCH_FILTER_L FILTER_L;
GLOBAL_REGISTER_L GHR_L;


void CACHE::llc_prefetcher_initialize() 
{

}

uint32_t CACHE::llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{

if (NUM_CPUS == 1) { // Single-Core Config
	// Ignoring any prefetch request which is not the last suggestion generated by L1D prefetcher
	if (type == PREFETCH && metadata_in !=6)
		return metadata_in;

  	uint64_t pf_addr = ((addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;
  	prefetch_line(ip, addr, pf_addr, FILL_LLC, 7);

  	return metadata_in;
}

if (NUM_CPUS == 4) { // Multi-core config
	
	if (type == PREFETCH)
		return metadata_in;

    uint64_t page = addr >> LOG2_PAGE_SIZE;
	uint32_t page_offset = (addr >> LOG2_BLOCK_SIZE) & (PAGE_SIZE / BLOCK_SIZE - 1),
             last_sig = 0,
             curr_sig = 0,
             confidence_q[LLC_MSHR_SIZE],
             depth = 0;

    int32_t  delta = 0,
             delta_q[LLC_MSHR_SIZE];

    for (uint32_t i = 0; i < LLC_MSHR_SIZE; i++){
        confidence_q[i] = 0;
        delta_q[i] = 0;
    }

#ifdef GHR_ON
	for (int i = PAGES_TRACKED-1; i>0; i--) { // N down to 1
		GHR_L.page_tracker[i][cpu] = GHR_L.page_tracker[i-1][cpu];
	}
	GHR_L.page_tracker[0][cpu] = page;

	int distinct_pages = 0;
	for (int i=0; i < PAGES_TRACKED; i++) {
		int j;
		for (j=0; j<i; j++) {
			if (GHR_L.page_tracker[i][cpu] == GHR_L.page_tracker[j][cpu])
				break;
		}
		if (i==j)
			distinct_pages++;
	}
	//cout << "Distinct Pages: " << distinct_pages << endl;
#else
	int distinct_pages = 1;
#endif

    confidence_q[0] = 100;
    GHR_L.global_accuracy[cpu] = GHR_L.pf_issued[cpu] ? ((100 * GHR_L.pf_useful[cpu]) / GHR_L.pf_issued[cpu])  : 0;
    
    SPP_DP (
        cout << endl << "[ChampSim] " << __func__ << " addr: " << hex << addr << " cache_line: " << (addr >> LOG2_BLOCK_SIZE);
        cout << " page: " << page << " page_offset: " << dec << page_offset << endl;
    );

    // Stage 1: Read and update a sig stored in ST
    // last_sig and delta are used to update (sig, delta) correlation in PT
    // curr_sig is used to read prefetch candidates in PT 
    ST_L.read_and_update_sig(page, page_offset, last_sig, curr_sig, delta, cpu);

    // Also check the prefetch filter in parallel to update global accuracy counters 
    FILTER_L.check(addr, LLC_DEMAND, cpu); 

    // Stage 2: Update delta patterns stored in PT
    if (last_sig) PT_L.update_pattern(last_sig, delta, cpu);

    uint32_t lookahead_conf = 100;
    // Stage 3: Start prefetching

    uint64_t base_addr = addr;
    uint32_t pf_q_head = 0, 
             pf_q_tail = 0;
    uint8_t  do_lookahead = 0;

	uint8_t num_pf = 0; 
#ifdef LOOKAHEAD_ON
    do {
#endif
        uint32_t lookahead_way = PT_WAY;
        PT_L.read_pattern(curr_sig, delta_q, confidence_q, lookahead_way, lookahead_conf, pf_q_tail, depth, PQ.occupancy, PQ.SIZE, MSHR.occupancy, MSHR.SIZE, cpu);

        do_lookahead = 0;
        for (uint32_t i = pf_q_head; i < pf_q_tail; i++) {
            if (confidence_q[i] >= PF_THRESHOLD_LLC) {
                uint64_t pf_addr = (base_addr & ~(BLOCK_SIZE - 1)) + (delta_q[i] << LOG2_BLOCK_SIZE);

                if ((addr & ~(PAGE_SIZE - 1)) == (pf_addr & ~(PAGE_SIZE - 1))) { // Prefetch request is in the same physical page
					if ( num_pf <((MSHR.SIZE)/distinct_pages)+1 ) {
                    	if (FILTER_L.check(pf_addr, SPP_LLC_PREFETCH, cpu)) {
							prefetch_line(ip, addr, pf_addr, FILL_LLC, 0); // Use addr (not base_addr) to obey the same physical page boundary
							num_pf++;

                        	if (confidence_q[i] >= PF_THRESHOLD_LLC) {
                            	GHR_L.pf_issued[cpu]++;
                            	if (GHR_L.pf_issued[cpu] > GLOBAL_COUNTER_MAX) {
                                	GHR_L.pf_issued[cpu] >>= 1;
                                	GHR_L.pf_useful[cpu] >>= 1;
                            	}
                            	SPP_DP (cout << "[ChampSim] SPP L2 prefetch issued GHR_L.pf_issued: " << GHR_L.pf_issued[cpu] << " GHR_L.pf_useful: " << GHR_L.pf_useful[cpu] << endl;);
                        	}

                        	SPP_DP (
                            	cout << "[ChampSim] " << __func__ << " base_addr: " << hex << base_addr << " pf_addr: " << pf_addr;
                            	cout << " pf_cache_line: " << (pf_addr >> LOG2_BLOCK_SIZE);
                            	cout << " prefetch_delta: " << dec << delta_q[i] << " confidence: " << confidence_q[i];
                            	cout << " depth: " << i << " fill_level: " << FILL_LLC << endl;
                        	);
						}
                    }
                } else { // Prefetch request is crossing the physical page boundary
#ifdef GHR_ON
                    // Store this prefetch request in GHR to bootstrap SPP learning when we see a ST miss (i.e., accessing a new page)
                    GHR_L.update_entry(curr_sig, confidence_q[i], (pf_addr >> LOG2_BLOCK_SIZE) & 0x3F, delta_q[i], cpu); 
#endif
                }

                do_lookahead = 1;
                pf_q_head++;
            }
        }

        // Update base_addr and curr_sig
        if (lookahead_way < PT_WAY) {
            uint32_t set = get_hash_l(curr_sig) % PT_SET;
            base_addr += (PT_L.delta[set][lookahead_way][cpu] << LOG2_BLOCK_SIZE);

            // PT.delta uses a 7-bit sign magnitude representation to generate sig_delta
            int sig_delta = (PT_L.delta[set][lookahead_way][cpu] < 0) ? (((-1) * PT_L.delta[set][lookahead_way][cpu]) + (1 << (SIG_DELTA_BIT - 1))) : PT_L.delta[set][lookahead_way][cpu];
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
}

uint32_t CACHE::llc_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{

if (NUM_CPUS == 1) {
  	return metadata_in;
}

if (NUM_CPUS == 4) {
#ifdef FILTER_ON
    SPP_DP (cout << endl;);
    FILTER_L.check(evicted_addr, LLC_EVICT, cpu);
#endif

    return metadata_in;
}

}

void CACHE::llc_prefetcher_final_stats()
{

}

/*
 * Below are SPP Specific Functions
 */

uint64_t get_hash_l(uint64_t key)
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

void SIGNATURE_TABLE_L::read_and_update_sig(uint64_t page, uint32_t page_offset, uint32_t &last_sig, uint32_t &curr_sig, int32_t &delta, int cpu)
{
    uint32_t set = get_hash_l(page) % ST_SET,
             match = ST_WAY,
             partial_page = page & ST_TAG_MASK;
    uint8_t  ST_hit = 0;
    int      sig_delta = 0;

    SPP_DP (cout << "[ST] " << __func__ << " page: " << hex << page << " partial_page: " << partial_page << dec << endl;);

    // Case 1: Hit
    for (match = 0; match < ST_WAY; match++) {
        if (valid[set][match][cpu] && (tag[set][match][cpu] == partial_page)) {
            last_sig = sig[set][match][cpu];
            delta = page_offset - last_offset[set][match][cpu];

            if (delta) {
                // Build a new sig based on 7-bit sign magnitude representation of delta
                //sig_delta = (delta < 0) ? ((((-1) * delta) & 0x3F) + 0x40) : delta;
                sig_delta = (delta < 0) ? (((-1) * delta) + (1 << (SIG_DELTA_BIT - 1))) : delta;
                sig[set][match][cpu] = ((last_sig << SIG_SHIFT) ^ sig_delta) & SIG_MASK;
                curr_sig = sig[set][match][cpu];
                last_offset[set][match][cpu] = page_offset;

                SPP_DP (
                    cout << "[ST] " << __func__ << " hit set: " << set << " way: " << match;
                    cout << " valid: " << valid[set][match][cpu] << " tag: " << hex << tag[set][match][cpu];
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
            if (valid[set][match][cpu] == 0) {
                valid[set][match][cpu] = 1;
                tag[set][match][cpu] = partial_page;
                sig[set][match][cpu] = 0;
                curr_sig = sig[set][match][cpu];
                last_offset[set][match][cpu] = page_offset;

                SPP_DP (
                    cout << "[ST] " << __func__ << " invalid set: " << set << " way: " << match;
                    cout << " valid: " << valid[set][match][cpu] << " tag: " << hex << partial_page;
                    cout << " sig: " << sig[set][match][cpu] << " last_offset: " << dec << page_offset << endl;
                );

                break;
            }
        }
    }

    // Case 3: Miss
    if (match == ST_WAY) {
        for (match = 0; match < ST_WAY; match++) {
            if (lru[set][match][cpu] == ST_WAY - 1) { // Find replacement victim
                tag[set][match][cpu] = partial_page;
                sig[set][match][cpu] = 0;
                curr_sig = sig[set][match][cpu];
                last_offset[set][match][cpu] = page_offset;

                SPP_DP (
                    cout << "[ST] " << __func__ << " miss set: " << set << " way: " << match;
                    cout << " valid: " << valid[set][match][cpu] << " victim tag: " << hex << tag[set][match][cpu] << " new tag: " << partial_page;
                    cout << " sig: " << sig[set][match][cpu] << " last_offset: " << dec << page_offset << endl;
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
        uint32_t GHR_found = GHR_L.check_entry(page_offset, cpu);
        if (GHR_found < MAX_GHR_ENTRY) {
            sig_delta = (GHR_L.delta[GHR_found][cpu] < 0) ? (((-1) * GHR_L.delta[GHR_found][cpu]) + (1 << (SIG_DELTA_BIT - 1))) : GHR_L.delta[GHR_found][cpu];
            sig[set][match][cpu] = ((GHR_L.sig[GHR_found][cpu] << SIG_SHIFT) ^ sig_delta) & SIG_MASK;
            curr_sig = sig[set][match][cpu];
        }
    }
#endif

    // Update LRU
    for (uint32_t way = 0; way < ST_WAY; way++) {
        if (lru[set][way][cpu] < lru[set][match][cpu]) {
            lru[set][way][cpu]++;

            #ifdef SPP_SANITY_CHECK
            // Assertion
            if (lru[set][way][cpu] >= ST_WAY) {
                cout << "[ST] LRU value is wrong! set: " << set << " way: " << way << " lru: " << lru[set][way] << endl;
                assert(0);
            }
            #endif
        }
    }
    lru[set][match][cpu] = 0; // Promote to the MRU position
}

void PATTERN_TABLE_L::update_pattern(uint32_t last_sig, int curr_delta, int cpu)
{
    // Update (sig, delta) correlation
    uint32_t set = get_hash_l(last_sig) % PT_SET,
             match = 0;

    // Case 1: Hit
    for (match = 0; match < PT_WAY; match++) {
        if (delta[set][match][cpu] == curr_delta) {
            c_delta[set][match][cpu]++;
            c_sig[set][cpu]++;
            if (c_sig[set][cpu] > C_SIG_MAX) {
                for (uint32_t way = 0; way < PT_WAY; way++)
                    c_delta[set][way][cpu] >>= 1;
                c_sig[set][cpu] >>= 1;
            }

            SPP_DP (
                cout << "[PT] " << __func__ << " hit sig: " << hex << last_sig << dec << " set: " << set << " way: " << match;
                cout << " delta: " << delta[set][match][cpu] << " c_delta: " << c_delta[set][match][cpu] << " c_sig: " << c_sig[set][cpu] << endl;
            );

            break;
        }
    }

    // Case 2: Miss
    if (match == PT_WAY) {
        uint32_t victim_way = PT_WAY,
                 min_counter = C_SIG_MAX;

        for (match = 0; match < PT_WAY; match++) {
            if (c_delta[set][match][cpu] < min_counter) { // Select an entry with the minimum c_delta
                victim_way = match;
                min_counter = c_delta[set][match][cpu];
            }
        }

        delta[set][victim_way][cpu] = curr_delta;
        c_delta[set][victim_way][cpu] = 0;
        c_sig[set][cpu]++;
        if (c_sig[set][cpu] > C_SIG_MAX) {
            for (uint32_t way = 0; way < PT_WAY; way++)
                c_delta[set][way][cpu] >>= 1;
            c_sig[set][cpu] >>= 1;
        }

        SPP_DP (
            cout << "[PT] " << __func__ << " miss sig: " << hex << last_sig << dec << " set: " << set << " way: " << victim_way;
            cout << " delta: " << delta[set][victim_way][cpu] << " c_delta: " << c_delta[set][victim_way][cpu] << " c_sig: " << c_sig[set][cpu] << endl;
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

void PATTERN_TABLE_L::read_pattern(uint32_t curr_sig, int *delta_q, uint32_t *confidence_q, uint32_t &lookahead_way, uint32_t &lookahead_conf, uint32_t &pf_q_tail, uint32_t &depth, uint32_t pq_occupancy, uint32_t pq_SIZE, uint32_t mshr_occupancy, uint32_t mshr_SIZE, int cpu)
{
    // Update (sig, delta) correlation
    uint32_t set = get_hash_l(curr_sig) % PT_SET,
             local_conf = 0,
             pf_conf = 0,
             max_conf = 0;

    if (c_sig[set][cpu]) {
        for (uint32_t way = 0; way < PT_WAY; way++) {
            local_conf = (100 * c_delta[set][way][cpu]) / c_sig[set][cpu];

            pf_conf = depth ? (GHR_L.global_accuracy[cpu] * c_delta[set][way][cpu] / c_sig[set][cpu] * lookahead_conf / 100) : local_conf;
			// debug
			if (pf_conf > 100) {
				cout << "[PT] " << __func__ << " CONF ERROR!! c_delta: " << c_delta[set][way][cpu] << " c_sig: " << c_sig[set][cpu] << " Alpha: " << GHR_L.global_accuracy[cpu] << endl;
			}

            if (pf_conf >= PF_THRESHOLD_LLC && pq_occupancy < pq_SIZE && mshr_occupancy < mshr_SIZE && pf_q_tail < LLC_MSHR_SIZE) {
                confidence_q[pf_q_tail] = pf_conf;
                delta_q[pf_q_tail] = delta[set][way][cpu];

                // Lookahead path follows the most confident entry
                if (pf_conf > max_conf) {
                    lookahead_way = way;
                    max_conf = pf_conf;
                }
                pf_q_tail++;

                SPP_DP (
                    cout << "[PT] " << __func__ << " HIGH CONF: " << pf_conf << " sig: " << hex << curr_sig << dec << " set: " << set << " way: " << way;
                    cout << " delta: " << delta[set][way][cpu] << " c_delta: " << c_delta[set][way][cpu] << " c_sig: " << c_sig[set][cpu];
                    cout << " conf: " << local_conf << " depth: " << depth << endl;
                );
            } else {
                SPP_DP (
                    cout << "[PT] " << __func__ << "  LOW CONF: " << pf_conf << " sig: " << hex << curr_sig << dec << " set: " << set << " way: " << way;
                    cout << " delta: " << delta[set][way][cpu] << " c_delta: " << c_delta[set][way][cpu] << " c_sig: " << c_sig[set][cpu];
                    cout << " conf: " << local_conf << " depth: " << depth << endl;
                );
            }
        }
        lookahead_conf = max_conf;
        if (lookahead_conf >= PF_THRESHOLD_LLC) depth++;

        SPP_DP (cout << "global_accuracy: " << GHR_L.global_accuracy[cpu] << " lookahead_conf: " << lookahead_conf << endl;);
    } else confidence_q[pf_q_tail] = 0;
}

bool PREFETCH_FILTER_L::check(uint64_t check_addr, FILTER_REQUEST filter_request, int cpu)
{
    uint64_t cache_line = check_addr >> LOG2_BLOCK_SIZE,
             hash = get_hash_l(cache_line),
             quotient = (hash >> REMAINDER_BIT) & ((1 << QUOTIENT_BIT) - 1),
             remainder = hash % (1 << REMAINDER_BIT);

    SPP_DP (
        cout << "[FILTER] check_addr: " << hex << check_addr << " check_cache_line: " << (check_addr >> LOG2_BLOCK_SIZE);
        cout << " hash: " << hash << dec << " quotient: " << quotient << " remainder: " << remainder << endl;
    );

    switch (filter_request) {
        case SPP_LLC_PREFETCH:
            if ((valid[quotient][cpu] || useful[quotient][cpu]) && remainder_tag[quotient][cpu] == remainder) { 
                SPP_DP (
                    cout << "[FILTER] " << __func__ << " line is already in the filter check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient][cpu] << " useful: " << useful[quotient][cpu] << endl; 
                );

                return false; // False return indicates "Do not prefetch"
            } else {
                valid[quotient][cpu] = 1;  // Mark as prefetched
                useful[quotient][cpu] = 0; // Reset useful bit
                remainder_tag[quotient][cpu] = remainder;

                SPP_DP (
                    cout << "[FILTER] " << __func__ << " set valid for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " remainder_tag: " << remainder_tag[quotient][cpu] << " valid: " << valid[quotient][cpu] << " useful: " << useful[quotient][cpu] << endl; 
                );
            }
            break;

        case LLC_DEMAND:
            if ((remainder_tag[quotient][cpu] == remainder) && (useful[quotient][cpu] == 0)) {
                useful[quotient][cpu] = 1;
                if (valid[quotient][cpu]) GHR_L.pf_useful[cpu]++; // This cache line was prefetched by SPP and actually used in the program

                SPP_DP (
                    cout << "[FILTER] " << __func__ << " set useful for check_addr: " << hex << check_addr << " cache_line: " << cache_line << dec;
                    cout << " quotient: " << quotient << " valid: " << valid[quotient][cpu] << " useful: " << useful[quotient][cpu];
                    cout << " GHR_L.pf_issued: " << GHR_L.pf_issued[cpu] << " GHR_L.pf_useful: " << GHR_L.pf_useful[cpu] << endl; 
                );
            }
            break;

        case LLC_EVICT:
            // Decrease global pf_useful counter when there is a useless prefetch (prefetched but not used)
            if (valid[quotient][cpu] && !useful[quotient][cpu] && GHR_L.pf_useful[cpu]) GHR_L.pf_useful[cpu]--;

            // Reset filter entry
            valid[quotient][cpu] = 0;
            useful[quotient][cpu] = 0;
            remainder_tag[quotient][cpu] = 0;
            break;

        default:
            // Assertion
            cout << "[FILTER] Invalid filter request type: " << filter_request << endl;
            assert(0);
    }

    return true;
}

void GLOBAL_REGISTER_L::update_entry(uint32_t pf_sig, uint32_t pf_confidence, uint32_t pf_offset, int pf_delta, int cpu) 
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
        if (valid[i][cpu] && (offset[i][cpu] == pf_offset)) {
            sig[i][cpu] = pf_sig;
            confidence[i][cpu] = pf_confidence;
            delta[i][cpu] = pf_delta;

            SPP_DP (cout << "[GHR] Found a matching index: " << i << endl;);

            return;
        }

        if (confidence[i][cpu] < min_conf) {
            min_conf = confidence[i][cpu];
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
        cout << " confidence: " << confidence[victim_way][cpu] << " pf_offset: " << offset[victim_way][cpu] << " pf_delta: " << delta[victim_way][cpu] << endl;
    );

    	  valid[victim_way][cpu] = 1;
    		sig[victim_way][cpu] = pf_sig;
     confidence[victim_way][cpu] = pf_confidence;
    	 offset[victim_way][cpu] = pf_offset;
    	  delta[victim_way][cpu] = pf_delta;
}

uint32_t GLOBAL_REGISTER_L::check_entry(uint32_t page_offset, int cpu)
{
    uint32_t max_conf = 0,
             max_conf_way = MAX_GHR_ENTRY;

    for (uint32_t i = 0; i < MAX_GHR_ENTRY; i++) {
        if ((offset[i][cpu] == page_offset) && (max_conf < confidence[i][cpu])) {
            max_conf = confidence[i][cpu];
            max_conf_way = i;
        }
    }

    return max_conf_way;
}
