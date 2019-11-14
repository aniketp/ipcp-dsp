/***************************************************************************
For the Third Data Prefetching Championship - DPC3

Paper ID: #4
Instruction Pointer Classifying Prefetcher - IPCP

Authors: 
Samuel Pakalapati - samuelpakalapati@gmail.com
Biswabandan Panda - biswap@cse.iitk.ac.in
***************************************************************************/

#include "cache.h"

//--------------------- Sangam --------------------------------
#define LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 0
#define NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 (1 << LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1)
#define NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1 40

#define IP_TABLE_TAG_MASK 0x3fffULL                 // For prefetch degree 4
#define LOG_NUM_SETS_IN_IP_TABLE_L1 7
#define NUM_SETS_IN_IP_TABLE_L1 (1 << LOG_NUM_SETS_IN_IP_TABLE_L1)
#define NUM_WAYS_IN_IP_TABLE_L1 15

#define IP_DELTA_TABLE_TAG_MASK 0xffffULL        // For prefetch degree 4
#define LOG_NUM_SETS_IN_IP_DELTA_TABLE_L1 8
#define NUM_SETS_IN_IP_DELTA_TABLE_L1 (1 << LOG_NUM_SETS_IN_IP_DELTA_TABLE_L1)
#define NUM_WAYS_IN_IP_DELTA_TABLE_L1 8

#define SATURATING_COUNTER_MAX_L1 3
#define PREDICTION_THRESHOLD_L1 ((NUM_CPUS == 1) ? 2 : 3)
#define BASE_PREFETCH_DEGREE_L1 4

#define PAGE_OFFSET_MASK ((1 << LOG2_PAGE_SIZE) - 1)

#define NUM_ENTRIES_IN_NL_BUFFER_L1 64
#define NL_THRESHOLD_NUMER_L1 1
#define NL_THRESHOLD_DENOM_L1 4

#define THROTTLE_LEVEL_MAX_L1 (2 + BASE_PREFETCH_DEGREE_L1)

#define POINTER_LAST true
#define POINTER_NON_LAST false
#define STRIDE_CONF_MAX 3
#define STRIDE_CONF_THRESHOLD 3

#define PARTIAL_IP_MASK 0x7fULL

#define NUM_STRIDES_IN_LONG_HIST_IP_TABLE 20
#define LONG_HIST_IP_TABLE_TAG_MASK 0x1fffffULL
#define NUM_ENTRIES_IN_LONG_HIST_IP_TABLE 32
#define LONG_HIST_MATCH_LENGTH 1

unsigned char throttle_level_L1;

// Next line buffer entry: 73 bits
class NLBufferL1 {
   public:
      // Tag
      uint64_t tag;          // 64 bits

      // LRU states
      uint64_t lru;          // 6 bits

      // Valid bit 
      bool valid;            // 1 bit

      // Degree
      unsigned char degree;  // 2 bits

      NLBufferL1 () {}
};

NLBufferL1 nlBufferL1[NUM_ENTRIES_IN_NL_BUFFER_L1];       // 64x73 bits = 4672 bits
unsigned degreeInsertionsL1[BASE_PREFETCH_DEGREE_L1];  // 32x4 bits = 128 bits
unsigned degreeHitsL1[BASE_PREFETCH_DEGREE_L1];        // 32x4 bits = 128 bits

// Recent access tag array entry: 71 bits
class RecentAccessTagArrayL1 {
   public:
      // Tag
      uint64_t tag;    // 64 bits

      // LRU states
      uint64_t lru;    // 6 bits

      // Valid bit
      bool valid;      // 1 bit

      RecentAccessTagArrayL1 () {}
};

RecentAccessTagArrayL1 recentAccessTagArrayL1[NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1][NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1];  // 1x40x71 bits = 2840 bits

// IP stride table entry: 63 bits
class IPtableL1 {
   public:
      // Tag
      uint64_t tag;					// 63 - (3+4+1+6+7*(BASE_PREFETCH_DEGREE_L1+1)) bits

      // LRU states
      uint64_t lru;					// 4 bits

      // Valid bit
      bool valid;					// 1 bit

      // Last seen offset in a page
      unsigned char offset;				// 6 bits

      // Last seen strides within a page
      char stride[BASE_PREFETCH_DEGREE_L1+1];		// 7*(BASE_PREFETCH_DEGREE_L1+1)

      // Stride confidence
      unsigned char conf;				// 2 bits

      // Stride confidence pointer
      bool confPointer;					// 1 bit

      IPtableL1 () {}
};

IPtableL1 ipTableL1[NUM_SETS_IN_IP_TABLE_L1][NUM_WAYS_IN_IP_TABLE_L1];  // 128x15x63 bits = 120960 bits
char ipTableStride[BASE_PREFETCH_DEGREE_L1+1];                               // 7*5 bits = 35 bits

// IP-delta correlation table entry: 64 bits
class IPDeltaTableL1 {
   public:
      // Tag
      uint64_t tag;					// 64 - (8+3+1+9*BASE_PREFETCH_DEGREE_L1) bits

      // LRU states
      uint64_t lru;					// 3 bits

      // Valid bit
      bool valid;					// 1 bit

      // Last seen strides within a page
      char stride[BASE_PREFETCH_DEGREE_L1];			// 7*BASE_PREFETCH_DEGREE_L1 bits

      // Confidence counters
      unsigned char counters[BASE_PREFETCH_DEGREE_L1];	// 2*BASE_PREFETCH_DEGREE_L1 bits

      // Partial IP for last prediction
      uint64_t partial_ip;				// 7 bits
      bool partial_ip_valid;				// 1 bit

      IPDeltaTableL1 () {}
};

char ipPrefetchStride[BASE_PREFETCH_DEGREE_L1+1];   // 7*5 bits = 35 bits
IPDeltaTableL1 ipDeltaTableL1[NUM_SETS_IN_IP_DELTA_TABLE_L1][NUM_WAYS_IN_IP_DELTA_TABLE_L1];  // 256x8x64 bits = 131072 bits

// Long history IP stride table entry: 225 bits
class LongHistIPtableL1 {
   public:
      // Tag
      uint64_t ip;                                      // 225 - (5+1+58+7*NUM_STRIDES_IN_LONG_HIST_TABLE) bits

      // LRU states
      uint64_t lru;                                     // 5 bits

      // Valid bit
      bool valid;                                       // 1 bit

      // Last seen block address
      uint64_t block_addr;				// 58 bits

      // Last seen many strides within a page
      char stride[NUM_STRIDES_IN_LONG_HIST_IP_TABLE];      // 7*NUM_STRIDES_IN_LONG_HIST_IP_TABLE

      LongHistIPtableL1 () {}
};

LongHistIPtableL1 longHistIPTableL1[NUM_ENTRIES_IN_LONG_HIST_IP_TABLE];		// 225*32 bits = 7200 bits
char longHistory[NUM_ENTRIES_IN_LONG_HIST_IP_TABLE+1];				// 7*(NUM_ENTRIES_IN_LONG_HIST_IP_TABLE+1) bits

// ----------------------- Bouquet ---------------------------------

#define NUM_IP_TABLE_L1_ENTRIES 1024                        // IP table entries 
#define NUM_GHB_ENTRIES 16                                  // Entries in the GHB
#define NUM_IP_INDEX_BITS 10                                // Bits to index into the IP table 
#define NUM_IP_TAG_BITS 6                                   // Tag bits per IP table entry
#define S_TYPE 1                                            // stream
#define CS_TYPE 2                                           // constant stride
#define CPLX_TYPE 3                                         // complex stride
#define NL_TYPE 4                                           // next line

// #define SIG_DEBUG_PRINT
#ifdef SIG_DEBUG_PRINT
#define SIG_DP(x) x
#else
#define SIG_DP(x)
#endif

class IP_TABLE_L1 {
  public:
    uint64_t ip_tag;
    uint64_t last_page;                                     // last page seen by IP
    uint64_t last_cl_offset;                                // last cl offset in the 4KB page
    int64_t last_stride;                                    // last delta observed
    uint16_t ip_valid;                                      // Valid IP or not   
    int conf;                                               // CS conf
    uint16_t signature;                                     // CPLX signature
    uint16_t str_dir;                                       // stream direction
    uint16_t str_valid;                                     // stream valid
    uint16_t str_strength;                                  // stream strength

    IP_TABLE_L1 () {
        ip_tag = 0;
        last_page = 0;
        last_cl_offset = 0;
        last_stride = 0;
        ip_valid = 0;
        signature = 0;
        conf = 0;
        str_dir = 0;
        str_valid = 0;
        str_strength = 0;
    };
};

class DELTA_PRED_TABLE {
public:
    int delta;
    int conf;

    DELTA_PRED_TABLE () {
        delta = 0;
        conf = 0;
    };        
};


IP_TABLE_L1 trackers_l1[NUM_CPUS][NUM_IP_TABLE_L1_ENTRIES];
DELTA_PRED_TABLE DPT_l1[NUM_CPUS][4096];
uint64_t ghb_l1[NUM_CPUS][NUM_GHB_ENTRIES];
uint64_t prev_cpu_cycle[NUM_CPUS];
uint64_t num_misses[NUM_CPUS];
float mpkc[NUM_CPUS] = {0};
int spec_nl[NUM_CPUS] = {0};


////#####################################################################////
//------------------- Helper Methods ----------------------------

// ------------------------ Sangam -------------------------------

static void NLBufferL1Insert (uint64_t cl_addr, int current_degree_index)
{
   int j;
   for (j=0; j<NUM_ENTRIES_IN_NL_BUFFER_L1; j++) {
      if (nlBufferL1[j].valid && (nlBufferL1[j].tag == cl_addr)) {
         if (nlBufferL1[j].degree > (current_degree_index+1)) {
            assert(degreeInsertionsL1[nlBufferL1[j].degree-1] > 0);
            degreeInsertionsL1[nlBufferL1[j].degree-1]--;
            nlBufferL1[j].degree = current_degree_index+1;   // Always favor smaller degree NL prefetcher
            degreeInsertionsL1[current_degree_index]++;
            for (int jj=0; jj<NUM_ENTRIES_IN_NL_BUFFER_L1; jj++) {
               nlBufferL1[jj].lru++;
            }
            nlBufferL1[j].lru = 0;
         }
         break;
      }
   }
   if (j == NUM_ENTRIES_IN_NL_BUFFER_L1) {  // MISS
      for (j=0; j<NUM_ENTRIES_IN_NL_BUFFER_L1; j++) {
         if (!nlBufferL1[j].valid) break;
      }
      if (j == NUM_ENTRIES_IN_NL_BUFFER_L1) {
         uint64_t maxlru = 0;
         int repl_index = -1;
         for (j=0; j<NUM_ENTRIES_IN_NL_BUFFER_L1; j++) {
            if (nlBufferL1[j].lru >= maxlru) {
               maxlru = nlBufferL1[j].lru;
               repl_index = j;
            }
         }
         j = repl_index;
      }
      nlBufferL1[j].tag = cl_addr;
      nlBufferL1[j].degree = current_degree_index+1;
      nlBufferL1[j].valid = true;
      degreeInsertionsL1[current_degree_index]++;
      for (int jj=0; jj<NUM_ENTRIES_IN_NL_BUFFER_L1; jj++) {
         nlBufferL1[jj].lru++;
      }
      nlBufferL1[j].lru = 0;
   }
}

static void RecentAccessTagArrayL1LookupAndInsertIfMiss (uint64_t cl_addr)
{
   int recentAccessTagArrayIndex = cl_addr & (NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 - 1);
   uint64_t recentAccessTagArrayTag = cl_addr >> LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1;
   int ii;
   for (ii=0; ii<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; ii++) {
      if (recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].valid && (recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].tag == recentAccessTagArrayTag)) {
         break;
      }
   }
   if (ii == NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1) {
      for (ii=0; ii<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; ii++) {
         if (!recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].valid) break;
      }
      if (ii == NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1) {
         uint64_t maxlru = 0;
         int repl_index = -1;
         for (ii=0; ii<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; ii++) {
            if (recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].lru > maxlru) {
               maxlru = recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].lru;
               repl_index = ii;
            }
         }
         ii = repl_index;
      }
      recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].tag = recentAccessTagArrayTag;
      recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].valid = true;
   }
   for (int jj=0; jj<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; jj++) recentAccessTagArrayL1[recentAccessTagArrayIndex][jj].lru++;
   recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].lru = 0;
}

static bool RecentAccessTagArrayL1DetermineHit (uint64_t cl_addr)
{
   int recentAccessTagArrayIndex = cl_addr & (NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 - 1);
   uint64_t recentAccessTagArrayTag = cl_addr >> LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1;
   int ii;
   for (ii=0; ii<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; ii++) {
      if (recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].valid && (recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].tag == recentAccessTagArrayTag)) {
         return true;
      }
   }
   return false;
}

static void RecentAccessTagArrayL1Insert (uint64_t cl_addr)
{
   int recentAccessTagArrayIndex = cl_addr & (NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 - 1);
   uint64_t recentAccessTagArrayTag = cl_addr >> LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1;

   int ii;
   
   for (ii=0; ii<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; ii++) {
      if (!recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].valid) break;
   }
   if (ii == NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1) {
      uint64_t maxlru = 0;
      int repl_index = -1;
      for (ii=0; ii<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; ii++) {
         if (recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].lru > maxlru) {
            maxlru = recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].lru;
            repl_index = ii;
         }
      }
      ii = repl_index;
   }
   recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].tag = recentAccessTagArrayTag;
   recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].valid = true;
   for (int jj=0; jj<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; jj++) recentAccessTagArrayL1[recentAccessTagArrayIndex][jj].lru++;
   recentAccessTagArrayL1[recentAccessTagArrayIndex][ii].lru = 0;
}


// ------------------------------ Bouquet ---------------------------------


/***************Updating the signature*************************************/ 
uint16_t update_sig_l1(uint16_t old_sig, int delta){                           
    uint16_t new_sig = 0;
    int sig_delta = 0;

// 7-bit sign magnitude form, since we need to track deltas from +63 to -63
    sig_delta = (delta < 0) ? (((-1) * delta) + (1 << 6)) : delta;
    new_sig = ((old_sig << 1) ^ sig_delta) & 0xFFF;                     // 12-bit signature

    return new_sig;
}



/****************Encoding the metadata***********************************/
uint32_t encode_metadata(int stride, uint16_t type, int spec_nl){

uint32_t metadata = 0;

// first encode stride in the last 8 bits of the metadata
if(stride > 0)
    metadata = stride;
else
    metadata = ((-1*stride) | 0b1000000);

// encode the type of IP in the next 4 bits
metadata = metadata | (type << 8);

// encode the speculative NL bit in the next 1 bit
metadata = metadata | (spec_nl << 12);

return metadata;

}


/*********************Checking for a global stream (GS class)***************/

void check_for_stream_l1(int index, uint64_t cl_addr, uint8_t cpu){
int pos_count=0, neg_count=0, count=0;
uint64_t check_addr = cl_addr;

// check for +ve stream
    for(int i=0; i<NUM_GHB_ENTRIES; i++){
        check_addr--;
        for(int j=0; j<NUM_GHB_ENTRIES; j++)
            if(check_addr == ghb_l1[cpu][j]){
                pos_count++;
                break;
            }
    }

check_addr = cl_addr;
// check for -ve stream
    for(int i=0; i<NUM_GHB_ENTRIES; i++){
        check_addr++;
        for(int j=0; j<NUM_GHB_ENTRIES; j++)
            if(check_addr == ghb_l1[cpu][j]){
                neg_count++;
                break;
            }
    }

    if(pos_count > neg_count){                                // stream direction is +ve
        trackers_l1[cpu][index].str_dir = 1;
        count = pos_count;
    }
    else{                                                     // stream direction is -ve
        trackers_l1[cpu][index].str_dir = 0;
        count = neg_count;
    }

if(count > NUM_GHB_ENTRIES/2){                                // stream is detected
    trackers_l1[cpu][index].str_valid = 1;
    if(count >= (NUM_GHB_ENTRIES*3)/4)                        // stream is classified as strong if more than 3/4th entries belong to stream
        trackers_l1[cpu][index].str_strength = 1;
}
else{
    if(trackers_l1[cpu][index].str_strength == 0)             // if identified as weak stream, we need to reset
        trackers_l1[cpu][index].str_valid = 0;
}

}

/**************************Updating confidence for the CS class****************/
int update_conf(int stride, int pred_stride, int conf){
    if(stride == pred_stride){             // use 2-bit saturating counter for confidence
        conf++;
        if(conf > 3)
            conf = 3;
    } else {
        conf--;
        if(conf < 0)
            conf = 0;
    }

return conf;
}




///////////////////////////////////////////////////////////////////////////
//////////////////////       MAIN ROUTINE        //////////////////////////
///////////////////////////////////////////////////////////////////////////



void CACHE::l1d_prefetcher_initialize() 
{
    cout << "CPU " << cpu << " MIX OF BOUQUET and SANGAM" << endl;

   for (int i=0; i<NUM_ENTRIES_IN_NL_BUFFER_L1; i++) {
      nlBufferL1[i].valid = false;
      nlBufferL1[i].lru = 0;
   }

   for (int i=0; i<BASE_PREFETCH_DEGREE_L1; i++) {
      degreeInsertionsL1[i] = 0;
      degreeHitsL1[i] = 0;
   }

   for (int i=0; i<NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1; i++) {
      for (int j=0; j<NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1; j++) {
         recentAccessTagArrayL1[i][j].valid = false;
         recentAccessTagArrayL1[i][j].lru = 0; 
      }
   }

   for (int i=0; i<NUM_SETS_IN_IP_TABLE_L1; i++) {
      for (int j=0; j<NUM_WAYS_IN_IP_TABLE_L1; j++) {
         ipTableL1[i][j].valid = false;
         ipTableL1[i][j].lru = 0;
      }
   }

   for (int i=0; i<NUM_SETS_IN_IP_DELTA_TABLE_L1; i++) {
      for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L1; j++) {
         ipDeltaTableL1[i][j].valid = false;
         ipDeltaTableL1[i][j].lru = 0;
      }
   }

   for (int i=0; i<NUM_ENTRIES_IN_LONG_HIST_IP_TABLE; i++) {
      longHistIPTableL1[i].valid = false;
      longHistIPTableL1[i].lru = 0;
   }

   throttle_level_L1 = 0;
}

void CACHE::l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type)
{
    int i = 0;
    if (type == PREFETCH) return;	// Don't update or inject new prefetches

    uint64_t curr_page = addr >> LOG2_PAGE_SIZE;
    uint64_t cl_addr = addr >> LOG2_BLOCK_SIZE;
    uint64_t cl_offset = (addr >> LOG2_BLOCK_SIZE) & 0x3F;
    uint16_t signature = 0, last_signature = 0;
    int prefetch_degree = 0;
    int spec_nl_threshold = 0;
    int num_prefs = 0;
    uint32_t metadata=0;
    uint16_t ip_tag = (ip >> NUM_IP_INDEX_BITS) & ((1 << NUM_IP_TAG_BITS)-1);

    uint64_t pageid = addr >> LOG2_PAGE_SIZE;
    unsigned char offset = (addr & PAGE_OFFSET_MASK) >> LOG2_BLOCK_SIZE;
    bool did_pref = false;
    bool current_delta_nonzero = false;

    if(NUM_CPUS == 1) {
        prefetch_degree = 3;
        spec_nl_threshold = 15; 
    } else {                                    // tightening the degree and MPKC constraints for multi-core
        prefetch_degree = 2;
        spec_nl_threshold = 5;
    }

    // update miss counter
    if(cache_hit == 0)
        num_misses[cpu] += 1;

    // update spec nl bit when num misses crosses certain threshold
    if(num_misses[cpu] == 256){
        mpkc[cpu] = ((float) num_misses[cpu]/(current_core_cycle[cpu]-prev_cpu_cycle[cpu]))*1000;
        prev_cpu_cycle[cpu] = current_core_cycle[cpu];
        if(mpkc[cpu] > spec_nl_threshold)
            spec_nl[cpu] = 0;
        else
            spec_nl[cpu] = 1;
        num_misses[cpu] = 0;
    }

    //
    //  Long history IP table lookup
    //
    char longHistIPTableNewDelta = 0;
    for (i=0; i<NUM_ENTRIES_IN_LONG_HIST_IP_TABLE; i++) {
        if (longHistIPTableL1[i].valid && (longHistIPTableL1[i].ip == (ip & LONG_HIST_IP_TABLE_TAG_MASK))) {
            for (int j=0; j<NUM_STRIDES_IN_LONG_HIST_IP_TABLE; j++) longHistory[j] = longHistIPTableL1[i].stride[j];
            if (pageid == (longHistIPTableL1[i].block_addr >> (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) {
                longHistory[NUM_STRIDES_IN_LONG_HIST_IP_TABLE] = offset - (((longHistIPTableL1[i].block_addr << LOG2_BLOCK_SIZE) & PAGE_OFFSET_MASK) >> LOG2_BLOCK_SIZE);
            }
            else longHistory[NUM_STRIDES_IN_LONG_HIST_IP_TABLE] = 0;
            longHistIPTableNewDelta = longHistory[NUM_STRIDES_IN_LONG_HIST_IP_TABLE];
            longHistIPTableL1[i].block_addr = cl_addr;
            if (longHistIPTableNewDelta != 0) {
                for (int j=0; j<NUM_STRIDES_IN_LONG_HIST_IP_TABLE-1; j++) longHistIPTableL1[i].stride[j] = longHistIPTableL1[i].stride[j+1];
                longHistIPTableL1[i].stride[NUM_STRIDES_IN_LONG_HIST_IP_TABLE-1] = longHistIPTableNewDelta;
            }
            break;
        }
    }
    if (i == NUM_ENTRIES_IN_LONG_HIST_IP_TABLE) {
        for (i=0; i<NUM_ENTRIES_IN_LONG_HIST_IP_TABLE; i++) {
            if (!longHistIPTableL1[i].valid) break;
        }
        if (i == NUM_ENTRIES_IN_LONG_HIST_IP_TABLE) {
            uint64_t maxlru = 0;
            int rep_index;
            for (i=0; i<NUM_ENTRIES_IN_LONG_HIST_IP_TABLE; i++) {
                if (longHistIPTableL1[i].lru >= maxlru) {
                maxlru = longHistIPTableL1[i].lru;
                rep_index = i;
            }
            }
            i = rep_index;
        }
        assert(i < NUM_ENTRIES_IN_LONG_HIST_IP_TABLE);
        longHistIPTableL1[i].ip = (ip & LONG_HIST_IP_TABLE_TAG_MASK);
        longHistIPTableL1[i].block_addr = cl_addr;
        for (int j=0; j<NUM_STRIDES_IN_LONG_HIST_IP_TABLE; j++) longHistIPTableL1[i].stride[j] = 0;
        longHistIPTableL1[i].valid = true;
    }
    for (int j=0; j<NUM_ENTRIES_IN_LONG_HIST_IP_TABLE; j++)
        longHistIPTableL1[j].lru++;
    longHistIPTableL1[i].lru = 0;

    //
    // NL Buffer Lookup
    //
    for (i=0; i<NUM_ENTRIES_IN_NL_BUFFER_L1; i++) {
        // Determine next line prefetcher's usefulness
        if (nlBufferL1[i].valid && (nlBufferL1[i].tag == cl_addr)) {
            degreeHitsL1[nlBufferL1[i].degree-1]++;
            nlBufferL1[i].valid = false;
            break;
        }
    }

    //
    // IP Table Lookup
    //
    bool constantStrideValid = false;
    char constantStride = 0;
    int ipTableIndex = (pageid) & (NUM_SETS_IN_IP_TABLE_L1 - 1);
    uint64_t ipTableTag = ((pageid) >> LOG_NUM_SETS_IN_IP_TABLE_L1) & IP_TABLE_TAG_MASK;
    int ii;
    for (ii=0; ii<NUM_WAYS_IN_IP_TABLE_L1; ii++) {
        if (ipTableL1[ipTableIndex][ii].valid && (ipTableL1[ipTableIndex][ii].tag == ipTableTag)) {
            if ((signed)(offset - ipTableL1[ipTableIndex][ii].offset) != 0) {
                current_delta_nonzero = true;
                for (i=0; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
                if (ipTableL1[ipTableIndex][ii].stride[i] == 0) {
                    ipTableL1[ipTableIndex][ii].stride[i] = (signed)(offset - ipTableL1[ipTableIndex][ii].offset);
                    break;
                }
                }
                if (i == BASE_PREFETCH_DEGREE_L1+1) {
                for (i=0; i<BASE_PREFETCH_DEGREE_L1; i++) {
                    ipTableL1[ipTableIndex][ii].stride[i] = ipTableL1[ipTableIndex][ii].stride[i+1];
                }
                ipTableL1[ipTableIndex][ii].stride[i] = (signed)(offset - ipTableL1[ipTableIndex][ii].offset);
                }
                
                if (i == 0) {
                ipTableL1[ipTableIndex][ii].conf = 0;
                ipTableL1[ipTableIndex][ii].confPointer = POINTER_LAST;
                }
                else if (ipTableL1[ipTableIndex][ii].stride[i] == ipTableL1[ipTableIndex][ii].stride[i-1]) {
                if (ipTableL1[ipTableIndex][ii].confPointer == POINTER_LAST) {
                    if (ipTableL1[ipTableIndex][ii].conf < STRIDE_CONF_MAX) ipTableL1[ipTableIndex][ii].conf++;
                }
                else {
                    ipTableL1[ipTableIndex][ii].conf = 1;
                    ipTableL1[ipTableIndex][ii].confPointer = POINTER_LAST;
                }
                }
                else {
                if (ipTableL1[ipTableIndex][ii].confPointer == POINTER_LAST) {
                    ipTableL1[ipTableIndex][ii].confPointer = POINTER_NON_LAST;
                if (ipTableL1[ipTableIndex][ii].conf > 0) ipTableL1[ipTableIndex][ii].conf--;
                }
                else {
                    assert(i > 1);
                    ipTableL1[ipTableIndex][ii].confPointer = POINTER_LAST;
                    if (ipTableL1[ipTableIndex][ii].stride[i] == ipTableL1[ipTableIndex][ii].stride[i-2]) {
                        if (ipTableL1[ipTableIndex][ii].conf < STRIDE_CONF_MAX) ipTableL1[ipTableIndex][ii].conf++;
                    }
                    else {
                        ipTableL1[ipTableIndex][ii].conf = 0;
                    }  
                }
                }

            if ((ipTableL1[ipTableIndex][ii].conf >= STRIDE_CONF_THRESHOLD) && (ipTableL1[ipTableIndex][ii].stride[i] == ipTableL1[ipTableIndex][ii].stride[i-1])) {
                constantStride = ipTableL1[ipTableIndex][ii].stride[i];
                constantStrideValid = true;
                }
                
                ipTableL1[ipTableIndex][ii].offset = offset;
            }
            break;
        }
    }
    if (ii == NUM_WAYS_IN_IP_TABLE_L1) {
        for (ii=0; ii<NUM_WAYS_IN_IP_TABLE_L1; ii++) {
            if (!ipTableL1[ipTableIndex][ii].valid) break;
        }
        if (ii == NUM_WAYS_IN_IP_TABLE_L1) {
            uint64_t maxlru = 0;
            int repl_index = -1;
            for (ii=0; ii<NUM_WAYS_IN_IP_TABLE_L1; ii++) {
                if (ipTableL1[ipTableIndex][ii].lru > maxlru) {
                    maxlru = ipTableL1[ipTableIndex][ii].lru;
                    repl_index = ii;
                }
            }
            ii = repl_index;
        }
        ipTableL1[ipTableIndex][ii].tag = ipTableTag;
        ipTableL1[ipTableIndex][ii].offset = offset;
        for (i=0; i<BASE_PREFETCH_DEGREE_L1+1; i++)
            ipTableL1[ipTableIndex][ii].stride[i] = 0;
        ipTableL1[ipTableIndex][ii].valid = true;
    }
    for (i=0; i<NUM_WAYS_IN_IP_TABLE_L1; i++) ipTableL1[ipTableIndex][i].lru++;
    ipTableL1[ipTableIndex][ii].lru = 0;
    int lastNonZeroIndex = -1;
    for (i=0; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
        ipTableStride[i] = ipTableL1[ipTableIndex][ii].stride[i];
        if (ipTableStride[i] != 0) lastNonZeroIndex = i;
    }

    //
    // IP delta table lookup and training
    //
    if (current_delta_nonzero) {
      for (i=0; i<lastNonZeroIndex; i++) {
         assert(ipTableStride[i] != 0);
         unsigned delta;
         delta = (ipTableStride[i] >= 0) ? ipTableStride[i] : ((-ipTableStride[i]) | (1 <<  (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)));
         int ipDeltaTableIndex = (((pageid) << 3) ^ delta) & (NUM_SETS_IN_IP_DELTA_TABLE_L1 - 1);
         uint64_t ipDeltaTableTag = ((((pageid) << 3) ^ delta) >> LOG_NUM_SETS_IN_IP_DELTA_TABLE_L1) & IP_DELTA_TABLE_TAG_MASK;
         for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L1; ii++) {
            if (ipDeltaTableL1[ipDeltaTableIndex][ii].valid && (ipDeltaTableL1[ipDeltaTableIndex][ii].tag == ipDeltaTableTag)) {
               if (ipTableStride[lastNonZeroIndex] == ipDeltaTableL1[ipDeltaTableIndex][ii].stride[lastNonZeroIndex-i-1]) {
                  if (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[lastNonZeroIndex-i-1] < SATURATING_COUNTER_MAX_L1) {
                     ipDeltaTableL1[ipDeltaTableIndex][ii].counters[lastNonZeroIndex-i-1]++;
                  }
               }
               else {
                  ipDeltaTableL1[ipDeltaTableIndex][ii].stride[lastNonZeroIndex-i-1] = ipTableStride[lastNonZeroIndex];
                  ipDeltaTableL1[ipDeltaTableIndex][ii].counters[lastNonZeroIndex-i-1] = 1;
               }
               break;
            }
         }
         if (ii == NUM_WAYS_IN_IP_DELTA_TABLE_L1) {
            for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L1; ii++) {
               if (!ipDeltaTableL1[ipDeltaTableIndex][ii].valid) break;
            }
            if (ii == NUM_WAYS_IN_IP_DELTA_TABLE_L1) {
               uint64_t maxlru = 0;
               int repl_index = -1;
               for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L1; ii++) {
                  if (ipDeltaTableL1[ipDeltaTableIndex][ii].lru > maxlru) {
                     maxlru = ipDeltaTableL1[ipDeltaTableIndex][ii].lru;
                     repl_index = ii;
                  }
               }
               ii = repl_index;
            }
            ipDeltaTableL1[ipDeltaTableIndex][ii].partial_ip_valid = false;
            ipDeltaTableL1[ipDeltaTableIndex][ii].tag = ipDeltaTableTag;
            for (int j=0; j<BASE_PREFETCH_DEGREE_L1; j++) {
               if (i+j+1 < BASE_PREFETCH_DEGREE_L1+1) {
                  ipDeltaTableL1[ipDeltaTableIndex][ii].stride[j] = ipTableStride[i+j+1];
                  if (ipTableStride[i+j+1] != 0) ipDeltaTableL1[ipDeltaTableIndex][ii].counters[j] = 1;
                  else ipDeltaTableL1[ipDeltaTableIndex][ii].counters[j] = 0;
               }
               else {
                  ipDeltaTableL1[ipDeltaTableIndex][ii].stride[j] = 0;
                  ipDeltaTableL1[ipDeltaTableIndex][ii].counters[j] = 0;
               }
            }
            ipDeltaTableL1[ipDeltaTableIndex][ii].valid = true;
         }
         for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L1; j++) ipDeltaTableL1[ipDeltaTableIndex][j].lru++;
         ipDeltaTableL1[ipDeltaTableIndex][ii].lru = 0;
      }
   }
    RecentAccessTagArrayL1LookupAndInsertIfMiss (cl_addr);

    
    //-----------------------------------------------------------------------
    //-------------------- Main Prefetching decision ------------------------
    //-----------------------------------------------------------------------


    // calculate the index bit
    int index = ip & ((1 << NUM_IP_INDEX_BITS)-1);
    if(trackers_l1[cpu][index].ip_tag != ip_tag){               // new/conflict IP
        if(trackers_l1[cpu][index].ip_valid == 0){              // if valid bit is zero, update with latest IP info
        trackers_l1[cpu][index].ip_tag = ip_tag;
        trackers_l1[cpu][index].last_page = curr_page;
        trackers_l1[cpu][index].last_cl_offset = cl_offset;
        trackers_l1[cpu][index].last_stride = 0;
        trackers_l1[cpu][index].signature = 0;
        trackers_l1[cpu][index].conf = 0;
        trackers_l1[cpu][index].str_valid = 0;
        trackers_l1[cpu][index].str_strength = 0;
        trackers_l1[cpu][index].str_dir = 0;
        trackers_l1[cpu][index].ip_valid = 1;
    } else {                                                    // otherwise, reset valid bit and leave the previous IP as it is
        trackers_l1[cpu][index].ip_valid = 0;
    }

    // issue a next line prefetch upon encountering new IP
        uint64_t pf_address = ((addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE; // BASE NL=1, changing it to 3
        metadata = encode_metadata(1, NL_TYPE, spec_nl[cpu]);
        prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
        return;
    }
    else {                                                     // if same IP encountered, set valid bit
        trackers_l1[cpu][index].ip_valid = 1;
    }
    

    // calculate the stride between the current address and the last address
    int64_t stride = 0;
    if (cl_offset > trackers_l1[cpu][index].last_cl_offset)
        stride = cl_offset - trackers_l1[cpu][index].last_cl_offset;
    else {
        stride = trackers_l1[cpu][index].last_cl_offset - cl_offset;
        stride *= -1;
    }

    // don't do anything if same address is seen twice in a row
    if (stride == 0)
        return;


    // page boundary learning
    if(curr_page != trackers_l1[cpu][index].last_page){
        if(stride < 0)
            stride += 64;
        else
            stride -= 64;
    }

    // update constant stride(CS) confidence
    trackers_l1[cpu][index].conf = update_conf(stride, trackers_l1[cpu][index].last_stride, trackers_l1[cpu][index].conf);

    // update CS only if confidence is zero
    if(trackers_l1[cpu][index].conf == 0)                      
        trackers_l1[cpu][index].last_stride = stride;

    last_signature = trackers_l1[cpu][index].signature;
    // update complex stride(CPLX) confidence
    DPT_l1[cpu][last_signature].conf = update_conf(stride, DPT_l1[cpu][last_signature].delta, DPT_l1[cpu][last_signature].conf);

    // update CPLX only if confidence is zero
    if(DPT_l1[cpu][last_signature].conf == 0)
        DPT_l1[cpu][last_signature].delta = stride;

    // calculate and update new signature in IP table
    signature = update_sig_l1(last_signature, stride);
    trackers_l1[cpu][index].signature = signature;

    // check GHB for stream IP
    check_for_stream_l1(index, cl_addr, cpu);           

    SIG_DP(
    cout << ip << ", " << cache_hit << ", " << cl_addr << ", " << addr << ", " << stride << "; ";
    cout << last_signature<< ", "  << DPT_l1[cpu][last_signature].delta<< ", "  << DPT_l1[cpu][last_signature].conf << "; ";
    cout << trackers_l1[cpu][index].last_stride << ", " << stride << ", " << trackers_l1[cpu][index].conf << ", " << "; ";
    );

    if(trackers_l1[cpu][index].str_valid == 1){                         // stream IP
        // for stream, prefetch with twice the usual degree
            prefetch_degree = prefetch_degree*2;
        for (int i=0; i<prefetch_degree; i++) {
            uint64_t pf_address = 0;

            if(trackers_l1[cpu][index].str_dir == 1){                   // +ve stream
                pf_address = (cl_addr + i + 1) << LOG2_BLOCK_SIZE;
                metadata = encode_metadata(1, S_TYPE, spec_nl[cpu]);    // stride is 1
            }
            else{                                                       // -ve stream
                pf_address = (cl_addr - i - 1) << LOG2_BLOCK_SIZE;
                metadata = encode_metadata(-1, S_TYPE, spec_nl[cpu]);   // stride is -1
            }

            // Check if prefetch address is in same 4 KB page
            if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)){
                break;
            }

            prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
            num_prefs++;
            SIG_DP(cout << "1, ");
            }

    } else if(trackers_l1[cpu][index].conf > 1 && trackers_l1[cpu][index].last_stride != 0){            // CS IP  
        for (int i=0; i<prefetch_degree; i++) {
            uint64_t pf_address = (cl_addr + (trackers_l1[cpu][index].last_stride*(i+1))) << LOG2_BLOCK_SIZE;

            // Check if prefetch address is in same 4 KB page
            if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)){
                break;
            }

            metadata = encode_metadata(trackers_l1[cpu][index].last_stride, CS_TYPE, spec_nl[cpu]);
            prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
            num_prefs++;
            SIG_DP(cout << trackers_l1[cpu][index].last_stride << ", ");
        }
    } else if(DPT_l1[cpu][signature].conf >= 0 && DPT_l1[cpu][signature].delta != 0) {  // if conf>=0, continue looking for delta
        int pref_offset = 0,i=0;                                                        // CPLX IP
        for (i=0; i<prefetch_degree; i++) {
            pref_offset += DPT_l1[cpu][signature].delta;
            uint64_t pf_address = ((cl_addr + pref_offset) << LOG2_BLOCK_SIZE);

            // Check if prefetch address is in same 4 KB page
            if (((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) || 
                    (DPT_l1[cpu][signature].conf == -1) ||
                    (DPT_l1[cpu][signature].delta == 0)){
                // if new entry in DPT or delta is zero, break
                break;
            }

            // we are not prefetching at L2 for CPLX type, so encode delta as 0
            metadata = encode_metadata(0, CPLX_TYPE, spec_nl[cpu]);
            if(DPT_l1[cpu][signature].conf > 0){                                 // prefetch only when conf>0 for CPLX
                prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
                num_prefs++;
                SIG_DP(cout << pref_offset << ", ");
            }
            signature = update_sig_l1(signature, DPT_l1[cpu][signature].delta);
        }
    } 

    // update the IP table entries
    trackers_l1[cpu][index].last_cl_offset = cl_offset;
    trackers_l1[cpu][index].last_page = curr_page;

    // update GHB
    // search for matching cl addr
    int ghb_index=0;
    for(ghb_index = 0; ghb_index < NUM_GHB_ENTRIES; ghb_index++)
        if(cl_addr == ghb_l1[cpu][ghb_index])
            break;
    // only update the GHB upon finding a new cl address
    if(ghb_index == NUM_GHB_ENTRIES){
        for(ghb_index=NUM_GHB_ENTRIES-1; ghb_index>0; ghb_index--)
            ghb_l1[cpu][ghb_index] = ghb_l1[cpu][ghb_index-1];
        ghb_l1[cpu][0] = cl_addr;
    }

    if ((lastNonZeroIndex == -1) || !current_delta_nonzero) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  LONG  HISTORY  IP  PREFETCHER  IF  NO  PREFETCH  YET                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (longHistIPTableNewDelta) {
         int j, length, chosen_j=-1;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    DETERMINE STRIDE PATTERN MATCH WITH MATCH_LENGTH HISTORY                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         for (j=0; j<NUM_STRIDES_IN_LONG_HIST_IP_TABLE+1-LONG_HIST_MATCH_LENGTH; j++) {
            length = 0;
            while (length != LONG_HIST_MATCH_LENGTH) {
               if (longHistory[NUM_STRIDES_IN_LONG_HIST_IP_TABLE+1-LONG_HIST_MATCH_LENGTH+length] != longHistory[j+length]) break;
               length++;
            }
            if (length == LONG_HIST_MATCH_LENGTH) {
               assert(chosen_j == -1);
               chosen_j = j;
               break;
            }
         }
         if (chosen_j != -1) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MATCHING LONGEST PATTERN DETERMINED                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            j = chosen_j + LONG_HIST_MATCH_LENGTH;
            if (throttle_level_L1 == 0) {
               while (j < NUM_STRIDES_IN_LONG_HIST_IP_TABLE+1) {
                  uint64_t pf_address = (cl_addr + longHistory[j]) << LOG2_BLOCK_SIZE;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if (!RecentAccessTagArrayL1DetermineHit(pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                     if (prefetch_line(ip, addr, pf_address, FILL_L1, 0)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                           INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                        RecentAccessTagArrayL1Insert(pf_address >> LOG2_BLOCK_SIZE);
                        did_pref = true;
                     }
                  }
                  j++;
               }
               if (did_pref) return;
            }
         }
      }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  NEXT  LINE  PREFETCHER  IF  IP-DELTA  TABLE  CANNOT  OFFER  PREDICTION                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      uint64_t pf_address = (cl_addr+1) << LOG2_BLOCK_SIZE;
      if (throttle_level_L1 == 0) {
         i = 0;
         while (i < BASE_PREFETCH_DEGREE_L1) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INJECT  PREFETCHES  FOR  APPROPRIATE  DEGREE                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (degreeHitsL1[i]*NL_THRESHOLD_DENOM_L1 > degreeInsertionsL1[i]*NL_THRESHOLD_NUMER_L1) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               if (!RecentAccessTagArrayL1DetermineHit (pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if (prefetch_line(ip, addr, pf_address, FILL_L1, 0)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                           INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                     RecentAccessTagArrayL1Insert (pf_address >> LOG2_BLOCK_SIZE);
                  }
               }
            }
            pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
            i++;
         }
      }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INSERT  POSSIBLE  NEXT  LINE  PREFETCH  CANDIDATES  IN  THE  NL  BUFFER  TO  MEASURE  USEFULNESS //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      i = 0;
      pf_address = (cl_addr + 1) << LOG2_BLOCK_SIZE;
      while (i < BASE_PREFETCH_DEGREE_L1) {
         if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;
         NLBufferL1Insert (pf_address >> LOG2_BLOCK_SIZE, i);
         pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
         i++;
      }
      return;
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     IP  DELTA  TABLE  LOOKUP  TO  DECIDE  PREFETCH  CANDIDATES                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   ipPrefetchStride[0] = ipTableStride[lastNonZeroIndex];
   for (i=1; i<BASE_PREFETCH_DEGREE_L1+1; i++) ipPrefetchStride[i] = 0;
   unsigned delta = (ipTableStride[lastNonZeroIndex] >= 0) ? ipTableStride[lastNonZeroIndex] : ((-ipTableStride[lastNonZeroIndex]) | (1 <<  (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)));
   int ipDeltaTableIndex = (((pageid) << 3) ^ delta) & (NUM_SETS_IN_IP_DELTA_TABLE_L1 - 1);
   uint64_t ipDeltaTableTag = ((((pageid) << 3) ^ delta) >> LOG_NUM_SETS_IN_IP_DELTA_TABLE_L1) & IP_DELTA_TABLE_TAG_MASK;
   for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L1; ii++) {
      if (ipDeltaTableL1[ipDeltaTableIndex][ii].valid && (ipDeltaTableL1[ipDeltaTableIndex][ii].tag == ipDeltaTableTag)) {
         for (i=1; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
	    if (NUM_CPUS == 1) {
               if (((i < BASE_PREFETCH_DEGREE_L1) && (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= PREDICTION_THRESHOLD_L1)) ||
                   ((i == BASE_PREFETCH_DEGREE_L1) && (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= (PREDICTION_THRESHOLD_L1 + 1)))) {
                  ipPrefetchStride[i] = ipDeltaTableL1[ipDeltaTableIndex][ii].stride[i-1];
               }
            }
            else {
               if (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= PREDICTION_THRESHOLD_L1) {
                  ipPrefetchStride[i] = ipDeltaTableL1[ipDeltaTableIndex][ii].stride[i-1];
               }
            }
         }
         break;
      }
   }
   if (ii < NUM_WAYS_IN_IP_DELTA_TABLE_L1) {
      for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L1; j++) ipDeltaTableL1[ipDeltaTableIndex][j].lru++;
      ipDeltaTableL1[ipDeltaTableIndex][ii].lru = 0;
      ipDeltaTableL1[ipDeltaTableIndex][ii].partial_ip = ip & PARTIAL_IP_MASK;
      ipDeltaTableL1[ipDeltaTableIndex][ii].partial_ip_valid = true;
   }
   else {
      for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L1; ii++) {
         if (ipDeltaTableL1[ipDeltaTableIndex][ii].valid && ipDeltaTableL1[ipDeltaTableIndex][ii].partial_ip_valid && (ipDeltaTableL1[ipDeltaTableIndex][ii].partial_ip == (ip & PARTIAL_IP_MASK))) {
            for (i=1; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
               if (NUM_CPUS == 1) {
                  if (((i < BASE_PREFETCH_DEGREE_L1) && (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= PREDICTION_THRESHOLD_L1)) ||
                      ((i == BASE_PREFETCH_DEGREE_L1) && (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= (PREDICTION_THRESHOLD_L1 + 1)))) {
                     ipPrefetchStride[i] = ipDeltaTableL1[ipDeltaTableIndex][ii].stride[i-1];
                  }
               }
               else {
                  if (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= PREDICTION_THRESHOLD_L1) {
                     ipPrefetchStride[i] = ipDeltaTableL1[ipDeltaTableIndex][ii].stride[i-1];
                  }
               }
            }
            break;
         }
      }
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     INJECT  PREFETCHES                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   uint64_t pf_address = cl_addr << LOG2_BLOCK_SIZE;
   bool stopPrefetching = false;
   int num_pref = 0;
   if (throttle_level_L1 < THROTTLE_LEVEL_MAX_L1) {
      for (i=1; i<((throttle_level_L1 > 2) ? BASE_PREFETCH_DEGREE_L1 - (throttle_level_L1 - 2) + 1 : BASE_PREFETCH_DEGREE_L1 + 1); i++) {
         if (ipPrefetchStride[i] == 0) break;
         pf_address = pf_address + (ipPrefetchStride[i] << LOG2_BLOCK_SIZE);
         if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) continue;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (!RecentAccessTagArrayL1DetermineHit (pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (PQ.occupancy < (PQ.SIZE - 1)) {
               assert(prefetch_line(ip, addr, pf_address, FILL_L1, 0));
               did_pref = true;
            }
            else if (PQ.occupancy < PQ.SIZE) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  REMAINING  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               uint32_t pfmetadata = 0;
               unsigned char residue = 0;
               for (int j=i+1; j<((throttle_level_L1 > 2) ? BASE_PREFETCH_DEGREE_L1 - (throttle_level_L1 - 2) + 1 : BASE_PREFETCH_DEGREE_L1 + 1); j++) {
                  if (ipPrefetchStride[j] == 0) break;
                  unsigned char delta = ((ipPrefetchStride[j] < 0) ? ((-ipPrefetchStride[j]) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : ipPrefetchStride[j]);
                  pfmetadata = pfmetadata | (delta << ((1 + LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)*residue));
                  residue++;
               }
               assert(prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata));
               did_pref = true;
            }
            else {
               assert(!prefetch_line(ip, addr, pf_address, FILL_L1, 0));
               stopPrefetching = true;
               break;
            }
            num_pref++;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                           INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            RecentAccessTagArrayL1Insert (pf_address >> LOG2_BLOCK_SIZE);
         }
      }
   }

   if (throttle_level_L1 < 2) {
      if ((BASE_PREFETCH_DEGREE_L1 > 1) && !stopPrefetching && (i<BASE_PREFETCH_DEGREE_L1+1)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//            INJECT  MORE  PREFETCHES  BASED  ON  LAST  TWO  IDENTICAL  DELTAS,  IF  ANY              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         if ((ipPrefetchStride[i-1] != 0) && (i > 1) && (ipPrefetchStride[i-1] == ipPrefetchStride[i-2])) {
            assert(num_pref < BASE_PREFETCH_DEGREE_L1);

            while (1) {
               pf_address = pf_address + (ipPrefetchStride[i-1] << LOG2_BLOCK_SIZE);
               if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               if (!RecentAccessTagArrayL1DetermineHit(pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//             MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if (PQ.occupancy < (PQ.SIZE - 1)) {
                     assert(prefetch_line(ip, addr, pf_address, FILL_L1, 0));
                     num_pref++;
                     did_pref = true;
                  }
                  else if (PQ.occupancy < PQ.SIZE) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  REMAINING  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                     uint32_t pfmetadata = 0x80000000U;
                     unsigned char delta = ((ipPrefetchStride[i-1] < 0) ? ((-ipPrefetchStride[i-1]) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : ipPrefetchStride[i-1]);
                     pfmetadata = pfmetadata | delta;
                     assert(prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata));
                     num_pref++;
                     did_pref = true;
                  }
                  else {
                     assert(!prefetch_line(ip, addr, pf_address, FILL_L1, 0));
                     break;
                  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                         INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  RecentAccessTagArrayL1Insert (pf_address >> LOG2_BLOCK_SIZE);
               }

               if (num_pref == (BASE_PREFETCH_DEGREE_L1)) break;
            }
            if ((num_pref == (BASE_PREFETCH_DEGREE_L1)) && (PQ.occupancy < PQ.SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  FURTHER  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               pf_address = pf_address + (ipPrefetchStride[i-1] << LOG2_BLOCK_SIZE);
               if ((pf_address >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)) {
                  uint32_t pfmetadata = 0x80000000U;
                  unsigned char delta = ((ipPrefetchStride[i-1] < 0) ? ((-ipPrefetchStride[i-1]) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : ipPrefetchStride[i-1]);
                  pfmetadata = pfmetadata | delta;
                  assert(prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata));
                  did_pref = true;
               }
            }
         }
         else if (constantStrideValid) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                   INJECT  MORE  PREFETCHES  BASED  ON  IP  STRIDE                                   //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            while (1) {
               pf_address = pf_address + (constantStride << LOG2_BLOCK_SIZE);
               if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               if (!RecentAccessTagArrayL1DetermineHit(pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//              MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if (PQ.occupancy < (PQ.SIZE - 1)) {
                     assert(prefetch_line(ip, addr, pf_address, FILL_L1, 0));
                     num_pref++;
                     did_pref = true;
                  }
                  else if (PQ.occupancy < PQ.SIZE) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  REMAINING  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                     uint32_t pfmetadata = 0x80000000U;
                     unsigned char delta = ((constantStride < 0) ? ((-constantStride) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : constantStride);
                     pfmetadata = pfmetadata | delta;
                     assert(prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata));
                     num_pref++;
                     did_pref = true;
                  }
                  else {
                     assert(!prefetch_line(ip, addr, pf_address, FILL_L1, 0));
                     break;
                  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                       INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  RecentAccessTagArrayL1Insert(pf_address >> LOG2_BLOCK_SIZE);
               }

               if (num_pref == (BASE_PREFETCH_DEGREE_L1)) break;
            }
            if ((num_pref == (BASE_PREFETCH_DEGREE_L1)) && (PQ.occupancy < PQ.SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  FURTHER  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               pf_address = pf_address + (constantStride << LOG2_BLOCK_SIZE);
               if ((pf_address >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)) {
                  uint32_t pfmetadata = 0x80000000U;
                  unsigned char delta = ((constantStride < 0) ? ((-constantStride) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : constantStride);
                  pfmetadata = pfmetadata | delta;
                  assert(prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata));
                  did_pref = true;
               }
            }
         }
      }
      else if (!stopPrefetching) {
         if ((ipPrefetchStride[i-1] != 0) && (i > 1) && (ipPrefetchStride[i-1] == ipPrefetchStride[i-2])) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  REMAINING  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            pf_address = pf_address + (ipPrefetchStride[i-1] << LOG2_BLOCK_SIZE);
            if ((pf_address >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)) {
               uint32_t pfmetadata = 0x80000000U;
               unsigned char delta = ((ipPrefetchStride[i-1] < 0) ? ((-ipPrefetchStride[i-1]) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : ipPrefetchStride[i-1]);
               pfmetadata = pfmetadata | delta;
               prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata);
               did_pref = true;
            }
         }
         else if (constantStrideValid) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  REMAINING  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            pf_address = pf_address + (constantStride << LOG2_BLOCK_SIZE);
            if ((pf_address >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)) {
               uint32_t pfmetadata = 0x80000000U;
               unsigned char delta = ((constantStride < 0) ? ((-constantStride) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : constantStride);
               pfmetadata = pfmetadata | delta;
               prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata);
               did_pref = true;
            }
         }
      }
   }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  LONG  HISTORY  IP  PREFETCHER  IF  NO  PREFETCH  YET                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (!did_pref && longHistIPTableNewDelta) {
      int j, length, chosen_j=-1;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    DETERMINE STRIDE PATTERN MATCH WITH MATCH_LENGTH HISTORY                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (j=0; j<NUM_STRIDES_IN_LONG_HIST_IP_TABLE+1-LONG_HIST_MATCH_LENGTH; j++) {
         length = 0;
         while (length != LONG_HIST_MATCH_LENGTH) {
            if (longHistory[NUM_STRIDES_IN_LONG_HIST_IP_TABLE+1-LONG_HIST_MATCH_LENGTH+length] != longHistory[j+length]) break;
            length++;
	 }
         if (length == LONG_HIST_MATCH_LENGTH) {
	    assert(chosen_j == -1);
            chosen_j = j;
            break;
         }
      }
      if (chosen_j != -1) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MATCHING LONGEST PATTERN DETERMINED                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         j = chosen_j + LONG_HIST_MATCH_LENGTH;
         if (throttle_level_L1 == 0) {
            while (j < NUM_STRIDES_IN_LONG_HIST_IP_TABLE+1) {
               uint64_t pf_address = (cl_addr + longHistory[j]) << LOG2_BLOCK_SIZE;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               if (!RecentAccessTagArrayL1DetermineHit(pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if (prefetch_line(ip, addr, pf_address, FILL_L1, 0)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                           INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                     RecentAccessTagArrayL1Insert(pf_address >> LOG2_BLOCK_SIZE);
                     did_pref = true;
                  }
               }
               j++;
            }
         }
      }
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  NEXT  LINE  PREFETCHER  IF  NO  PREFETCH  YET                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (!did_pref) {
      uint64_t pf_address = (cl_addr+1) << LOG2_BLOCK_SIZE;
      if (throttle_level_L1 == 0) {
         i = 0;
         while (i < BASE_PREFETCH_DEGREE_L1) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INJECT  PREFETCHES  FOR  APPROPRIATE  DEGREE                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (degreeHitsL1[i]*NL_THRESHOLD_DENOM_L1 > degreeInsertionsL1[i]*NL_THRESHOLD_NUMER_L1) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECENT  ACCESS  TAG  ARRAY  LOOKUP  TO  DETERMINE  RECENT  ACCESS  TO  PREFETCH  CANDIDATE       //
////////////////////////////////////////////////////////////////////////////////////////////////////////
               if (!RecentAccessTagArrayL1DetermineHit(pf_address >> LOG2_BLOCK_SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MISS  IN  RECENT  ACCESS  TAG  ARRAY, INJECT  PREFETCH                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if (prefetch_line(ip, addr, pf_address, FILL_L1, 0)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                           INSERT  IN  RECENT  ACCESS  TAG  ARRAY                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////
                     RecentAccessTagArrayL1Insert(pf_address >> LOG2_BLOCK_SIZE);
                  }
               }
            }
            pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
            i++;
         }
      }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INSERT  POSSIBLE  NEXT  LINE  PREFETCH  CANDIDATES  INTO  NL  BUFFER TO  MEASURE  USEFULNESS     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      i = 0;
      pf_address = (cl_addr+1) << LOG2_BLOCK_SIZE;
      while (i < BASE_PREFETCH_DEGREE_L1) {
         if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;
         NLBufferL1Insert(pf_address >> LOG2_BLOCK_SIZE, i);
         pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
         i++;
      }
   }

    // NOTE: The following if condition is eternally useless

    // if no prefetches are issued till now, speculatively issue a next_line prefetch
    // if(num_prefs == 0 && spec_nl[cpu] == 1){                                        // NL IP
    //     // uint64_t pf_address = ((addr>>LOG2_BLOCK_SIZE)+1) << LOG2_BLOCK_SIZE;  
    //     // metadata = encode_metadata(1, NL_TYPE, spec_nl[cpu]);
    //     // prefetch_line(ip, addr, pf_address, FILL_L1, metadata);
    //     // SIG_DP(cout << "1, ");
    

    SIG_DP(cout << endl);

    return;
}

void CACHE::l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{

}
void CACHE::l1d_prefetcher_final_stats()
{
    cout << endl;
}

