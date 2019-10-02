//
// Written for Data Prefetching Championship Simulator 3
// Mainak Chaudhuri, Nayan Deshmukh. Sangam: A Multi-component Core Cache Prefetcher.
//

#include "cache.h"

#define LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 0
#define NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1 (1 << LOG_NUM_SETS_IN_RECENT_ACCESS_TAG_ARRAY_L1)
#define NUM_WAYS_IN_RECENT_ACCESS_TAG_ARRAY_L1 40

#define IP_TABLE_TAG_MASK 0x1ffffULL                 // For prefetch degree 4
#define LOG_NUM_SETS_IN_IP_TABLE_L1 7
#define NUM_SETS_IN_IP_TABLE_L1 (1 << LOG_NUM_SETS_IN_IP_TABLE_L1)
#define NUM_WAYS_IN_IP_TABLE_L1 15

#define IP_DELTA_TABLE_TAG_MASK 0xffffffULL        // For prefetch degree 4
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

// IP stride table entry: 63 bits
class IPtableL1 {
   public:
      // Tag
      uint64_t tag;					// 63 - (4+1+6+7*(BASE_PREFETCH_DEGREE_L1+1)) bits

      // LRU states
      uint64_t lru;					// 4 bits

      // Valid bit
      bool valid;					// 1 bit

      // Last seen offset in a page
      unsigned char offset;				// 6 bits

      // Last seen strides within a page
      char stride[BASE_PREFETCH_DEGREE_L1+1];		// 7*(BASE_PREFETCH_DEGREE_L1+1)

      IPtableL1 () {}
};

IPtableL1 ipTableL1[NUM_SETS_IN_IP_TABLE_L1][NUM_WAYS_IN_IP_TABLE_L1];  // 128x15x63 bits = 120960 bits
char ipTableStride[BASE_PREFETCH_DEGREE_L1+1];                               // 7*5 bits = 35 bits

// IP-delta correlation table entry: 64 bits
class IPDeltaTableL1 {
   public:
      // Tag
      uint64_t tag;					// 64 - (3+1+9*BASE_PREFETCH_DEGREE_L1) bits

      // LRU states
      uint64_t lru;					// 3 bits

      // Valid bit
      bool valid;					// 1 bit

      // Last seen strides within a page
      char stride[BASE_PREFETCH_DEGREE_L1];			// 7*BASE_PREFETCH_DEGREE_L1 bits

      // Confidence counters
      unsigned char counters[BASE_PREFETCH_DEGREE_L1];	// 2*BASE_PREFETCH_DEGREE_L1 bits

      IPDeltaTableL1 () {}
};

char ipPrefetchStride[BASE_PREFETCH_DEGREE_L1+1];   // 7*5 bits = 35 bits
IPDeltaTableL1 ipDeltaTableL1[NUM_SETS_IN_IP_DELTA_TABLE_L1][NUM_WAYS_IN_IP_DELTA_TABLE_L1];  // 256x8x64 bits = 131072 bits

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  TOTAL L1 CACHE PREFETCHER STORAGE BUDGET =                                                              //
//  4672+128+128+2840+120960+35+35+131072 BITS = 259870 BITS                                                //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CACHE::l1d_prefetcher_initialize() 
{
   cout << "CPU " << cpu << " L1d IP-delta,IP stride,NL prefetcher" << endl;

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
}

void CACHE::l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type)
{
   int i;

   if (type == PREFETCH) return;	// Don't update or inject new prefetches on prefetch events

   uint64_t cl_addr = addr >> LOG2_BLOCK_SIZE;
   unsigned char offset = (addr & PAGE_OFFSET_MASK) >> LOG2_BLOCK_SIZE;

   bool did_pref = false;
   bool current_delta_nonzero = false;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     NL  BUFFER  LOOKUP                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   for (i=0; i<NUM_ENTRIES_IN_NL_BUFFER_L1; i++) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    DETERMINE  NEXT  LINE  PREFETCHER'S  APPROPRIATE  DEGREE  BY MEASURING  USEFULNESS              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (nlBufferL1[i].valid && (nlBufferL1[i].tag == cl_addr)) {
         degreeHitsL1[nlBufferL1[i].degree-1]++;
         nlBufferL1[i].valid = false;
         break;
      }
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     IP  TABLE  LOOKUP                                                               //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   int ipTableIndex = (ip) & (NUM_SETS_IN_IP_TABLE_L1 - 1);
   uint64_t ipTableTag = ((ip) >> LOG_NUM_SETS_IN_IP_TABLE_L1) & IP_TABLE_TAG_MASK;
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
      for (i=0; i<BASE_PREFETCH_DEGREE_L1+1; i++) ipTableL1[ipTableIndex][ii].stride[i] = 0;
      ipTableL1[ipTableIndex][ii].valid = true;
   }
   for (i=0; i<NUM_WAYS_IN_IP_TABLE_L1; i++) ipTableL1[ipTableIndex][i].lru++;
   ipTableL1[ipTableIndex][ii].lru = 0;
   int lastNonZeroIndex = -1;
   for (i=0; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
      ipTableStride[i] = ipTableL1[ipTableIndex][ii].stride[i];
      if (ipTableStride[i] != 0) lastNonZeroIndex = i;
   }
   char fallBackPrefetchStride;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECORD  IF  LAST  TWO  NON-ZERO  STRIDES  ARE  EQUAL; WILL  BE  USED  FOR  FALL-BACK  STRIDE  PREFETCHING  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if ((lastNonZeroIndex > 0) && (ipTableStride[lastNonZeroIndex] == ipTableStride[lastNonZeroIndex-1])) fallBackPrefetchStride = ipTableStride[lastNonZeroIndex];
   else fallBackPrefetchStride = 0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     IP  DELTA  TABLE  LOOKUP  AND  TRAINING                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (current_delta_nonzero) {
      for (i=0; i<lastNonZeroIndex; i++) {
         assert(ipTableStride[i] != 0);
         unsigned delta;
         delta = (ipTableStride[i] >= 0) ? ipTableStride[i] : ((-ipTableStride[i]) | (1 <<  (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)));
         int ipDeltaTableIndex = ((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) & (NUM_SETS_IN_IP_DELTA_TABLE_L1 - 1);
         uint64_t ipDeltaTableTag = (((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) >> LOG_NUM_SETS_IN_IP_DELTA_TABLE_L1) & IP_DELTA_TABLE_TAG_MASK;
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     RECENT  ACCESS  TAG  ARRAY  LOOKUP  AND  INSERT                                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   RecentAccessTagArrayL1LookupAndInsertIfMiss (cl_addr);

   if ((lastNonZeroIndex == -1) || !current_delta_nonzero) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  NEXT  LINE  PREFETCHER  IF  IP-DELTA  TABLE  CANNOT  OFFER  PREDICTION                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      uint64_t pf_address = (cl_addr+1) << LOG2_BLOCK_SIZE;
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
   int ipDeltaTableIndex = ((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) & (NUM_SETS_IN_IP_DELTA_TABLE_L1 - 1);
   uint64_t ipDeltaTableTag = (((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) >> LOG_NUM_SETS_IN_IP_DELTA_TABLE_L1) & IP_DELTA_TABLE_TAG_MASK;
   for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L1; ii++) {
      if (ipDeltaTableL1[ipDeltaTableIndex][ii].valid && (ipDeltaTableL1[ipDeltaTableIndex][ii].tag == ipDeltaTableTag)) {
         for (i=1; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
            if (ipDeltaTableL1[ipDeltaTableIndex][ii].counters[i-1] >= PREDICTION_THRESHOLD_L1) {
               ipPrefetchStride[i] = ipDeltaTableL1[ipDeltaTableIndex][ii].stride[i-1];
            }
         }
         break;
      }
   }
   if (ii < NUM_WAYS_IN_IP_DELTA_TABLE_L1) {
      for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L1; j++) ipDeltaTableL1[ipDeltaTableIndex][j].lru++;
      ipDeltaTableL1[ipDeltaTableIndex][ii].lru = 0;
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     INJECT  PREFETCHES                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   uint64_t pf_address = cl_addr << LOG2_BLOCK_SIZE;
   bool stopPrefetching = false;
   int num_pref = 0;
   for (i=1; i<BASE_PREFETCH_DEGREE_L1+1; i++) {
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
            for (int j=i+1; j<BASE_PREFETCH_DEGREE_L1+1; j++) {
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

            if (num_pref == BASE_PREFETCH_DEGREE_L1) break;
         }
         if ((num_pref == BASE_PREFETCH_DEGREE_L1) && (PQ.occupancy < PQ.SIZE)) {
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
      else if (fallBackPrefetchStride != 0) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                   INJECT  MORE  PREFETCHES  BASED  ON  IP  STRIDE                                   //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         while (1) {
            pf_address = pf_address + (fallBackPrefetchStride << LOG2_BLOCK_SIZE);
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
                  unsigned char delta = ((fallBackPrefetchStride < 0) ? ((-fallBackPrefetchStride) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : fallBackPrefetchStride);
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

            if (num_pref == BASE_PREFETCH_DEGREE_L1) break;
         }
         if ((num_pref == BASE_PREFETCH_DEGREE_L1) && (PQ.occupancy < PQ.SIZE)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  FURTHER  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
            pf_address = pf_address + (fallBackPrefetchStride << LOG2_BLOCK_SIZE);
            if ((pf_address >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)) {
               uint32_t pfmetadata = 0x80000000U;
               unsigned char delta = ((fallBackPrefetchStride < 0) ? ((-fallBackPrefetchStride) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : fallBackPrefetchStride);
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
      else if (fallBackPrefetchStride != 0) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    HAND  OVER  REMAINING  PREFETCH  INJECTION  TO  L2  CACHE  BY  ENCODING  THE  METADATA           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         pf_address = pf_address + (fallBackPrefetchStride << LOG2_BLOCK_SIZE);
         if ((pf_address >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)) {
            uint32_t pfmetadata = 0x80000000U;
            unsigned char delta = ((fallBackPrefetchStride < 0) ? ((-fallBackPrefetchStride) | (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE))) : fallBackPrefetchStride);
            pfmetadata = pfmetadata | delta;
            prefetch_line(ip, addr, pf_address, FILL_L1, pfmetadata);
            did_pref = true;
         }
      }
   }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  NEXT  LINE  PREFETCHER  IF  NO  PREFETCH  YET                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (!did_pref) {
      uint64_t pf_address = (cl_addr+1) << LOG2_BLOCK_SIZE;
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
}

void CACHE::l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
}

void CACHE::l1d_prefetcher_final_stats()
{
}
