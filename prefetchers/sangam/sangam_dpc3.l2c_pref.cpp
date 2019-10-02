//
// Written for Data Prefetching Championship Simulator 3
// Mainak Chaudhuri, Nayan Deshmukh. Sangam: A Multi-component Core Cache Prefetcher.
//

#include "cache.h"

#define BASE_PREFETCH_DEGREE_L2 ((NUM_CPUS == 1) ? 3 : 2)

#define IP_TABLE_TAG_MASK_L2 ((NUM_CPUS == 1) ? 0xffffffULL : 0x7fffffffULL)   // For prefetch degree 3 or 2
#define LOG_NUM_SETS_IN_IP_TABLE_L2 7
#define NUM_SETS_IN_IP_TABLE_L2 (1 << LOG_NUM_SETS_IN_IP_TABLE_L2)
#define NUM_WAYS_IN_IP_TABLE_L2 15

#define IP_DELTA_TABLE_TAG_MASK_L2 ((NUM_CPUS == 1) ? 0x1ffffffffULL : 0x3ffffffffffULL)   // For prefetch degree 3 or 2
#define LOG_NUM_SETS_IN_IP_DELTA_TABLE_L2 8
#define NUM_SETS_IN_IP_DELTA_TABLE_L2 (1 << LOG_NUM_SETS_IN_IP_DELTA_TABLE_L2)
#define NUM_WAYS_IN_IP_DELTA_TABLE_L2 8

#define SATURATING_COUNTER_MAX_L2 3
#define PREDICTION_THRESHOLD_L2 2

#define PAGE_OFFSET_MASK_L2 ((1 << LOG2_PAGE_SIZE) - 1)

#define NUM_ENTRIES_IN_NL_BUFFER_L2 64
#define NL_THRESHOLD_NUMER_L2 1
#define NL_THRESHOLD_DENOM_L2 2

// Next line buffer entry: 73 bits
class NLBufferL2 {
   public:
      // Tag
      uint64_t tag;   // 64 bits

      // LRU states
      uint64_t lru;   // 6 bits

      // Valid bit
      bool valid;     // 1 bit

      // Degree
      unsigned char degree;  // 2 bits

      NLBufferL2 () {}
};

NLBufferL2 nlBufferL2[NUM_ENTRIES_IN_NL_BUFFER_L2];           // 64x73 bits = 4672 bits
unsigned degreeInsertionsL2[BASE_PREFETCH_DEGREE_L2]; // 32x3 bits = 96 bits
unsigned degreeHitsL2[BASE_PREFETCH_DEGREE_L2];       // 32x3 bits = 96 bits

static void NLBufferL2Insert (uint64_t cl_addr, int current_degree_index)
{
   int j;
   for (j=0; j<NUM_ENTRIES_IN_NL_BUFFER_L2; j++) {
      if (nlBufferL2[j].valid && (nlBufferL2[j].tag == cl_addr)) {
         if (nlBufferL2[j].degree > (current_degree_index+1)) {
            assert(degreeInsertionsL2[nlBufferL2[j].degree-1] > 0);
            degreeInsertionsL2[nlBufferL2[j].degree-1]--;
            nlBufferL2[j].degree = current_degree_index+1;   // Always favor smaller degree NL prefetcher
            degreeInsertionsL2[current_degree_index]++;
            for (int jj=0; jj<NUM_ENTRIES_IN_NL_BUFFER_L2; jj++) {
               nlBufferL2[jj].lru++;
            }
            nlBufferL2[j].lru = 0;
         }
         break;
      }
   }
   if (j == NUM_ENTRIES_IN_NL_BUFFER_L2) {  // MISS
      for (j=0; j<NUM_ENTRIES_IN_NL_BUFFER_L2; j++) {
         if (!nlBufferL2[j].valid) break;
      }
      if (j == NUM_ENTRIES_IN_NL_BUFFER_L2) {
         uint64_t maxlru = 0;
         int repl_index = -1;
         for (j=0; j<NUM_ENTRIES_IN_NL_BUFFER_L2; j++) {
            if (nlBufferL2[j].lru >= maxlru) {
               maxlru = nlBufferL2[j].lru;
               repl_index = j;
            }
         }
         j = repl_index;
      }
      nlBufferL2[j].tag = cl_addr;
      nlBufferL2[j].degree = current_degree_index+1;
      nlBufferL2[j].valid = true;
      degreeInsertionsL2[current_degree_index]++;
      for (int jj=0; jj<NUM_ENTRIES_IN_NL_BUFFER_L2; jj++) {
         nlBufferL2[jj].lru++;
      }
      nlBufferL2[j].lru = 0;
   }
}

// IP stride table entry: 63 bits
class IPtableL2 {
   public:
      // Tag
      uint64_t tag;                                     // 63 - (4+1+6+7*(BASE_PREFETCH_DEGREE_L2+1)) bits

      // LRU states
      uint64_t lru;                                     // 4 bits

      // Valid bit
      bool valid;                                       // 1 bit

      // Last seen offset in a page
      unsigned char offset;                             // 6 bits

      // Last seen strides within a page
      char stride[BASE_PREFETCH_DEGREE_L2+1];           // 7*(BASE_PREFETCH_DEGREE_L2+1) bits

      IPtableL2 () {}
};

IPtableL2 ipTableL2[NUM_SETS_IN_IP_TABLE_L2][NUM_WAYS_IN_IP_TABLE_L2];  // 128x15x63 bits = 120960 bits
char ipTableStrideL2[BASE_PREFETCH_DEGREE_L2+1];                        // 7x4 bits = 28 bits

// IP-delta correlation table entry: 64 bits
class IPDeltaTableL2 {
   public:
      // Tag
      uint64_t tag;                                     // 64 - (3+1+9*BASE_PREFETCH_DEGREE_L2) bits

      // LRU states
      uint64_t lru;                                     // 3 bits

      // Valid bit
      bool valid;                                       // 1 bit

      // Last seen strides within a page
      char stride[BASE_PREFETCH_DEGREE_L2];                  // 7*BASE_PREFETCH_DEGREE_L2 bits

      // Confidence counters
      unsigned char counters[BASE_PREFETCH_DEGREE_L2];       // 2*BASE_PREFETCH_DEGREE_L2 bits

      IPDeltaTableL2 () {}
};

char ipPrefetchStrideL2[BASE_PREFETCH_DEGREE_L2+1];   // 7x4 bits = 28 bits
IPDeltaTableL2 ipDeltaTableL2[NUM_SETS_IN_IP_DELTA_TABLE_L2][NUM_WAYS_IN_IP_DELTA_TABLE_L2];  // 256x8x64 bits = 131072 bits

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  TOTAL L2 CACHE PREFETCHER STORAGE BUDGET =                                                              //
//  4672+96+96+120960+28+28+131072 BITS = 256952 BITS                                                       //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CACHE::l2c_prefetcher_initialize() 
{
    cout << "CPU " << cpu << " L2C IP-delta,IP stride,NL prefetcher" << endl;

    for (int i=0; i<NUM_ENTRIES_IN_NL_BUFFER_L2; i++) {
      nlBufferL2[i].valid = false;
      nlBufferL2[i].lru = 0;
    }

    for (int i=0; i<BASE_PREFETCH_DEGREE_L2; i++) {
      degreeInsertionsL2[i] = 0;
      degreeHitsL2[i] = 0;
    }

    for (int i=0; i<NUM_SETS_IN_IP_TABLE_L2; i++) {
      for (int j=0; j<NUM_WAYS_IN_IP_TABLE_L2; j++) {
         ipTableL2[i][j].valid = false;
         ipTableL2[i][j].lru = 0;
      }
   }

   for (int i=0; i<NUM_SETS_IN_IP_DELTA_TABLE_L2; i++) {
      for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L2; j++) {
         ipDeltaTableL2[i][j].valid = false;
         ipDeltaTableL2[i][j].lru = 0;
      }
   }
}

uint32_t CACHE::l2c_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{
    if (type == PREFETCH) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    FINISH  OFF  RESIDUAL  PREFETCH  INJECTION  AS ENCODED  BY  THE  L2  CACHE  PREFETCHER           //
////////////////////////////////////////////////////////////////////////////////////////////////////////
       if (metadata_in != 0) {
          uint64_t pf_addr = (addr>>LOG2_BLOCK_SIZE)<< LOG2_BLOCK_SIZE;   // Aligned base
          if ((metadata_in & 0x80000000U) != 0) {  // Repeated stride encoding
             char delta = metadata_in & ((1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)) - 1);
             assert(delta != 0);
             unsigned char signbit = metadata_in & (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE));
             if (signbit) delta = -delta;
             int i;
             for (i=0; i<BASE_PREFETCH_DEGREE_L2; i++) {
                pf_addr = pf_addr + (delta<<LOG2_BLOCK_SIZE);
                prefetch_line(ip, addr, pf_addr, FILL_L2, 0);
             }
          }
          else {   // Delta sequence encoding
             int num_pref = 0;
             char lastdelta = 0;
             char lastlastdelta = 0;
             while (metadata_in != 0) {
                char delta = metadata_in & ((1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)) - 1);
                assert (delta != 0);
                unsigned char signbit = metadata_in & (1 << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE));
                if (signbit) delta = -delta;
                pf_addr = pf_addr + (delta<<LOG2_BLOCK_SIZE);
                if (prefetch_line(ip, addr, pf_addr, FILL_L2, 0)) num_pref++;
                metadata_in = metadata_in >> (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1);
                lastlastdelta = lastdelta;
                lastdelta = delta;
             }
             if ((metadata_in == 0) && (num_pref < BASE_PREFETCH_DEGREE_L2) && (lastlastdelta == lastdelta)) {
                int i;
                for (i=num_pref; i<BASE_PREFETCH_DEGREE_L2; i++) {
                   pf_addr = pf_addr + (lastdelta<<LOG2_BLOCK_SIZE);
                   prefetch_line(ip, addr, pf_addr, FILL_L2, 0);
                }
             }
          }
       }
    }

    int i;
    bool did_pref = false;
    bool current_delta_nonzero = false;

    uint64_t cl_addr = addr >> LOG2_BLOCK_SIZE;
    unsigned char offset = (addr & PAGE_OFFSET_MASK_L2) >> LOG2_BLOCK_SIZE;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     NL  BUFFER  LOOKUP                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   for (i=0; i<NUM_ENTRIES_IN_NL_BUFFER_L2; i++) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    DETERMINE  NEXT  LINE  PREFETCHER'S  APPROPRIATE  DEGREE  BY MEASURING  USEFULNESS               //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (nlBufferL2[i].valid && (nlBufferL2[i].tag == cl_addr)) {
         degreeHitsL2[nlBufferL2[i].degree-1]++;
         nlBufferL2[i].valid = false;
         break;
      }
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     IP  TABLE  LOOKUP                                                               //
////////////////////////////////////////////////////////////////////////////////////////////////////////
    int ipTableIndex = (ip) & (NUM_SETS_IN_IP_TABLE_L2 - 1);
    uint64_t ipTableTag = ((ip) >> LOG_NUM_SETS_IN_IP_TABLE_L2) & IP_TABLE_TAG_MASK_L2;
    int ii;
    for (ii=0; ii<NUM_WAYS_IN_IP_TABLE_L2; ii++) {
      if (ipTableL2[ipTableIndex][ii].valid && (ipTableL2[ipTableIndex][ii].tag == ipTableTag)) {
         if ((signed)(offset - ipTableL2[ipTableIndex][ii].offset) != 0) {
            current_delta_nonzero = true;
            for (i=0; i<BASE_PREFETCH_DEGREE_L2+1; i++) {
               if (ipTableL2[ipTableIndex][ii].stride[i] == 0) {
                  ipTableL2[ipTableIndex][ii].stride[i] = (signed)(offset - ipTableL2[ipTableIndex][ii].offset);
                  break;
               }
            }
            if (i == BASE_PREFETCH_DEGREE_L2+1) {
               for (i=0; i<BASE_PREFETCH_DEGREE_L2; i++) {
                  ipTableL2[ipTableIndex][ii].stride[i] = ipTableL2[ipTableIndex][ii].stride[i+1];
               }
               ipTableL2[ipTableIndex][ii].stride[i] = (signed)(offset - ipTableL2[ipTableIndex][ii].offset);
            }
            ipTableL2[ipTableIndex][ii].offset = offset;
         }
         break;
      }
   }
   if (ii == NUM_WAYS_IN_IP_TABLE_L2) {
      for (ii=0; ii<NUM_WAYS_IN_IP_TABLE_L2; ii++) {
         if (!ipTableL2[ipTableIndex][ii].valid) break;
      }
      if (ii == NUM_WAYS_IN_IP_TABLE_L2) {
         uint64_t maxlru = 0;
         int repl_index = -1;
         for (ii=0; ii<NUM_WAYS_IN_IP_TABLE_L2; ii++) {
            if (ipTableL2[ipTableIndex][ii].lru > maxlru) {
               maxlru = ipTableL2[ipTableIndex][ii].lru;
               repl_index = ii;
            }
         }
         ii = repl_index;
      }
      ipTableL2[ipTableIndex][ii].tag = ipTableTag;
      ipTableL2[ipTableIndex][ii].offset = offset;
      for (i=0; i<BASE_PREFETCH_DEGREE_L2+1; i++) ipTableL2[ipTableIndex][ii].stride[i] = 0;
      ipTableL2[ipTableIndex][ii].valid = true;
   }
   for (i=0; i<NUM_WAYS_IN_IP_TABLE_L2; i++) ipTableL2[ipTableIndex][i].lru++;
   ipTableL2[ipTableIndex][ii].lru = 0;
   int lastNonZeroIndex = -1;
   for (i=0; i<BASE_PREFETCH_DEGREE_L2+1; i++) {
      ipTableStrideL2[i] = ipTableL2[ipTableIndex][ii].stride[i];
      if (ipTableStrideL2[i] != 0) lastNonZeroIndex = i;
   }
   char fallBackPrefetchStride;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    RECORD  IF  LAST  TWO  NON-ZERO  STRIDES  ARE  EQUAL; WILL  BE  USED  FOR  FALL-BACK  STRIDE  PREFETCHING  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if ((lastNonZeroIndex > 0) && (ipTableStrideL2[lastNonZeroIndex] == ipTableStrideL2[lastNonZeroIndex-1])) fallBackPrefetchStride = ipTableStrideL2[lastNonZeroIndex];
   else fallBackPrefetchStride = 0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     IP  DELTA  TABLE  LOOKUP  AND  TRAINING                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (current_delta_nonzero) {
      for (i=0; i<lastNonZeroIndex; i++) {
         assert(ipTableStrideL2[i] != 0);
         unsigned delta;
         delta = (ipTableStrideL2[i] >= 0) ? ipTableStrideL2[i] : ((-ipTableStrideL2[i]) | (1 <<  (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)));
         int ipDeltaTableIndex = ((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) & (NUM_SETS_IN_IP_DELTA_TABLE_L2 - 1);
         uint64_t ipDeltaTableTag = (((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) >> LOG_NUM_SETS_IN_IP_DELTA_TABLE_L2) & IP_DELTA_TABLE_TAG_MASK_L2;
         for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L2; ii++) {
            if (ipDeltaTableL2[ipDeltaTableIndex][ii].valid && (ipDeltaTableL2[ipDeltaTableIndex][ii].tag == ipDeltaTableTag)) {
               if (ipTableStrideL2[lastNonZeroIndex] == ipDeltaTableL2[ipDeltaTableIndex][ii].stride[lastNonZeroIndex-i-1]) {
                  if (ipDeltaTableL2[ipDeltaTableIndex][ii].counters[lastNonZeroIndex-i-1] < SATURATING_COUNTER_MAX_L2) {
                     ipDeltaTableL2[ipDeltaTableIndex][ii].counters[lastNonZeroIndex-i-1]++;
                  }
               }
               else {
                  ipDeltaTableL2[ipDeltaTableIndex][ii].stride[lastNonZeroIndex-i-1] = ipTableStrideL2[lastNonZeroIndex];
                  ipDeltaTableL2[ipDeltaTableIndex][ii].counters[lastNonZeroIndex-i-1] = 1;
               }
               break;
            }
         }
         if (ii == NUM_WAYS_IN_IP_DELTA_TABLE_L2) {
            for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L2; ii++) {
               if (!ipDeltaTableL2[ipDeltaTableIndex][ii].valid) break;
            }
            if (ii == NUM_WAYS_IN_IP_DELTA_TABLE_L2) {
               uint64_t maxlru = 0;
               int repl_index = -1;
               for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L2; ii++) {
                  if (ipDeltaTableL2[ipDeltaTableIndex][ii].lru > maxlru) {
                     maxlru = ipDeltaTableL2[ipDeltaTableIndex][ii].lru;
                     repl_index = ii;
                  }
               }
               ii = repl_index;
            }
            ipDeltaTableL2[ipDeltaTableIndex][ii].tag = ipDeltaTableTag;
            for (int j=0; j<BASE_PREFETCH_DEGREE_L2; j++) {
               if (i+j+1 < BASE_PREFETCH_DEGREE_L2+1) {
                  ipDeltaTableL2[ipDeltaTableIndex][ii].stride[j] = ipTableStrideL2[i+j+1];
                  if (ipTableStrideL2[i+j+1] != 0) ipDeltaTableL2[ipDeltaTableIndex][ii].counters[j] = 1;
                  else ipDeltaTableL2[ipDeltaTableIndex][ii].counters[j] = 0;
               }
               else {
                  ipDeltaTableL2[ipDeltaTableIndex][ii].stride[j] = 0;
                  ipDeltaTableL2[ipDeltaTableIndex][ii].counters[j] = 0;
               }
            }
            ipDeltaTableL2[ipDeltaTableIndex][ii].valid = true;
         }
         for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L2; j++) ipDeltaTableL2[ipDeltaTableIndex][j].lru++;
         ipDeltaTableL2[ipDeltaTableIndex][ii].lru = 0;
      }
   }


   if ((lastNonZeroIndex == -1) || !current_delta_nonzero) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  NEXT  LINE  PREFETCHER  IF  IP-DELTA  TABLE  CANNOT  OFFER  PREDICTION                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      uint64_t pf_address = (cl_addr + 1) << LOG2_BLOCK_SIZE;
      i = 0;
      while (i < BASE_PREFETCH_DEGREE_L2) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INJECT  PREFETCHES  FOR  APPROPRIATE  DEGREE                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (degreeHitsL2[i]*NL_THRESHOLD_DENOM_L2 > degreeInsertionsL2[i]*NL_THRESHOLD_NUMER_L2) {
            prefetch_line(ip, addr, pf_address, FILL_L2, 0);
         }
         pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
         i++;
      }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INSERT  POSSIBLE  NEXT  LINE  PREFETCH  CANDIDATES  IN  THE  NL  BUFFER  TO  MEASURE  USEFULNESS //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      i = 0;
      pf_address = (cl_addr + 1) << LOG2_BLOCK_SIZE;
      while (i < BASE_PREFETCH_DEGREE_L2) {
         if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;
         NLBufferL2Insert(pf_address >> LOG2_BLOCK_SIZE, i);
         pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
         i++;
      }
      return 0;
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     IP  DELTA  TABLE  LOOKUP  TO  DECIDE  PREFETCH  CANDIDATES                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   ipPrefetchStrideL2[0] = ipTableStrideL2[lastNonZeroIndex];
   for (i=1; i<BASE_PREFETCH_DEGREE_L2+1; i++) ipPrefetchStrideL2[i] = 0;
   unsigned delta = (ipTableStrideL2[lastNonZeroIndex] >= 0) ? ipTableStrideL2[lastNonZeroIndex] : ((-ipTableStrideL2[lastNonZeroIndex]) | (1 <<  (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE)));
   int ipDeltaTableIndex = ((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) & (NUM_SETS_IN_IP_DELTA_TABLE_L2 - 1);
   uint64_t ipDeltaTableTag = (((ip << (LOG2_PAGE_SIZE - LOG2_BLOCK_SIZE + 1)) | delta) >> LOG_NUM_SETS_IN_IP_DELTA_TABLE_L2) & IP_DELTA_TABLE_TAG_MASK_L2;
   for (ii=0; ii<NUM_WAYS_IN_IP_DELTA_TABLE_L2; ii++) {
      if (ipDeltaTableL2[ipDeltaTableIndex][ii].valid && (ipDeltaTableL2[ipDeltaTableIndex][ii].tag == ipDeltaTableTag)) {
         for (i=1; i<BASE_PREFETCH_DEGREE_L2+1; i++) {
            if (ipDeltaTableL2[ipDeltaTableIndex][ii].counters[i-1] >= PREDICTION_THRESHOLD_L2) {
               ipPrefetchStrideL2[i] = ipDeltaTableL2[ipDeltaTableIndex][ii].stride[i-1];
            }
         }
         break;
      }
   }
   if (ii < NUM_WAYS_IN_IP_DELTA_TABLE_L2) {
      for (int j=0; j<NUM_WAYS_IN_IP_DELTA_TABLE_L2; j++) ipDeltaTableL2[ipDeltaTableIndex][j].lru++;
      ipDeltaTableL2[ipDeltaTableIndex][ii].lru = 0;
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                     INJECT  PREFETCHES                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   uint64_t pf_address = cl_addr << LOG2_BLOCK_SIZE;
   bool stopPrefetching = false;
   int num_pref = 0;
   for (i=1; i<BASE_PREFETCH_DEGREE_L2+1; i++) {
      if (ipPrefetchStrideL2[i] == 0) break;
      pf_address = pf_address + (ipPrefetchStrideL2[i] << LOG2_BLOCK_SIZE);
      if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) continue;

      if (PQ.occupancy < PQ.SIZE) {
         assert(prefetch_line(ip, addr, pf_address, FILL_L2, 0));
         did_pref = true;
      }
      else {
         assert(!prefetch_line(ip, addr, pf_address, FILL_L2, 0));
         stopPrefetching = true;
         break;
      }
      num_pref++;
   }

   if ((BASE_PREFETCH_DEGREE_L2 > 1) && !stopPrefetching && (i<BASE_PREFETCH_DEGREE_L2+1)) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//            INJECT  MORE  PREFETCHES  BASED  ON  LAST  TWO  IDENTICAL  DELTAS,  IF  ANY              //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      if ((ipPrefetchStrideL2[i-1] != 0) && (i > 1) && (ipPrefetchStrideL2[i-1] == ipPrefetchStrideL2[i-2])) {
         assert(num_pref < BASE_PREFETCH_DEGREE_L2);

         while (1) {
            pf_address = pf_address + (ipPrefetchStrideL2[i-1] << LOG2_BLOCK_SIZE);
            if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;

            if (PQ.occupancy < PQ.SIZE) {
               assert(prefetch_line(ip, addr, pf_address, FILL_L2, 0));
               num_pref++;
               did_pref = true;
            }
            else {
               assert(!prefetch_line(ip, addr, pf_address, FILL_L2, 0));
               break;
            }
            if (num_pref == BASE_PREFETCH_DEGREE_L2) break;
         }
      }
      else if (fallBackPrefetchStride != 0) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                   INJECT  MORE  PREFETCHES  BASED  ON  IP  STRIDE                                   //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         while (1) {
            pf_address = pf_address + (fallBackPrefetchStride << LOG2_BLOCK_SIZE);
            if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;

            if (PQ.occupancy < PQ.SIZE) {
               assert(prefetch_line(ip, addr, pf_address, FILL_L2, 0));
               num_pref++;
               did_pref = true;
            }
            else {
               assert(!prefetch_line(ip, addr, pf_address, FILL_L2, 0));
               break;
            }
            if (num_pref == BASE_PREFETCH_DEGREE_L2) break;
         }
      }
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INVOKE  NEXT  LINE  PREFETCHER  IF  NO  PREFETCH  YET                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (!did_pref) {
      uint64_t pf_address = (cl_addr + 1) << LOG2_BLOCK_SIZE;
      i = 0;
      while (i < BASE_PREFETCH_DEGREE_L2) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INJECT  PREFETCHES  FOR  APPROPRIATE  DEGREE                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (degreeHitsL2[i]*NL_THRESHOLD_DENOM_L2 > degreeInsertionsL2[i]*NL_THRESHOLD_NUMER_L2) {
            prefetch_line(ip, addr, pf_address, FILL_L2, 0);
         }
         pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
         i++;
      }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    INSERT  POSSIBLE  NEXT  LINE  PREFETCH  CANDIDATES  INTO  NL  BUFFER TO  MEASURE  USEFULNESS     //
////////////////////////////////////////////////////////////////////////////////////////////////////////
      i = 0;
      pf_address = (cl_addr+1) << LOG2_BLOCK_SIZE;
      while (i < BASE_PREFETCH_DEGREE_L2) {
         if ((pf_address >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) break;
         NLBufferL2Insert (pf_address >> LOG2_BLOCK_SIZE, i);
         pf_address = pf_address + (1 << LOG2_BLOCK_SIZE);
         i++;
      }
   }

   return 0;
}

uint32_t CACHE::l2c_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
  return metadata_in;
}

void CACHE::l2c_prefetcher_final_stats()
{
    cout << "CPU " << cpu << " L2C prefetcher final stats" << endl;
}
