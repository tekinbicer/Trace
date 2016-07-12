#ifndef BGPMProfiler
#define BGPMProfiler

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>
#include "bgpm.h"

//#undef __INLINE__

#define MAX_COUNTERS 24
#define PUEventSetSize 5
#define numL2Evts 5

int hL2EvtSet;

unsigned int EventsPu[PUEventSetSize] = {
  PEVT_CYCLES, 
  PEVT_INST_ALL,  
  PEVT_AXU_FP_EXCEPT,
  PEVT_INST_QFPU_FMUL,
  PEVT_AXU_DENORM_FLUSH};

unsigned int EventsL2[numL2Evts] = {
  PEVT_L2_HITS,
  PEVT_L2_MISSES,
  PEVT_L2_PREFETCH,
  PEVT_L2_FETCH_LINE,
  PEVT_L2_STORE_LINE};

extern "C" {
//Array to hold overflow count for each event in the event set
int ovfArray[MAX_COUNTERS];

// Calculate total number of overflow for the eventset.
int GetTotalOvfs(unsigned int hEvtSet, int *ovfs);

// Fresh start of Array holding overlow count for each event.
void reset_ofarray();

//This is the overflow handler.
void OvfHandler(int hEvtSet, uint64_t address, uint64_t ovfVector, const ucontext_t *pContext);

void bgpm_profiler_init();

void bgpm_profiler_finalize();

void bgpm_profiler_start();

void bgpm_profiler_stop();

void bgpm_profiler_print();
}

#endif
