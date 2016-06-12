#ifdef BGPM_PROFILER

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>
#include "bgpm/include/bgpm.h"

class BGPMProfiler
{
  private: 

    const int MAX_COUNTERS 24;
    const int PUEventSetSize = 5;
    const int numL2Evts = 5;

    int hL2EvtSet;

    unsigned EventsPu[EventSetSize] = {
      PEVT_CYCLES, 
      PEVT_INST_ALL,  
      PEVT_AXU_FP_EXCEPT,
      PEVT_INST_QFPU_FMUL,
      PEVT_AXU_DENORM_FLUSH};
    unsigned EventsL2[numL2Evts] = {
      PEVT_L2_HITS,
      PEVT_L2_MISSES,
      PEVT_L2_PREFETCH,
      PEVT_L2_FETCH_LINE,
      PEVT_L2_STORE_LINE};

    //Array to hold overflow count for each event in the event set
    int ovfArray[MAX_COUNTERS];

    //This is the overflow handler.
    void OvfHandler(int hEvtSet, uint64_t address, uint64_t ovfVector, const ucontext_t *pContext)
    {
      unsigned ovfIdxs[BGPM_MAX_OVERFLOW_EVENTS];
      unsigned len = BGPM_MAX_OVERFLOW_EVENTS;
      //Within an overflow event handler, convert the passed opaque ovfVector into 
      //a list of indicies into the event set to the events which overflowed
      //parameters: 
      //  hEvtSet    =>  handle to event set
      //  ovfVector  =>  eventset unique opaque mask needed by Bgpm_GetOverflowEventIndices() to identify which events have overflowed. 
      //  *pIndicies =>  user allocated array to receive output indicies
      //  *pLen      =>  Input array length / output number of filled indicies

      Bgpm_GetOverflowEventIndices(hEvtSet, ovfVector, ovfIdxs, &len);
      for (unsigned i=0; i<len; i++) {
        unsigned idx = ovfIdxs[i];
        ovfArray[idx]++;
      }

    }

    // Calculate total number of overflow for the eventset.
    int GetTotalOvfs(unsigned hEvtSet, int ovfs[MAX_COUNTERS])
    {
      int numEvts = Bgpm_NumEvents(hEvtSet);
      int total = 0;
      for (int idx=0; idx<numEvts; idx++) {
        total += ovfs[idx];
      }
      return total;
    }

    // Fresh start of Array holding overlow count for each event.
    void reset_ofarray()
    {
      for(int i=0;i<MAX_COUNTERS;i++){
        ovfArray[i]=0;
      }
    }


  public:
    void init() 
    {
      reset_ofarray();
      Bgpm_Init(BGPM_MODE_SWDISTRIB);

      hL2EvtSet = Bgpm_CreateEventSet();
      Bgpm_AddEventList(hL2EvtSet, evtL2List, sizeof(evtL2List)/sizeof(unsigned));

      Bgpm_SetOverflow(hL2EvtSet, Bgpm_GetEventIndex(hL2EvtSet, PEVT_L2_MISSES, 0), 1000);
      Bgpm_SetOverflow(hL2EvtSet, Bgpm_GetEventIndex(hL2EvtSet, PEVT_L2_HITS, 0), 5000);
      Bgpm_SetOverflow(hL2EvtSet, Bgpm_GetEventIndex(hL2EvtSet, PEVT_L2_STORE_LINE, 0), 1000);
      //If an event has Overflow set, then the handler will be scheduled to be called as a child of an BGPM signal handler.
      Bgpm_SetOverflowHandler(hL2EvtSet, OvfHandler); 

      // Apply eventset
      Bgpm_Apply(hL2EvtSet);
    }

    void finalize() {
      Bgpm_Disable();
    }

    void start()
    {
      Bgpm_Start(hL2EvtSet);
    }

    void stop()
    {
      Bgpm_Stop(hL2EvtSet);
    }

    void print()
    {
      int i;

      int numEvts = Bgpm_NumEvents(hEvtSet);
      uint64_t cnt;
      printf("  %18s %20s %10s %10s\n\n","EVENT COUNT","EVENT LABEL","PERIOD","OVERFLOWS");
      for (i=0; i<numEvts; i++) {
        Read counter for given event, and current thread. 
          //Counts will be availible in "cnt" array.
          Bgpm_ReadEvent(hEvtSet, i, &cnt);
        uint64_t period;
        Bgpm_OverflowHandler_t handler;
        Bgpm_GetOverflow(hEvtSet, i, &period, &handler);
        if(period){
          printf("  0x%018lx %30s \t%10ld %5d\n", cnt, Bgpm_GetEventLabel(hEvtSet, i),period,ovfArray[i]);
        }
        else{
          printf("  0x%018lx %30s\n", cnt, Bgpm_GetEventLabel(hEvtSet, i));
        }
      }
    }
};

#endif // BGPM_PROFILER
