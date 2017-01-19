#ifndef TRACE_COMMON_TRACEENGINE_H
#define TRACE_COMMON_TRACEENGINE_H

#include "TraceData.h"
#include "TraceRuntimeConfig.h"
#include "data_region_base.h"
#include "disp_engine_reduction.h"

class TraceEngine {
  private:
    TraceData &trace_data;
    TraceRuntimeConfig &config;

    std::unique_ptr<DISPEngineBase<AReconSpace, float>> engine = nullptr;
    AReconSpace *main_recon_space = nullptr;

  public:
    TraceEngine(TraceData &trace_data, TraceRuntimeConfig &conf);

    void IterativeReconstruction(TraceData &trace_data, int iteration);
    void IterativeReconstruction();
};

#endif // TRACE_COMMON_TRACEENGINE_H
