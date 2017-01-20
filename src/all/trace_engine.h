#ifndef TRACE_COMMON_TRACEENGINE_H
#define TRACE_COMMON_TRACEENGINE_H

#include "trace_engine_data.h"
#include "trace_data.h"
#include "trace_runtime_config.h"
#include "data_region_base.h"
#include "disp_engine_reduction.h"
#include "recon_space.h"
#include "disp_comm_base.h"

class TraceEngine {
  private:
    TraceData &trace_data;
    DISPCommBase<float> &comm;
    TraceRuntimeConfig &config;

    std::unique_ptr<DISPEngineBase<AReconSpace, float>> engine = nullptr;
    AReconSpace *main_recon_space = nullptr;

  public:
    TraceEngine(TraceData &trace_data, DISPCommBase<float> &comm, TraceRuntimeConfig &conf);

    void IterativeReconstruction(TraceData &trace_data, int iteration);
    void IterativeReconstruction();
};

#endif // TRACE_COMMON_TRACEENGINE_H