#include "mpi.h"
#include "tclap/CmdLine.h"
#include "TraceRuntimeConfig.h"
#include "TraceData.h"
#include "TracePH5IO.h"
#include "data_region_base.h"
#include "disp_comm_mpi.h"
#include "disp_engine_reduction.h"
#include "sirt.h"
#include "mlem.h"
#include "pml.h"

int main(int argc, char **argv)
{
  /* Parse command line parameters */
  TraceRuntimeConfig config(argc, argv);

  /* Initiate middleware's communication layer */
  std::shared_ptr<DISPCommBase<float>> comm = 
    std::make_shared<DISPCommMPI<float>>(&argc, &argv);

  TracePH5IO trace_io(config);

  TraceData trace_data = trace_io.Read();

  TraceEngine trace_engine(trace_data, config);

  trace_engine.IterativeReconstruction();

  trace_io.Write(trace_engine.image(), trace_data);
}

