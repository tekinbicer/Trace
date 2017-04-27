#include "trace_runtime_config.h"
#include "trace_data.h"
#include "trace_ph5io.h"
#include "data_region_base.h"
#include "disp_comm_mpi.h"
#include "disp_engine_reduction.h"
#include "trace_engine.h"
#include "sirt.h"
#include "mlem.h"
#include "pml.h"

int main(int argc, char **argv)
{
  /* Parse command line parameters */
  TraceRuntimeConfig config(argc, argv);
  config.Print();

  /* Initiate middleware's communication layer */
  DISPCommMPI<float> comm(&argc, &argv);

  TracePH5IO trace_io(comm, config);

  TraceData trace_data = trace_io.Read();

  TraceEngine trace_engine(trace_data, comm, config);

  trace_engine.IterativeReconstruction();

  trace_io.Write(trace_data);
}

