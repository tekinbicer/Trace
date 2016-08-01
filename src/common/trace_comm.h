#ifndef DISP_APPS_RECONSTRUCTION_COMMON_TRACE_COMM_H_
#define DISP_APPS_RECONSTRUCTION_COMMON_TRACE_COMM_H_

#include "data_region_a.h"
#include "data_region_2d_bare_base.h"
#include "mpi.h"

namespace trace_comm {
  void GlobalNeighborUpdate(
      ADataRegion<float> &recon_a,
      int num_slices,
      size_t slice_size);

  void SetupSliceCommGroups(
      int slice_id, 
      MPI_Comm &newcomm);

  void UpdateSliceCommGroups(
      ADataRegion<float> &recon_a,
      size_t slice_size,
      MPI_Comm &slice_comm);

  void UpdateReconReplicasCommGroups(
      DataRegion2DBareBase<float> &dr,
      MPI_Comm &group_comm);
}

#endif /// DISP_APPS_RECONSTRUCTION_COMMON_TRACE_COMM_H_
