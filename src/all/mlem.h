#ifndef TRACE_SRC_MLEM_H_
#define TRACE_SRC_MLEM_H_

#include "recon_space.h"

class MLEMReconSpace final : public AReconSpace
{
  public:
    MLEMReconSpace(int rows, int cols) : 
      AReconSpace(rows, cols) {};

    // Backprojection
    void UpdateRecon(
        ADataRegion<float> &recon,                  // Reconstruction object
        DataRegion2DBareBase<float> &comb_replica);  // Locally combined replica

    void UpdateReconReplica(
        float simdata,
        float ray,
        int curr_slice,
        int const * const indi,
        float *leng, 
        int len);
};

#endif
