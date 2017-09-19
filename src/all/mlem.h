#ifndef TRACE_SRC_MLEM_H_
#define TRACE_SRC_MLEM_H_

#include "recon_space.h"

class MLEMReconSpace final : public AReconSpace
{
  private:
    void UpdateReconReplica(
        float simdata,
        float ray,
        int curr_slice,
        int const * const indi,
        float *leng, 
        int len);

  public:
    MLEMReconSpace(int rows, int cols) : 
      AReconSpace(rows, cols) {};

    // Backprojection
    void UpdateRecon(
        TraceData &trace_data,
        DataRegion2DBareBase<float> &comb_replica);  // Locally combined replica

    void UpdateReconReplica();

    void PartialBackProjection();

    virtual MLEMReconSpace* Clone();
};

#endif
