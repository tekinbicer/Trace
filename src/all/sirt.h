#ifndef TRACE_SRC_SIRT_H_
#define TRACE_SRC_SIRT_H_

#include "recon_space.h"

class SIRTReconSpace final : public AReconSpace {
  public:
    SIRTReconSpace(int rows, int cols);

    void Finalize();

    /// Backprojection
    void UpdateRecon(
        TraceData &trace_data,
        DataRegion2DBareBase<float> &comb_replica); /// Locally combined replica

    void UpdateReconReplica(
        float /* simdata */,
        float /* ray */,
        int /* curr_slice */,
        int const * const /* indi */,
        float * /* leng2 */,
        float * /* leng */, 
        int /* len */);

    void PartialBackProjection();

    virtual SIRTReconSpace* Clone();
};

#endif
