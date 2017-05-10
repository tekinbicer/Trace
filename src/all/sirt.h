#ifndef TRACE_SRC_SIRT_H_
#define TRACE_SRC_SIRT_H_

#include "recon_space.h"

class SIRTReconSpace final : public AReconSpace {
  private:

  protected:
    float *  leng2 = nullptr;

  public:
    SIRTReconSpace(int rows, int cols);

    void Initialize(int n_grids);
    void Finalize();

    /// Backprojection
    void UpdateRecon(
        ADataRegion<float> &recon,                  /// Reconstruction object
        DataRegion2DBareBase<float> &comb_replica); /// Locally combined replica

    void UpdateReconReplica(
        float /* simdata */,
        float /* ray */,
        int /* curr_slice */,
        int const * const /* indi */,
        float * /* leng2 */,
        float * /* leng */, 
        int /* len */);

    void Reduce(MirroredRegionBareBase<float> &input);

    virtual SIRTReconSpace* Clone();

    virtual ~SIRTReconSpace();
};

#endif
