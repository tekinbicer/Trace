#ifndef TRACE_SRC_PML_H_
#define TRACE_SRC_PML_H_

#include "recon_space.h"
#include "trace_data.h"
#include "reduction_space_a.h"
#include "data_region_base.h"

class PMLDataRegion final : public DataRegionBase<float, TraceMetadata>{
  private:
    float *F_ = nullptr;
    float *G_ = nullptr;

  public:
    const float kWeightIn[8] {
      0.1464466094, 0.1464466094, 0.1464466094, 0.1464466094,
        0.10355339059, 0.10355339059, 0.10355339059, 0.10355339059 };
    const float kWeightEdges[5] {
      0.226540919667, 0.226540919667, 0.226540919667,
        0.1601886205, 0.1601886205 };
    const float kWeightCorners[3] {
      0.36939806251, 0.36939806251, 0.26120387496 };

    PMLDataRegion(
        float *data, 
        size_t count, 
        TraceMetadata *metadata);

    PMLDataRegion(DataRegionBase<float, TraceMetadata> &region);

    virtual ~PMLDataRegion();

    void SetFG(float val);

    float* F() const;
    float* G() const;
};

class PMLReconSpace final : public AReconSpace
{
  private:
    void UpdateReconReplica(
        float simdata,
        float ray,
        float *recon,
        int curr_slice,
        int const * const indi,
        float *norms, 
        int len);


  protected:
    void PartialBackProjection();


  public:
    PMLReconSpace(int rows, int cols);

    void UpdateRecon(
        TraceData &trace_data,
        DataRegion2DBareBase<float> &comb_replica);  // Locally combined replica

    void CalculateFG(
        ADataRegion<float> &slices_,
        float beta);

    virtual PMLReconSpace* Clone();
};

#endif
