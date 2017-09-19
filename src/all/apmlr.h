#ifndef TRACE_SRC_APMLR_H
#define TRACE_SRC_APMLR_H 

#include "recon_space.h"
#include "trace_data.h"
#include "reduction_space_a.h"
#include "data_region_base.h"

class APMLRDataRegion : public DataRegionBase<float, TraceMetadata>{
  private:
    float *F_ = nullptr;
    float *G_ = nullptr;

  public:

    APMLRDataRegion(
        float *data, 
        size_t count, 
        TraceMetadata *metadata);

    APMLRDataRegion(DataRegionBase<float, TraceMetadata> &region);

    virtual ~APMLRDataRegion();

    void SetFG(float val);

    float* F() const;
    float* G() const;

};

class APMLRReconSpace final : public AReconSpace
{
  private:
    void CalculateFGInner(
        float *recon, float *F, float *G,
        float beta, float beta1, float delta, float delta1, float regw,
        int num_slices, int num_grids);
    void CalculateFGTop(
        float * recon, float *F, float *G,
        float beta, float beta1, float delta, float delta1, float regw, 
        int num_grids);
    void CalculateFGBottom(
        float * recon, float *F, float *G,
        float beta, float beta1, float delta, float delta1, float regw, 
        int num_grids);

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
    APMLRReconSpace(int rows, int cols);

    void UpdateRecon(
        TraceData &trace_data, // Reconstruction object
        DataRegion2DBareBase<float> &comb_replica); // Locally combined replica

    void CalculateFG(
        ADataRegion<float> &slices_,
        float beta, float beta1,
        float delta, float delta1,
        float regw);

    virtual APMLRReconSpace* Clone();
};

#endif    // TRACE_SRC_APMLR_H
