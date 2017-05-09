#ifndef TRACE_RECON_SPACE_H
#define TRACE_RECON_SPACE_H 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <chrono>
#include "hdf5.h"
#include "string.h"
#include "trace_data.h"
#include "trace_utils.h"
#include "reduction_space_a.h"
#include "data_region_base.h"

class AReconSpace : 
  public AReductionSpaceBase<AReconSpace, float>
{
  protected:
    float *coordx = nullptr;
    float *coordy = nullptr;
    float *ax = nullptr;
    float *ay = nullptr;
    float *bx = nullptr;
    float *by = nullptr;
    float *coorx = nullptr;
    float *coory = nullptr;
    float *leng = nullptr;
    int *indi = nullptr;

    int num_grids;

    // Forward projection
    virtual float CalculateSimdata(
        float *recon,
        int len,
        int *indi,
        float *leng);

    /// SIRT
    virtual void UpdateReconReplica(
        float /* simdata */,
        float /* ray */,
        int /* curr_slice */,
        int const * const /* indi */,
        float * /* leng2 */,
        float * /* leng */, 
        int /* len */);

    /// MLEM
    virtual void UpdateReconReplica(
        float /* simdata */,
        float /* ray */,
        int /* curr_slice */,
        int const * const /* indi */,
        float * /* leng */, 
        int /* len */);

    /// PML
    virtual void UpdateReconReplica(
        float /* simdata */,
        float /* ray */,
        float * /* recon */,
        int /* curr_slice */,
        int const * const /* indi */,
        float * /* leng */, 
        int /* len */);

  public:
    AReconSpace(int rows, int cols);

    virtual ~AReconSpace();

    virtual void Reduce(MirroredRegionBareBase<float> &input);

    virtual void UpdateRecon(
        ADataRegion<float> & /* recon */,                   // Reconstruction object
        DataRegion2DBareBase<float> & /*comb_replica */);   // Locally combined replica

    virtual void Initialize(int n_grids);

    virtual void CopyTo(AReconSpace &target);

    virtual AReconSpace* Clone()=0;

    virtual void Finalize();
};

#endif    /// TRACE_RECON_SPACE_H

