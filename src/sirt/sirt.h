#ifndef DISP_APPS_RECONSTRUCTION_SIRT_SIRT_H
#define DISP_APPS_RECONSTRUCTION_SIRT_SIRT_H

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

#define SIMD
#ifdef SIMD
#include <immintrin.h>
#endif


class SIRTReconSpace : 
  public AReductionSpaceBase<SIRTReconSpace, float>
{
  private:
#ifdef SIMD //512-bit alignment
    __declspec(align(64)) float *coordx = nullptr;
    __declspec(align(64)) float *coordy = nullptr;
    __declspec(align(64)) float *ax = nullptr;
    __declspec(align(64)) float *ay = nullptr;
    __declspec(align(64)) float *bx = nullptr;
    __declspec(align(64)) float *by = nullptr;
    __declspec(align(64)) float *coorx = nullptr;
    __declspec(align(64)) float *coory = nullptr;
    __declspec(align(64)) float *leng = nullptr;
    __declspec(align(64)) float *leng2 = nullptr;
    __declspec(align(64)) int *indi = nullptr;

    __declspec(align(64)) int num_grids;
#else
    float *coordx = nullptr;
    float *coordy = nullptr;
    float *ax = nullptr;
    float *ay = nullptr;
    float *bx = nullptr;
    float *by = nullptr;
    float *coorx = nullptr;
    float *coory = nullptr;
    float *leng = nullptr;
    float *leng2 = nullptr;
    int *indi = nullptr;

    int num_grids;
#endif

  protected:
    // Forward projection
    float CalculateSimdata(
        float *recon,
        int len,
        int *indi,
        float *leng);

    void UpdateReconReplica(
        float simdata,
        float ray,
        int curr_slice,
        int const * const indi,
        float *leng2,
        float *leng, 
        int len);

  public:
    SIRTReconSpace(int rows, int cols) : 
      AReductionSpaceBase<SIRTReconSpace, float>(rows, cols) {}

    virtual ~SIRTReconSpace(){
      Finalize();
    }

    void Reduce(MirroredRegionBareBase<float> &input);
    // Backward Projection
    void UpdateRecon(
        ADataRegion<float> &recon,                  // Reconstruction object
        DataRegion2DBareBase<float> &comb_replica); // Locally combined replica


    void Initialize(int n_grids);
    virtual void CopyTo(SIRTReconSpace &target){
      target.Initialize(num_grids);
    }
    void Finalize();
};

#endif    // DISP_APPS_RECONSTRUCTION_SIRT_SIRT_H
