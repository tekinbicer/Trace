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
        float *leng) 
    {
      float simdata = 0.;
      for(int i=0; i<len-1; ++i){
        #ifdef PREFETCHON
        size_t index = indi[i+32];
        __builtin_prefetch(&(recon[index]),1,0);
        #endif
        simdata += recon[indi[i]]*leng[i];
      }
      return simdata;
    }

    /// SIRT
    virtual void UpdateReconReplica(
        float simdata,
        float ray,
        int curr_slice,
        int const * const indi,
        float *leng2,
        float *leng, 
        int len) {};

    /// MLEM
    virtual void UpdateReconReplica(
        float simdata,
        float ray,
        int curr_slice,
        int const * const indi,
        float *leng, 
        int len) {};

    /// PML
    virtual void UpdateReconReplica(
        float simdata,
        float ray,
        float *recon,
        int curr_slice,
        int const * const indi,
        float *leng, 
        int len) {};

  public:
    AReconSpace(int rows, int cols) : 
      AReductionSpaceBase<AReconSpace, float>(rows, cols) {}

    virtual ~AReconSpace(){
      Finalize();
    };

    virtual void Reduce(MirroredRegionBareBase<float> &input)
    {
      auto &rays = *(static_cast<MirroredRegionBase<float, TraceMetadata>*>(&input));
      auto &metadata = rays.metadata();

      const float *theta = metadata.theta();
      const float *gridx = metadata.gridx();
      const float *gridy = metadata.gridy();
      float mov = metadata.mov();

      /* In-memory values */
      int num_cols = metadata.num_cols();
      int num_grids = metadata.num_cols();

      int curr_proj = metadata.RayProjection(rays.index());
      int count_projs = 
        metadata.RayProjection(rays.index()+rays.count()-1) - curr_proj;

      /* Reconstruction start */
      for (int proj = curr_proj; proj<=(curr_proj+count_projs); ++proj) {
        float theta_q = theta[proj];
        int quadrant = trace_utils::CalculateQuadrant(theta_q);
        float sinq = sinf(theta_q);
        float cosq = cosf(theta_q);

        int curr_slice = metadata.RaySlice(rays.index());
        int curr_slice_offset = curr_slice*num_grids*num_grids;
        float *recon = (&(metadata.recon()[0])+curr_slice_offset);

        for (int curr_col=0; curr_col<num_cols; ++curr_col) {
          /// Calculate coordinates
          float xi = -1e6;
          float yi = (1-num_cols)/2. + curr_col+mov;
          trace_utils::CalculateCoordinates(
              num_grids, 
              xi, yi, sinq, cosq, 
              gridx, gridy, 
              coordx, coordy);  /// Outputs coordx and coordy

          /// Merge the (coordx, gridy) and (gridx, coordy)
          /// Output alen and after
          int alen, blen;
          trace_utils::MergeTrimCoordinates(
              num_grids, 
              coordx, coordy, 
              gridx, gridy, 
              &alen, &blen, 
              ax, ay, bx, by);

          /// Sort the array of intersection points (ax, ay)
          /// The new sorted intersection points are
          /// stored in (coorx, coory).
          /// if quadrant=1 then a_ind = i; if 0 then a_ind = (alen-1-i)
          trace_utils::SortIntersectionPoints(
              quadrant, 
              alen, blen, 
              ax, ay, bx, by, 
              coorx, coory);

          /// Calculate the distances (leng) between the
          /// intersection points (coorx, coory). Find
          /// the indices of the pixels on the
          /// reconstruction grid (ind_recon).
          int len = alen + blen;
          trace_utils::CalculateDistanceLengths(
              len, 
              num_grids, 
              coorx, coory, 
              leng,
              indi);

          /*******************************************************/
          /* Below is for updating the reconstruction grid and
           * is algorithm specific part.
           */
          /// Forward projection
          float simdata = CalculateSimdata(recon, len, indi, leng);

          /// Update recon 
          UpdateReconReplica(
              simdata, 
              rays[curr_col], 
              curr_slice, 
              indi, 
              leng,
              len);
          /*******************************************************/
        }
      }
    }

    virtual void UpdateRecon(
        ADataRegion<float> &recon,                   // Reconstruction object
        DataRegion2DBareBase<float> &comb_replica){} // Locally combined replica

    virtual void Initialize(int n_grids){
      num_grids = n_grids; 

      coordx = new float[num_grids+1]; 
      coordy = new float[num_grids+1];
      ax = new float[num_grids+1];
      ay = new float[num_grids+1];
      bx = new float[num_grids+1];
      by = new float[num_grids+1];
      coorx = new float[2*num_grids];
      coory = new float[2*num_grids];
      leng = new float[2*num_grids];
      indi = new int[2*num_grids];
    }

    virtual void CopyTo(AReconSpace &target){
      target.Initialize(num_grids);
    }

    virtual void Finalize(){
      delete [] coordx;
      delete [] coordy;
      delete [] ax;
      delete [] ay;
      delete [] bx;
      delete [] by;
      delete [] coorx;
      delete [] coory;
      delete [] leng;
      delete [] indi;
    }
};

#endif    /// TRACE_RECON_SPACE_H

