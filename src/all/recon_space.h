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
#include "reduction_space_a.h"
#include "data_region_base.h"
#include "trace_engine_data.h"

constexpr float kPI = 3.14159265358979f;

/// Note: this is a mutable struct.
struct LocalReconParams {

  void Initialize(int n_grids){
    num_grids = n_grids; 

    coordx = new float[num_grids+1]; 
    coordy = new float[num_grids+1];
    ax = new float[num_grids+1];
    ay = new float[num_grids+1];
    bx = new float[num_grids+1];
    by = new float[num_grids+1];
    coorx = new float[2*num_grids];
    coory = new float[2*num_grids];
    leng2 = new float[2*num_grids];
    norms = new float[2*num_grids];
    indi = new int[2*num_grids];
  }

  ~LocalReconParams(){
    delete [] coordx;
    delete [] coordy;
    delete [] ax;
    delete [] ay;
    delete [] bx;
    delete [] by;
    delete [] coorx;
    delete [] coory;
    delete [] leng2;
    delete [] norms;
    delete [] indi;
  }

  /// Calculate indices and distances params
  float *coordx = nullptr;
  float *coordy = nullptr;
  float *ax = nullptr;
  float *ay = nullptr;
  float *bx = nullptr;
  float *by = nullptr;
  float *coorx = nullptr;
  float *coory = nullptr;
  float *leng2 = nullptr;
  float *norms = nullptr;
  int *indi = nullptr;
  int len;


  /// Input data params
  MirroredRegionBase<float, TraceMetadata> *rays;
  TraceMetadata *metadata;

  const float *gridx;
  const float *gridy;
  const float *theta;
  float mov;

  int num_cols;
  int num_grids;

  int curr_proj;
  int count_projs;


  float theta_q;
  int quadrant;
  float sinq;
  float cosq;

  int curr_slice;
  int curr_slice_offset;
  int curr_col;

  float *recon;
  float simdata;

  void InitInputParams(MirroredRegionBareBase<float> &input){
    rays = static_cast<MirroredRegionBase<float, TraceMetadata>*>(&input);
    metadata = &(rays->metadata());
    gridx = metadata->gridx();
    gridy = metadata->gridy();
    theta = metadata->theta();
    mov = metadata->mov();
    num_cols = metadata->num_cols();
    num_grids = metadata->num_cols();
    curr_proj = metadata->RayProjection(rays->index());
    count_projs = metadata->RayProjection(rays->index()+rays->count()-1) - curr_proj;
  }

  void InitLoopParams(int proj) {
    theta_q = theta[proj];
    quadrant = ((theta_q >= 0 && theta_q < kPI/2) ||
        (theta_q >= kPI && theta_q < 3*kPI/2)) ? 1 : 0;
    sinq = sinf(theta_q);
    cosq = cosf(theta_q);
    curr_slice = metadata->RaySlice(rays->index());
    curr_slice_offset = curr_slice*pow(num_grids,2);
    recon = &(metadata->recon()[0])+curr_slice_offset;
  }
};

class AReconSpace : 
  public AReductionSpaceBase<AReconSpace, float>
{
  protected:
    struct LocalReconParams reconparams;

    /// Local Loop Parameters

    /// 
    /// Computes visited indices and ray lengths
    virtual void ComputeIndexDistance();

    /// Computes forward projection, i.e. simulated data
    virtual void ForwardProjection();

    /// Computes backprojection and performs partial updates
    /// on replicated reconstruction objects
    virtual void PartialBackProjection() = 0;



    /// Utility functions
    
    /// Finds coordinates
    virtual void CalculateCoordinates(
      int num_grid,
      float xi, float yi, float sinq, float cosq,
      const float *gridx, const float *gridy,
      float *coordx, float *coordy); /// Outputs

    /// Merge the (coordx, gridy) and (gridx, coordy)
    virtual void MergeTrimCoordinates(
      int num_grid,
      float *coordx, float *coordy,
      const float *gridx, const float *gridy,
      int *alen, int *blen,  /// Outputs
      float *ax, float *ay,  /// Outputs
      float *bx, float *by); /// Outputs

    /// Sort the array of intersection points (ax, ay). The new sorted 
    /// intersection points are stored in (coorx, coory). 
    /// If quadrant=1 then a_ind = i; if 0 then a_ind = (alen-1-i)
    virtual void SortIntersectionPoints(
      int ind_cond,
      int alen, int blen,
      float *ax, float *ay,
      float *bx, float *by,
      float *coorx, float *coory); /// Outputs

    /// Calculate the distance squares (leng) between the intersection points 
    /// (coorx, coory). Find the indices of the pixels on the reconstruction
    /// grid (ind_recon).
    virtual void CalculateDistanceLengths(
      int len, int num_grids,
      float *coorx, float *coory,
      float *leng2, int *indi); /// Outputs

    /// Default function for calculating forward projection simulated data. 
    virtual float CalculateSimdata(
        float *recon,
        int len,
        int *indi,
        float *norms);


  public:
    AReconSpace(int rows, int cols);

    virtual ~AReconSpace();

    virtual void Reduce(MirroredRegionBareBase<float> &input);

    virtual void UpdateRecon(
        TraceData &trace_data,
        DataRegion2DBareBase<float> &comb_replica) = 0;   // Locally combined replica

    virtual void Initialize(int n_grids);

    virtual void CopyTo(AReconSpace &target);

    virtual AReconSpace* Clone()=0;
};

#endif    /// TRACE_RECON_SPACE_H

