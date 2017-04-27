#include "recon_space.h"
#include "sirt.h"

void SIRTReconSpace::Initialize(int n_grids)
{
  AReconSpace::Initialize(n_grids); 
  leng2= new float[2*num_grids];
}

void SIRTReconSpace::Finalize()
{
  AReconSpace::Finalize();
  delete [] leng2;
}

// Backprojection
void SIRTReconSpace::UpdateRecon(
    ADataRegion<float> &recon,                  // Reconstruction object
    DataRegion2DBareBase<float> &comb_replica) {// Locally combined replica
  size_t rows = comb_replica.rows();
  size_t cols = comb_replica.cols()/2;
  for(size_t i=0; i<rows; ++i){
    auto replica = comb_replica[i];
    for(size_t j=0; j<cols; ++j)
      recon[i*cols + j] +=
        replica[j*2] / replica[j*2+1];
  }
}


void SIRTReconSpace::UpdateReconReplica(
    float simdata,
    float ray,
    int curr_slice,
    int const * const indi,
    float *leng2,
    float *leng, 
    int len)
{
  float upd=0., a2=0.;

  auto &slice_t = reduction_objects()[curr_slice];
  auto slice = &slice_t[0];

  for (int i=0; i<len-1; ++i)
    a2 += leng2[i];

  upd = (ray-simdata) / a2;

  int i=0;
  for (; i<(len-1); ++i) {
#ifdef PREFETCHON
    size_t index2 = indi[i+32]*2;
    __builtin_prefetch(slice+index2,1,0);
#endif
    size_t index = indi[i]*2;
    slice[index] += leng[i]*upd; 
    slice[index+1] += leng[i];
  }
}

void SIRTReconSpace::Reduce(MirroredRegionBareBase<float> &input)
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
          leng, leng2, 
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
          leng2, leng,
          len);
      /*******************************************************/
    }
  }
}
