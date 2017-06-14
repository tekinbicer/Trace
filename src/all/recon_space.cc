#include "recon_space.h"

// Forward projection
float AReconSpace::CalculateSimdata(
    float *recon,
    int len,
    int *indi,
    float *leng) 
{
  float simdata = 0.;
  for(int i=0; i<len-1; ++i){
    simdata += recon[indi[i]]*leng[i];
  }
  return simdata;
}

/// SIRT
void AReconSpace::UpdateReconReplica(
  float /* simdata */,
  float /* ray */,
  int /* curr_slice */,
  int const * const /* indi */,
  float * /* leng2 */,
  float * /* leng */, 
  int /* len */)
{}

/// MLEM
void AReconSpace::UpdateReconReplica(
  float /* simdata */,
  float /* ray */,
  int /* curr_slice */,
  int const * const /* indi */,
  float * /* leng */, 
  int /* len */)
{}

/// PML
void AReconSpace::UpdateReconReplica(
  float /* simdata */,
  float /* ray */,
  float * /* recon */,
  int /* curr_slice */,
  int const * const /* indi */,
  float * /* leng */, 
  int /* len */)
{}

void AReconSpace::UpdateRecon(
    ADataRegion<float> & /* recon */,
    DataRegion2DBareBase<float> & /*comb_replica */)
{}


AReconSpace::AReconSpace(int rows, int cols) : 
  AReductionSpaceBase<AReconSpace, float>(rows, cols) 
{}

AReconSpace::~AReconSpace()
{
  #if defined(__AVX512F__) && defined(T_KNL_OPTIMIZED)
  _mm_free(coordx);
  _mm_free(coordy);
  _mm_free(ax);
  _mm_free(ay);
  _mm_free(bx);
  _mm_free(by);
  _mm_free(coorx);
  _mm_free(coory);
  _mm_free(leng);
  _mm_free(indi);
  #else
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
  #endif
}

void AReconSpace::Reduce(MirroredRegionBareBase<float> &input)
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

void AReconSpace::Initialize(int n_grids)
{
  num_grids = n_grids; 

  #if defined(__AVX512F__) && defined(T_KNL_OPTIMIZED)
  coordx = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  coordy = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);  
  ax = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  ay = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  bx = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  by = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  coorx = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
  coory = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
  leng = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
  indi = (int *)_mm_malloc((2*num_grids) * sizeof(int), 64);
  #else
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
  #endif
}

void AReconSpace::CopyTo(AReconSpace &target)
{
  target.Initialize(num_grids);
}
