#include "sirt.h"

#ifdef SIMD

const int MY_SIMD_WIDTH = 16;

/// Forward Projection SIMD
float SIRTReconSpace::CalculateSimdata(
    float *recon,
    int len,
    int *indi,
    float *leng)
{
  float simdata = 0.;
  int i = 0;
  
  //SIMD processing, note: change computation order
  __m512 vec_simdata = _mm512_setzero_ps();  
  for(; i < len - 1 - MY_SIMD_WIDTH; i += MY_SIMD_WIDTH){
    //simdata += recon[indi[i]]*leng[i];
    __m512 vec_leng = _mm512_load_ps((__m512*)(leng + i));
    __m512i vec_indi = _mm512_load_epi32((__m512i*)(indi + i));
    __m512 vec_recon = _mm512_i32gather_ps(vec_indi, recon, 4);
    vec_simdata = _mm512_add_ps(vec_simdata, _mm512_mul_ps(vec_recon, vec_leng));
  }

  simdata = _mm512_reduce_add_ps(vec_simdata); 

  //Process the rest
  for(; i < len-1; ++i){
    simdata += recon[indi[i]]*leng[i];
  }

  return simdata;
}

//SIMD
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

  int i = 0;
  //SIMD processing, note: change computation order
  __m512 vec_a2 = _mm512_setzero_ps();  
  for(; i < len - 1 - MY_SIMD_WIDTH; i += MY_SIMD_WIDTH){
    //a2 += leng2[i];
    __m512 vec_leng2 = _mm512_load_ps((__m512*)(leng2 + i));
    vec_a2 = _mm512_add_ps(vec_a2, vec_leng2);
  }

  a2 = _mm512_reduce_add_ps(vec_a2); 

  //Process the rest
  for (; i<len-1; ++i)
    a2 += leng2[i];

  upd = (ray-simdata) / a2;

  i=0;

  //TODO: OPT
  //SIMD processing
  for(; i < len - 1 - MY_SIMD_WIDTH; i += MY_SIMD_WIDTH){
    //size_t index = indi[i]*2;
    //slice[index] += leng[i]*upd; 
    //slice[index+1] += leng[i];
    __m512i vec_indi = _mm512_load_epi32((__m512i*)(indi + i));
    __m512 vec_leng = _mm512_load_ps((__m512*)(leng + i));
    //TODO: scatter with conflict detection
    vec_indi = _mm512_add_epi32(vec_indi, vec_indi);
    __m512 vec_slice_i = _mm512_i32gather_ps(vec_indi, slice, 4);
    __m512 vec_slice_i_1 = _mm512_i32gather_ps(_mm512_add_epi32(vec_indi, _mm512_set1_epi32(1)), slice, 4);

    vec_slice_i = _mm512_add_ps(vec_slice_i, _mm512_mul_ps(vec_leng, _mm512_set1_ps(upd)));
    vec_slice_i_1 = _mm512_add_ps(vec_slice_i_1, vec_leng);
    _mm512_i32scatter_ps(slice, vec_indi, vec_slice_i, 4);
    _mm512_i32scatter_ps(slice, _mm512_add_epi32(vec_indi, _mm512_set1_epi32(1)), vec_slice_i_1, 4);
  }

  //Process the rest
  //#pragma vector aligned
  for (; i<(len-1); ++i) {
    size_t index = indi[i]*2;
    slice[index] += leng[i]*upd; 
    slice[index+1] += leng[i];
  }
}



#else //Scalar

/// Forward Projection
float SIRTReconSpace::CalculateSimdata(
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


#endif


void SIRTReconSpace::UpdateRecon(
    ADataRegion<float> &recon,                  // Reconstruction object
    DataRegion2DBareBase<float> &comb_replica)  // Locally combined replica
{
  size_t rows = comb_replica.rows();
  size_t cols = comb_replica.cols()/2;
  for(size_t i=0; i<rows; ++i){
    auto replica = comb_replica[i];
    for(size_t j=0; j<cols; ++j)
      recon[i*cols + j] +=
        replica[j*2] / replica[j*2+1];
  }
}


void SIRTReconSpace::Initialize(int n_grids){
  num_grids = n_grids; 

#ifdef SIMD//512-bit alignment
  coordx = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  coordy = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);  
  ax = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  ay = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  bx = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  by = (float *)_mm_malloc((num_grids+1) * sizeof(float), 64);
  coorx = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
  coory = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
  leng = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
  leng2 = (float *)_mm_malloc((2*num_grids) * sizeof(float), 64);
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
  leng2 = new float[2*num_grids];
  indi = new int[2*num_grids];
#endif
}

void SIRTReconSpace::Finalize(){
#ifdef SIMD
  _mm_free(coordx);
  _mm_free(coordy);
  _mm_free(ax);
  _mm_free(ay);
  _mm_free(bx);
  _mm_free(by);
  _mm_free(coorx);
  _mm_free(coory);
  _mm_free(leng);
  _mm_free(leng2);
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
  delete [] leng2;
  delete [] indi;
#endif
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
