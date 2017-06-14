#include "mlem.h"

void MLEMReconSpace::UpdateRecon(
    ADataRegion<float> &recon,                  // Reconstruction object
    DataRegion2DBareBase<float> &comb_replica)  // Locally combined replica
{
  size_t rows = comb_replica.rows();
  size_t cols = comb_replica.cols()/2;
  for(size_t i=0; i<rows; ++i){
    auto replica = comb_replica[i];
    for(size_t j=0; j<cols; ++j){
      recon[i*cols + j] *=
        replica[j*2] / replica[j*2+1];
    }
  }
//  size_t rows = comb_replica.rows()/2;
//  size_t cols = comb_replica.cols();
//  for(size_t i=0; i<rows; ++i){
//    auto weights = comb_replica[i*2];
//    auto lengths = comb_replica[i*2+1];
//    for(size_t j=0; j<cols; ++j){
//      recon[i*cols + j] *= weights[j] / lengths[j];
//    }
//  }
}

void MLEMReconSpace::UpdateReconReplica(
    float simdata,
    float ray,
    int curr_slice,
    int const * const restrict indi,
    float *restrict leng, 
    int len)
{

  auto &slice_t = reduction_objects()[curr_slice];
  float * restrict slice = &slice_t[0];
  size_t off = reduction_objects().cols()/2;
  float * restrict slice2 = slice+off;

//  float * restrict slice = &(reduction_objects()[curr_slice][0]);
//  float * restrict slice2 = &(reduction_objects()[curr_slice+1][0]);

  UpdateReconReplicaV(
    simdata,
    ray,
    indi,
    leng, 
    len,
    slice,
    slice2);
}

void MLEMReconSpace::UpdateReconReplicaV(
    float simdata,
    float ray,
    int const * const restrict indi,
    float *restrict leng, 
    int len,
    float *restrict slice,
    float *restrict slice2)
{
  float upd;
  upd = (ray/simdata);

  #if defined(__AVX512F__) && defined(T_KNL_OPTIMIZED)
  #pragma message("KNL code is being used")
  __assume_aligned(indi, 64);
  __assume_aligned(slice, 64);
  __assume_aligned(slice2, 64);
  __assume_aligned(leng, 64);

  //#pragma prefetch slice:1:3  // prefetch to L2 cache 3 iterations ahead
  #pragma ivdep
  #endif
  for (int i=0; i <len-1; ++i) {
    size_t index = indi[i]*2;
    slice[index] += leng[i]*upd;
    slice[index+1] += leng[i];
  }

  //#pragma vector always assert
//  #pragma ivdep
//  for (int i=0; i <len-1; ++i) {
//    size_t index = indi[i];
//    slice[index] += leng[i]*upd;
//    slice2[index] += leng[i];
//  }
}

MLEMReconSpace* MLEMReconSpace::Clone()
{
  auto &red_objs = reduction_objects();

  MLEMReconSpace *cloned_obj = new MLEMReconSpace(red_objs.rows(), red_objs.cols());
  (*cloned_obj).reduction_objects() = red_objs;

  static_cast<MLEMReconSpace*>(this)->CopyTo(*cloned_obj);

  return cloned_obj;
}
