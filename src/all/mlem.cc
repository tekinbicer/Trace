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
}

void MLEMReconSpace::UpdateReconReplica(
    float simdata,
    float ray,
    int curr_slice,
    int const * const indi,
    float *leng, 
    int len)
{
  float upd;

  auto &slice_t = reduction_objects()[curr_slice];
  auto slice = &slice_t[0];

  upd = (ray/simdata);

  for (int i=0; i <len-1; ++i) {
#ifdef PREFETCHON
    size_t index2 = indi[i+32]*2;
    __builtin_prefetch(slice+index2,1,0);
#endif
    size_t index = indi[i]*2;
    slice[index] += leng[i]*upd;
    slice[index+1] += leng[i];
  }
}


