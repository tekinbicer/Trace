#include "mlem.h"

void MLEMReconSpace::UpdateRecon(
    TraceData &trace_data,
    DataRegion2DBareBase<float> &comb_replica)  // Locally combined replica
{
  ADataRegion<float> &recon = trace_data.metadata().recon();
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
    float *norms, 
    int len)
{
  float upd;

  auto &slice_t = reduction_objects()[curr_slice];
  auto slice = &slice_t[0];

  upd = (ray/simdata);

  for (int i=0; i <len-1; ++i) {
    size_t index = indi[i]*2;
    slice[index] += norms[i]*upd;
    slice[index+1] += norms[i];
  }
}

void MLEMReconSpace::PartialBackProjection()
{
  UpdateReconReplica(
      reconparams.simdata,
      (*reconparams.rays)[reconparams.curr_col],
      reconparams.curr_slice,
      reconparams.indi,
      reconparams.norms,
      reconparams.len);
}


MLEMReconSpace* MLEMReconSpace::Clone()
{
  auto &red_objs = reduction_objects();

  MLEMReconSpace *cloned_obj = new MLEMReconSpace(red_objs.rows(), red_objs.cols());
  (*cloned_obj).reduction_objects() = red_objs;

  static_cast<MLEMReconSpace*>(this)->CopyTo(*cloned_obj);

  return cloned_obj;
}
