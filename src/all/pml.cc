#include "pml.h"

/** PMLDataRegion */

PMLDataRegion::PMLDataRegion(
    float *data, 
    size_t count, 
    TraceMetadata *metadata): 
  DataRegionBase<float, TraceMetadata>(data, count, metadata) 
{
  size_t count_i = 
    metadata->num_slices()*metadata->num_grids()*metadata->num_grids();
  F_ = new float[count_i];
  G_ = new float[count_i];
  SetFG(0.);
}

PMLDataRegion::PMLDataRegion(DataRegionBase<float, TraceMetadata> &region) :
  PMLDataRegion(&region[0], region.count(), &region.metadata())
{ }

PMLDataRegion::~PMLDataRegion() {
  delete [] F_;
  delete [] G_;
}

void PMLDataRegion::SetFG(float val)
{
  size_t count = 
    metadata().num_slices()*
      metadata().num_grids()*
      metadata().num_grids();
  for(size_t i=0; i<count; ++i){
    F_[i]=val;
    G_[i]=val;
  }
}

float* PMLDataRegion::F() const { return F_; }
float* PMLDataRegion::G() const { return G_; }



/** PMLReconSpace */

PMLReconSpace::PMLReconSpace(int rows, int cols) :
  AReconSpace(rows, cols) {}

PMLReconSpace* PMLReconSpace::Clone()
{
  auto &red_objs = reduction_objects();

  PMLReconSpace *cloned_obj = new PMLReconSpace(red_objs.rows(), red_objs.cols());
  (*cloned_obj).reduction_objects() = red_objs;

  static_cast<PMLReconSpace*>(this)->CopyTo(*cloned_obj);

  return cloned_obj;
}

void PMLReconSpace::CalculateFG(
    ADataRegion<float> &slices_,
    float beta)
{
  PMLDataRegion &slices = dynamic_cast<PMLDataRegion&>(slices_); 
  int num_slices = slices.metadata().num_slices();
  int num_grids = slices.metadata().num_grids();

  ADataRegion<float> &recon = slices.metadata().recon();
  float *F = slices.F();
  float *G = slices.G();

  int k, n, m, q, i;
  int ind0, indg[8];

  /// (inner region)
  const float *wg = slices.kWeightIn;
  for (k = 0; k < num_slices; k++) {
    for (n = 1; n < num_grids - 1; n++) {
      for (m = 1; m < num_grids - 1; m++) {
        ind0 = m + n * num_grids + k * num_grids * num_grids;

        indg[0] = ind0 + 1;
        indg[1] = ind0 - 1;
        indg[2] = ind0 + num_grids;
        indg[3] = ind0 - num_grids;
        indg[4] = ind0 + num_grids + 1;
        indg[5] = ind0 + num_grids - 1;
        indg[6] = ind0 - num_grids + 1;
        indg[7] = ind0 - num_grids - 1;

        for (q = 0; q < 8; q++) {
          F[ind0] += 2 * beta * wg[q];
          G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
        }
      }
    }
  }

  /// top
  wg = slices.kWeightEdges;
  for (k = 0; k < num_slices; k++) {
    for (m = 1; m < num_grids - 1; m++) {
      ind0 = m + k * num_grids * num_grids;

      indg[0] = ind0 + 1;
      indg[1] = ind0 - 1;
      indg[2] = ind0 + num_grids;
      indg[3] = ind0 + num_grids + 1;
      indg[4] = ind0 + num_grids - 1;

      for (q = 0; q < 5; q++) {
        F[ind0] += 2 * beta * wg[q];
        G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
      }
    }
  }

  // (bottom)
  for (k = 0; k < num_slices; k++) {
    for (m = 1; m < num_grids - 1; m++) {
      ind0 = m + (num_grids - 1) * num_grids + k * num_grids * num_grids;

      indg[0] = ind0 + 1;
      indg[1] = ind0 - 1;
      indg[2] = ind0 - num_grids;
      indg[3] = ind0 - num_grids + 1;
      indg[4] = ind0 - num_grids - 1;

      for (q = 0; q < 5; q++) {
        F[ind0] += 2 * beta * wg[q];
        G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
      }
    }
  }

  // (left)
  for (k = 0; k < num_slices; k++) {
    for (n = 1; n < num_grids - 1; n++) {
      ind0 = n * num_grids + k * num_grids * num_grids;

      indg[0] = ind0 + 1;
      indg[1] = ind0 + num_grids;
      indg[2] = ind0 - num_grids;
      indg[3] = ind0 + num_grids + 1;
      indg[4] = ind0 - num_grids + 1;

      for (q = 0; q < 5; q++) {
        F[ind0] += 2 * beta * wg[q];
        G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
      }
    }
  }

  // (right)
  for (k = 0; k < num_slices; k++) {
    for (n = 1; n < num_grids - 1; n++) {
      ind0 = (num_grids - 1) + n * num_grids + k * num_grids * num_grids;

      indg[0] = ind0 - 1;
      indg[1] = ind0 + num_grids;
      indg[2] = ind0 - num_grids;
      indg[3] = ind0 + num_grids - 1;
      indg[4] = ind0 - num_grids - 1;

      for (q = 0; q < 5; q++) {
        F[ind0] += 2 * beta * wg[q];
        G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
      }
    }
  }

  // (top-left)
  wg = slices.kWeightCorners;
  for (k = 0; k < num_slices; k++) {
    ind0 = k * num_grids * num_grids;

    indg[0] = ind0 + 1;
    indg[1] = ind0 + num_grids;
    indg[2] = ind0 + num_grids + 1;

    for (q = 0; q < 3; q++) {
      F[ind0] += 2 * beta * wg[q];
      G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
    }
  }

  // (top-right)
  for (k = 0; k < num_slices; k++) {
    ind0 = (num_grids - 1) + k * num_grids * num_grids;

    indg[0] = ind0 - 1;
    indg[1] = ind0 + num_grids;
    indg[2] = ind0 + num_grids - 1;

    for (q = 0; q < 3; q++) {
      F[ind0] += 2 * beta * wg[q];
      G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
    }
  }

  // (bottom-left)
  for (k = 0; k < num_slices; k++) {
    ind0 = (num_grids - 1) * num_grids + k * num_grids * num_grids;

    indg[0] = ind0 + 1;
    indg[1] = ind0 - num_grids;
    indg[2] = ind0 - num_grids + 1;

    for (q = 0; q < 3; q++) {
      F[ind0] += 2 * beta * wg[q];
      G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
    }
  }

  // (bottom-right)
  for (k = 0; k < num_slices; k++) {
    ind0 = (num_grids - 1) + (num_grids - 1) * num_grids +
      k * num_grids * num_grids;

    indg[0] = ind0 - 1;
    indg[1] = ind0 - num_grids;
    indg[2] = ind0 - num_grids - 1;

    for (q = 0; q < 3; q++) {
      F[ind0] += 2 * beta * wg[q];
      G[ind0] -= 2 * beta * wg[q] * (recon[ind0] + recon[indg[q]]);
    }
  }

  int count = num_grids*num_grids;
  for (i=0; i<num_slices; i++) {
    float *suma = &reduction_objects()[i][0];
    for (int j=0; j<count; j++) {
      G[i*count+j] += suma[j*2+1];
    }
  }
}

void PMLReconSpace::UpdateRecon(
    TraceData &trace_data,
    DataRegion2DBareBase<float> &comb_replica)  // Locally combined replica
{
  PMLDataRegion &slices = dynamic_cast<PMLDataRegion&>(trace_data.sinograms());
  auto &recon = trace_data.metadata().recon();
  size_t rows = comb_replica.rows();
  size_t cols = comb_replica.cols()/2;

  float *F = slices.F();
  float *G = slices.G();

  for(size_t i=0; i<rows; ++i){
    auto replica = comb_replica[i];
    for(size_t j=0; j<cols; ++j){
      size_t index = (i*cols) + j;
      recon[index] =
        (-G[index] + sqrt(G[index]*G[index] - 8*replica[j*2]*F[index])) /
        (4*F[index]);
    }
  }
}

void PMLReconSpace::PartialBackProjection()
{
  UpdateReconReplica(
    reconparams.simdata,
    (*reconparams.rays)[reconparams.curr_col],
    reconparams.recon,
    reconparams.curr_slice,
    reconparams.indi,
    reconparams.norms,
    reconparams.len);
}


void PMLReconSpace::UpdateReconReplica(
    float simdata,
    float ray,
    float *recon,
    int curr_slice,
    int const * const indi,
    float *norms, 
    int len)
{
  auto &slice_t = reduction_objects()[curr_slice];
  auto slice = &slice_t[0];

  float upd = ray/simdata;
  for (int i=0; i <len-1; ++i) {
    size_t index = indi[i]*2;
    slice[index] -= recon[indi[i]]*norms[i]*upd;
    slice[index+1] += norms[i];
  }
}


