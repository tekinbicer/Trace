#include "apmlr.h"

APMLRDataRegion::APMLRDataRegion(
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

APMLRDataRegion::APMLRDataRegion(DataRegionBase<float, TraceMetadata> &region) :
  APMLRDataRegion(&region[0], region.count(), &region.metadata())
{ }

APMLRDataRegion::~APMLRDataRegion(){
  delete [] F_;
  delete [] G_;
}

void APMLRDataRegion::SetFG(float val)
{
  size_t count = 
    metadata().num_slices()*metadata().num_grids()*metadata().num_grids();
  for(size_t i=0; i<count; ++i){
    F_[i]=val;
    G_[i]=val;
  }
}

float* APMLRDataRegion::F() const { return F_; };
float* APMLRDataRegion::G() const { return G_; };

APMLRReconSpace* APMLRReconSpace::Clone()
{
  auto &red_objs = reduction_objects();

  APMLRReconSpace *cloned_obj = new APMLRReconSpace(red_objs.rows(), red_objs.cols());
  (*cloned_obj).reduction_objects() = red_objs;

  static_cast<APMLRReconSpace*>(this)->CopyTo(*cloned_obj);

  return cloned_obj;
}





void APMLRReconSpace::UpdateRecon(
    TraceData &trace_data,                    // Input slices, metadata, recon
    DataRegion2DBareBase<float> &comb_replica)  // Locally combined replica
{
  APMLRDataRegion &slices = dynamic_cast<APMLRDataRegion&>(trace_data.sinograms());
  int slice_count = 
    trace_data.metadata().num_neighbor_recon_slices() *
    trace_data.metadata().num_grids() * trace_data.metadata().num_grids();

  float *recon = &(trace_data.metadata().recon()[slice_count]);
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

void APMLRReconSpace::CalculateFGInner(
  float * recon, float *F, float *G, 
  float beta, float beta1, float delta, float delta1, float regw, 
  int num_slices, int num_grids)
{
  int k, n, m, q;
  int ind0, indg[10];
  float rg[10];
  float mg[10];
  float wg[10];
  float gammag[10];
   
  // Weights for inner neighborhoods.
  float totalwg = 4+4/sqrt(2)+2*regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/totalwg;
  wg[3] = 1/totalwg;
  wg[4] = 1/sqrt(2)/totalwg;
  wg[5] = 1/sqrt(2)/totalwg;
  wg[6] = 1/sqrt(2)/totalwg;
  wg[7] = 1/sqrt(2)/totalwg;
  wg[8] = regw/totalwg;
  wg[9] = regw/totalwg;

  // (inner region)
  //for (k = 1; k < num_slices-1; k++) {
  for (k = 0; k < num_slices; k++) {
    for (n = 1; n < num_grids-1; n++) {
      for (m = 1; m < num_grids-1; m++) {
        ind0 = m + n*num_grids + k*num_grids*num_grids;

        indg[0] = ind0+1;
        indg[1] = ind0-1;
        indg[2] = ind0+num_grids;
        indg[3] = ind0-num_grids;
        indg[4] = ind0+num_grids+1; 
        indg[5] = ind0+num_grids-1;
        indg[6] = ind0-num_grids+1;
        indg[7] = ind0-num_grids-1;
        indg[8] = ind0+num_grids*num_grids;
        indg[9] = ind0-num_grids*num_grids;

        for (q = 0; q < 8; q++) {
          mg[q] = recon[ind0]+recon[indg[q]];
          rg[q] = recon[ind0]-recon[indg[q]];
          gammag[q] = 1/(1+fabs(rg[q]/delta));
          F[ind0] += 2*beta*wg[q]*gammag[q];
          G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
        }

        for (q = 8; q < 10; q++) {
          mg[q] = recon[ind0]+recon[indg[q]];
          rg[q] = recon[ind0]-recon[indg[q]];
          gammag[q] = 1/(1+fabs(rg[q]/delta1));
          F[ind0] += 2*beta1*wg[q]*gammag[q];
          G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
        }
      }
    }
  }

  // Weights for edges.
  totalwg = 3+2/sqrt(2)+2*regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/totalwg;
  wg[3] = 1/sqrt(2)/totalwg;
  wg[4] = 1/sqrt(2)/totalwg;
  wg[5] = regw/totalwg;
  wg[6] = regw/totalwg;

  // (top)
  //for (k = 1; k < num_slices-1; k++) {
  for (k = 0; k < num_slices; k++) {
    for (m = 1; m < num_grids-1; m++) {
      ind0 = m + k*num_grids*num_grids;

      indg[0] = ind0+1;
      indg[1] = ind0-1;
      indg[2] = ind0+num_grids;
      indg[3] = ind0+num_grids+1; 
      indg[4] = ind0+num_grids-1;
      indg[5] = ind0+num_grids*num_grids;
      indg[6] = ind0-num_grids*num_grids;

      for (q = 0; q < 5; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta));
        F[ind0] += 2*beta*wg[q]*gammag[q];
        G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
      }

      for (q = 5; q < 7; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta1));
        F[ind0] += 2*beta1*wg[q]*gammag[q];
        G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
      }
    }
  }

  // (bottom)
  //for (k = 1; k < num_slices-1; k++) {
  for (k = 0; k < num_slices; k++) {
    for (m = 1; m < num_grids-1; m++) {
      ind0 = m + (num_grids-1)*num_grids + k*num_grids*num_grids;

      indg[0] = ind0+1;
      indg[1] = ind0-1;
      indg[2] = ind0-num_grids;
      indg[3] = ind0-num_grids+1;
      indg[4] = ind0-num_grids-1;
      indg[5] = ind0+num_grids*num_grids;
      indg[6] = ind0-num_grids*num_grids;

      for (q = 0; q < 5; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta));
        F[ind0] += 2*beta*wg[q]*gammag[q];
        G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
      }

      for (q = 5; q < 7; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta1));
        F[ind0] += 2*beta1*wg[q]*gammag[q];
        G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
      }
    }
  }

  // (left)  
  //for (k = 1; k < num_slices-1; k++) {
  for (k = 0; k < num_slices; k++) {
    for (n = 1; n < num_grids-1; n++) {
      ind0 = n*num_grids + k*num_grids*num_grids;

      indg[0] = ind0+1;
      indg[1] = ind0+num_grids;
      indg[2] = ind0-num_grids;
      indg[3] = ind0+num_grids+1; 
      indg[4] = ind0-num_grids+1;
      indg[5] = ind0+num_grids*num_grids;
      indg[6] = ind0-num_grids*num_grids;

      for (q = 0; q < 5; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta));
        F[ind0] += 2*beta*wg[q]*gammag[q];
        G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
      }

      for (q = 5; q < 7; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta1));
        F[ind0] += 2*beta1*wg[q]*gammag[q];
        G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
      }
    }
  }

  // (right)                
  //for (k = 1; k < num_slices-1; k++) {
  for (k = 0; k < num_slices; k++) {
    for (n = 1; n < num_grids-1; n++) {
      ind0 = (num_grids-1) + n*num_grids + k*num_grids*num_grids;

      indg[0] = ind0-1;
      indg[1] = ind0+num_grids;
      indg[2] = ind0-num_grids;
      indg[3] = ind0+num_grids-1;
      indg[4] = ind0-num_grids-1;
      indg[5] = ind0+num_grids*num_grids;
      indg[6] = ind0-num_grids*num_grids;

      for (q = 0; q < 5; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta));
        F[ind0] += 2*beta*wg[q]*gammag[q];
        G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
      }

      for (q = 5; q < 7; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta1));
        F[ind0] += 2*beta1*wg[q]*gammag[q];
        G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
      }
    }
  }

  // Weights for corners.
  totalwg = 2+1/sqrt(2)+2*regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/sqrt(2)/totalwg;
  wg[3] = regw/totalwg;
  wg[4] = regw/totalwg;

  // (top-left)
  //for (k = 1; k < num_slices-1; k++) {     
  for (k = 0; k < num_slices; k++) {     
    ind0 = k*num_grids*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0+num_grids;
    indg[2] = ind0+num_grids+1; 
    indg[3] = ind0+num_grids*num_grids;
    indg[4] = ind0-num_grids*num_grids;

    for (q = 0; q < 3; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 3; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (top-right)
  //for (k = 1; k < num_slices-1; k++) {     
  for (k = 0; k < num_slices; k++) {     
    ind0 = (num_grids-1) + k*num_grids*num_grids;

    indg[0] = ind0-1;
    indg[1] = ind0+num_grids;
    indg[2] = ind0+num_grids-1;
    indg[3] = ind0+num_grids*num_grids;
    indg[4] = ind0-num_grids*num_grids;

    for (q = 0; q < 3; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 3; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (bottom-left)
  //for (k = 1; k < num_slices-1; k++) {     
  for (k = 0; k < num_slices; k++) {     
    ind0 = (num_grids-1)*num_grids + k*num_grids*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0-num_grids;
    indg[2] = ind0-num_grids+1;
    indg[3] = ind0+num_grids*num_grids;
    indg[4] = ind0-num_grids*num_grids;

    for (q = 0; q < 3; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 3; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (bottom-right)        
  //for (k = 1; k < num_slices-1; k++) {     
  for (k = 0; k < num_slices; k++) {     
    ind0 = (num_grids-1) + (num_grids-1)*num_grids + k*num_grids*num_grids;

    indg[0] = ind0-1;
    indg[1] = ind0-num_grids;
    indg[2] = ind0-num_grids-1;
    indg[3] = ind0+num_grids*num_grids;
    indg[4] = ind0-num_grids*num_grids;

    for (q = 0; q < 3; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 3; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }
}

void APMLRReconSpace::CalculateFGTop(
   float *recon, float *F, float *G,
   float beta, float beta1, float delta, float delta1, float regw,
   int num_grids)
{
  int n, m, q;
  int ind0, indg[10];
  float rg[10];
  float mg[10];
  float wg[10];
  float gammag[10];

  // Weights for inner neighborhoods.
  float totalwg = 4+4/sqrt(2)+regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/totalwg;
  wg[3] = 1/totalwg;
  wg[4] = 1/sqrt(2)/totalwg;
  wg[5] = 1/sqrt(2)/totalwg;
  wg[6] = 1/sqrt(2)/totalwg;
  wg[7] = 1/sqrt(2)/totalwg;
  wg[8] = regw/totalwg;

  // Upper-most (only for rank=0)
  // (inner region)
  for (n = 1; n < num_grids-1; n++) {
    for (m = 1; m < num_grids-1; m++) {
      ind0 = m + n*num_grids;

      indg[0] = ind0+1;
      indg[1] = ind0-1;
      indg[2] = ind0+num_grids;
      indg[3] = ind0-num_grids;
      indg[4] = ind0+num_grids+1; 
      indg[5] = ind0+num_grids-1;
      indg[6] = ind0-num_grids+1;
      indg[7] = ind0-num_grids-1;
      indg[8] = ind0+num_grids*num_grids;

      for (q = 0; q < 8; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta));
        F[ind0] += 2*beta*wg[q]*gammag[q];
        G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
      }

      for (q = 8; q < 9; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta1));
        F[ind0] += 2*beta1*wg[q]*gammag[q];
        G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
      }
    }
  }

  // Weights for edges.
  totalwg = 3+2/sqrt(2)+regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/totalwg;
  wg[3] = 1/sqrt(2)/totalwg;
  wg[4] = 1/sqrt(2)/totalwg;
  wg[5] = regw/totalwg;

  // (top)
  for (m = 1; m < num_grids-1; m++) {
    ind0 = m;

    indg[0] = ind0+1;
    indg[1] = ind0-1;
    indg[2] = ind0+num_grids;
    indg[3] = ind0+num_grids+1; 
    indg[4] = ind0+num_grids-1;
    indg[5] = ind0+num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (bottom)
  for (m = 1; m < num_grids-1; m++) {
    ind0 = m + (num_grids-1)*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0-1;
    indg[2] = ind0-num_grids;
    indg[3] = ind0-num_grids+1;
    indg[4] = ind0-num_grids-1;
    indg[5] = ind0+num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (left)  
  for (n = 1; n < num_grids-1; n++) {
    ind0 = n*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0+num_grids;
    indg[2] = ind0-num_grids;
    indg[3] = ind0+num_grids+1; 
    indg[4] = ind0-num_grids+1;
    indg[5] = ind0+num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (right)                
  for (n = 1; n < num_grids-1; n++) {
    ind0 = (num_grids-1) + n*num_grids;

    indg[0] = ind0-1;
    indg[1] = ind0+num_grids;
    indg[2] = ind0-num_grids;
    indg[3] = ind0+num_grids-1;
    indg[4] = ind0-num_grids-1;
    indg[5] = ind0+num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // Weights for corners.
  totalwg = 2+1/sqrt(2)+regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/sqrt(2)/totalwg;
  wg[3] = regw/totalwg;

  // (top-left)   
  ind0 = 0;

  indg[0] = ind0+1;
  indg[1] = ind0+num_grids;
  indg[2] = ind0+num_grids+1; 
  indg[3] = ind0+num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }

  // (top-right)    
  ind0 = (num_grids-1);

  indg[0] = ind0-1;
  indg[1] = ind0+num_grids;
  indg[2] = ind0+num_grids-1;
  indg[3] = ind0+num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }

  // (bottom-left) 
  ind0 = (num_grids-1)*num_grids;

  indg[0] = ind0+1;
  indg[1] = ind0-num_grids;
  indg[2] = ind0-num_grids+1;
  indg[3] = ind0+num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }

  // (bottom-right)         
  ind0 = (num_grids-1) + (num_grids-1)*num_grids;

  indg[0] = ind0-1;
  indg[1] = ind0-num_grids;
  indg[2] = ind0-num_grids-1;
  indg[3] = ind0+num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }
}

void APMLRReconSpace::CalculateFGBottom(
    float *recon, float *F, float *G,
    float beta, float beta1, float delta, float delta1, float regw, 
    int num_grids)
{
  int n, m, q;
  int ind0, indg[10];
  float rg[10];
  float mg[10];
  float wg[10];
  float gammag[10];

  // Weights for inner neighborhoods.
  float totalwg = 4+4/sqrt(2)+regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/totalwg;
  wg[3] = 1/totalwg;
  wg[4] = 1/sqrt(2)/totalwg;
  wg[5] = 1/sqrt(2)/totalwg;
  wg[6] = 1/sqrt(2)/totalwg;
  wg[7] = 1/sqrt(2)/totalwg;
  wg[8] = regw/totalwg;

  // Down-most (only for rank=mpi_size)
  // (inner region)
  for (n = 1; n < num_grids-1; n++) {
    for (m = 1; m < num_grids-1; m++) {
      ind0 = m + n*num_grids; //+ (num_slices-1)*num_grids*num_grids;

      indg[0] = ind0+1;
      indg[1] = ind0-1;
      indg[2] = ind0+num_grids;
      indg[3] = ind0-num_grids;
      indg[4] = ind0+num_grids+1; 
      indg[5] = ind0+num_grids-1;
      indg[6] = ind0-num_grids+1;
      indg[7] = ind0-num_grids-1;
      indg[8] = ind0-num_grids*num_grids;

      for (q = 0; q < 8; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta));
        F[ind0] += 2*beta*wg[q]*gammag[q];
        G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
      }

      for (q = 8; q < 9; q++) {
        mg[q] = recon[ind0]+recon[indg[q]];
        rg[q] = recon[ind0]-recon[indg[q]];
        gammag[q] = 1/(1+fabs(rg[q]/delta1));
        F[ind0] += 2*beta1*wg[q]*gammag[q];
        G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
      }
    }
  }

  // Weights for edges.
  totalwg = 3+2/sqrt(2)+regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/totalwg;
  wg[3] = 1/sqrt(2)/totalwg;
  wg[4] = 1/sqrt(2)/totalwg;
  wg[5] = regw/totalwg;

  // (top)
  for (m = 1; m < num_grids-1; m++) {
    ind0 = m; //+ (num_slices-1)*num_grids*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0-1;
    indg[2] = ind0+num_grids;
    indg[3] = ind0+num_grids+1; 
    indg[4] = ind0+num_grids-1;
    indg[5] = ind0-num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (bottom)
  for (m = 1; m < num_grids-1; m++) {
    ind0 = m + (num_grids-1)*num_grids; //+ (num_slices-1)*num_grids*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0-1;
    indg[2] = ind0-num_grids;
    indg[3] = ind0-num_grids+1;
    indg[4] = ind0-num_grids-1;
    indg[5] = ind0-num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (left)  
  for (n = 1; n < num_grids-1; n++) {
    ind0 = n*num_grids; //+ (num_slices-1)*num_grids*num_grids;

    indg[0] = ind0+1;
    indg[1] = ind0+num_grids;
    indg[2] = ind0-num_grids;
    indg[3] = ind0+num_grids+1; 
    indg[4] = ind0-num_grids+1;
    indg[5] = ind0-num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // (right)                
  for (n = 1; n < num_grids-1; n++) {
    ind0 = (num_grids-1) + n*num_grids; //+ (num_slices-1)*num_grids*num_grids;

    indg[0] = ind0-1;
    indg[1] = ind0+num_grids;
    indg[2] = ind0-num_grids;
    indg[3] = ind0+num_grids-1;
    indg[4] = ind0-num_grids-1;
    indg[5] = ind0-num_grids*num_grids;

    for (q = 0; q < 5; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta));
      F[ind0] += 2*beta*wg[q]*gammag[q];
      G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
    }

    for (q = 5; q < 6; q++) {
      mg[q] = recon[ind0]+recon[indg[q]];
      rg[q] = recon[ind0]-recon[indg[q]];
      gammag[q] = 1/(1+fabs(rg[q]/delta1));
      F[ind0] += 2*beta1*wg[q]*gammag[q];
      G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
    }
  }

  // Weights for corners.
  totalwg = 2+1/sqrt(2)+regw;
  wg[0] = 1/totalwg;
  wg[1] = 1/totalwg;
  wg[2] = 1/sqrt(2)/totalwg;
  wg[3] = regw/totalwg;

  // (top-left)
  ind0 = 0;// +(num_slices-1)*num_grids*num_grids;

  indg[0] = ind0+1;
  indg[1] = ind0+num_grids;
  indg[2] = ind0+num_grids+1; 
  indg[3] = ind0-num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }

  // (top-right)   
  ind0 = (num_grids-1); //+ (num_slices-1)*num_grids*num_grids;

  indg[0] = ind0-1;
  indg[1] = ind0+num_grids;
  indg[2] = ind0+num_grids-1;
  indg[3] = ind0-num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }

  // (bottom-left)   
  ind0 = (num_grids-1)*num_grids; //+ (num_slices-1)*num_grids*num_grids;

  indg[0] = ind0+1;
  indg[1] = ind0-num_grids;
  indg[2] = ind0-num_grids+1;
  indg[3] = ind0-num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }

  // (bottom-right)          
  ind0 = (num_grids-1) + (num_grids-1)*num_grids; //+ (num_slices-1)*num_grids*num_grids;

  indg[0] = ind0-1;
  indg[1] = ind0-num_grids;
  indg[2] = ind0-num_grids-1;
  indg[3] = ind0-num_grids*num_grids;

  for (q = 0; q < 3; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta));
    F[ind0] += 2*beta*wg[q]*gammag[q];
    G[ind0] -= 2*beta*wg[q]*gammag[q]*mg[q];
  }

  for (q = 3; q < 4; q++) {
    mg[q] = recon[ind0]+recon[indg[q]];
    rg[q] = recon[ind0]-recon[indg[q]];
    gammag[q] = 1/(1+fabs(rg[q]/delta1));
    F[ind0] += 2*beta1*wg[q]*gammag[q];
    G[ind0] -= 2*beta1*wg[q]*gammag[q]*mg[q];
  }


}

/* If there is any performance bottleneck here, then turn it 
 * into OpenMP code. */
void APMLRReconSpace::CalculateFG(
    ADataRegion<float> &slices_,
    float beta, float beta1, 
    float delta, float delta1, 
    float regw)
{
  APMLRDataRegion &slices = dynamic_cast<APMLRDataRegion&>(slices_);
  TraceMetadata &metadata = slices.metadata();

  int num_slices = metadata.num_slices();
  int slice_beg_index = metadata.slice_id();    
  int slice_end_index = metadata.slice_id() + num_slices-1;
  /// Final index of the last slice in whole dataset
  int slice_final_index = slice_beg_index + metadata.num_total_slices()-1;
  int num_grids = metadata.num_grids();

  float *F = slices.F();
  float *G = slices.G();

  int recon_slice_data_index =
    metadata.num_neighbor_recon_slices()*
    metadata.num_grids() * metadata.num_grids();
  float *recon = &metadata.recon()[recon_slice_data_index];

  float *init_recon_ptr = recon;
  float *init_F_ptr = F;
  float *init_G_ptr = G;
  int remaining_slices = num_slices;
  if(slice_beg_index==0){
    CalculateFGTop(
        init_recon_ptr, init_F_ptr, init_G_ptr, 
        beta, beta1, delta, delta1, regw, 
        num_grids);
    init_recon_ptr = recon + num_grids*num_grids;
    init_F_ptr = F + num_grids*num_grids;
    init_G_ptr = G + num_grids*num_grids;

    --remaining_slices;
  }
  if(slice_final_index==slice_end_index && slice_final_index!=0){
    float *final_slice_ptr = recon + (num_slices-1)*num_grids*num_grids;
    float *final_F_ptr = F + (num_slices-1)*num_grids*num_grids;
    float *final_G_ptr = G + (num_slices-1)*num_grids*num_grids;
    CalculateFGBottom(
        final_slice_ptr, final_F_ptr, final_G_ptr, 
        beta, beta1, delta, delta1, regw, 
        num_grids);
    --remaining_slices;
  }
  CalculateFGInner(
      init_recon_ptr, init_F_ptr, init_G_ptr, 
      beta, beta1, delta, delta1, regw, 
      remaining_slices, num_grids);

  // Calculate G
  int count = num_grids*num_grids;
  for (int i=0; i<num_slices; i++) {
    float *suma = &reduction_objects()[i][0];
    for (int j=0; j<count; j++) {
      G[i*count+j] += suma[j*2+1];
    }
  }

}

APMLRReconSpace::APMLRReconSpace(int rows, int cols) : 
      AReconSpace(rows, cols) {}


void APMLRReconSpace::UpdateReconReplica(
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


void APMLRReconSpace::PartialBackProjection()
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
