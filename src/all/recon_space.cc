#include "recon_space.h"

void AReconSpace::Initialize(int n_grids)
{
  reconparams.Initialize(n_grids);
}

void AReconSpace::CopyTo(AReconSpace &target)
{
  target.Initialize(reconparams.num_grids);
}

// Forward projection
float AReconSpace::CalculateSimdata(
    float *recon,
    int len,
    int *indi,
    float *norms) 
{
  float simdata = 0.;
  for(int i=0; i<len-1; ++i){
    simdata += recon[indi[i]]*norms[i];
  }
  return simdata;
}

AReconSpace::AReconSpace(int rows, int cols) : 
  AReductionSpaceBase<AReconSpace, float>(rows, cols) 
{}

AReconSpace::~AReconSpace() { }


void AReconSpace::MergeTrimCoordinates(
    int num_grid,
    float *coordx, float *coordy,
    const float *gridx, const float *gridy,
    int *alen, int *blen,
    float *ax, float *ay,
    float *bx, float *by) 
{
  /// Merge the (coordx, gridy) and (gridx, coordy)
  /// on a single array of points (ax, ay) and trim
  /// the coordinates that are outside the
  /// reconstruction grid. 
  *alen = 0;
  *blen = 0;
  for (int i = 0; i <= num_grid; ++i) {
    if (coordx[i] > gridx[0] &&
        coordx[i] < gridx[num_grid]) {
      ax[*alen] = coordx[i];
      ay[*alen] = gridy[i];
      (*alen)++;
    }
    if (coordy[i] > gridy[0] &&
        coordy[i] < gridy[num_grid]) {
      bx[*blen] = gridx[i];
      by[*blen] = coordy[i];
      (*blen)++;
    }
  }
}


void AReconSpace::CalculateCoordinates(
    int num_grid,
    float xi, float yi, float sinq, float cosq,
    const float *gridx, const float *gridy,
    float *coordx, float *coordy)
{
  float srcx, srcy, detx, dety;
  float slope, islope;
  int n;

  /// Find the corresponding source and 
  /// detector locations for a given line
  /// trajectory of a projection (Projection
  /// is specified by sinp and cosp).  
  
  srcx = xi*cosq-yi*sinq;
  srcy = xi*sinq+yi*cosq;
  detx = -1 * (xi*cosq+yi*sinq);
  dety = -xi*sinq+yi*cosq;

  /// Find the intersection points of the
  /// line connecting the source and the detector
  /// points with the reconstruction grid. The 
  /// intersection points are then defined as: 
  /// (coordx, gridy) and (gridx, coordy)
  slope = (srcy-dety)/(srcx-detx);
  islope = 1/slope;
  for (n = 0; n <= num_grid; n++) {
    coordx[n] = islope*(gridy[n]-srcy)+srcx;
    coordy[n] = slope*(gridx[n]-srcx)+srcy;
  }
}

void AReconSpace::SortIntersectionPoints(
    int ind_cond,
    int alen, int blen,
    float *ax, float *ay,
    float *bx, float *by,
    float *coorx, float *coory)
{
  int i=0, j=0, k=0;
  int a_ind;

  for(int i=0, j=0; i<alen; ++i){
    a_ind = (ind_cond) ? i : (alen-1-i);
    coorx[j] = ax[a_ind];
    coory[j] = ay[a_ind];
    i++; j++;
  }

  while (i < alen && j < blen)
  {
    a_ind = (ind_cond) ? i : (alen-1-i);
    if (ax[a_ind] < bx[j]) {
      coorx[k] = ax[a_ind];
      coory[k] = ay[a_ind];
      i++;
      k++;
    } else {
      coorx[k] = bx[j];
      coory[k] = by[j];
      j++;
      k++;
    }
  }
  while (i < alen) {
    a_ind = (ind_cond) ? i : (alen-1-i);
    coorx[k] = ax[a_ind];
    coory[k] = ay[a_ind];
    i++;
    k++;
  }
  while (j < blen) {
    coorx[k] = bx[j];
    coory[k] = by[j];
    j++;
    k++;
  }
}

void AReconSpace::CalculateDistanceLengths(
    int len, int num_grids,
    float *coorx, float *coory,
    float *leng2, int *indi)
{
  int x1, x2, i1, i2;
  float diffx, diffy, midx, midy;
  int indx, indy;

  float mgrids = num_grids/2.;

  for (int i=0; i<len-1; ++i) {
    diffx = coorx[i+1]-coorx[i];
    midx = (coorx[i+1]+coorx[i])/2.;
    diffy = coory[i+1]-coory[i];
    midy = (coory[i+1]+coory[i])/2.;
    leng2[i] = diffx*diffx + diffy*diffy;

    x1 = midx + mgrids;
    i1 = static_cast<int>(x1);
    indx = i1 - (i1>x1);
    x2 = midy + mgrids;
    i2 = static_cast<int>(x2);
    indy = i2 - (i2>x2);
    indi[i] = indx+(indy*num_grids);
  }
}



void AReconSpace::ComputeIndexDistance()
{
      /// Calculate coordinates
      float xi = -1e6;
      float yi = (1-reconparams.num_cols)/2. + 
                  reconparams.curr_col + reconparams.mov;
      CalculateCoordinates(
          reconparams.num_grids, 
          xi, yi, reconparams.sinq, reconparams.cosq, 
          reconparams.gridx, reconparams.gridy, 
          reconparams.coordx, reconparams.coordy);  /// Outputs coordx and coordy

      int alen, blen;
      MergeTrimCoordinates(
          reconparams.num_grids, 
          reconparams.coordx, reconparams.coordy, 
          reconparams.gridx, reconparams.gridy, 
          &alen, &blen, 
          reconparams.ax, reconparams.ay, reconparams.bx, reconparams.by);

      SortIntersectionPoints(
          reconparams.quadrant, 
          alen, blen, 
          reconparams.ax, reconparams.ay, reconparams.bx, reconparams.by, 
          reconparams.coorx, reconparams.coory);

      reconparams.len = alen + blen;
      CalculateDistanceLengths(
          reconparams.len, 
          reconparams.num_grids, 
          reconparams.coorx, reconparams.coory, 
          reconparams.leng2,
          reconparams.indi);
}

void AReconSpace::ForwardProjection()
{
  /// Default normalization values are sqrt of leng2 array
  for(int i=0; i<reconparams.len; ++i) 
    reconparams.norms[i] = sqrt(reconparams.leng2[i]);

  reconparams.simdata = 
    CalculateSimdata(
        reconparams.recon, reconparams.len, 
        reconparams.indi, reconparams.norms);
}

void AReconSpace::Reduce(MirroredRegionBareBase<float> &input)
{
  reconparams.InitInputParams(input);

  /// Partial reconstruction loop
  for (
    int proj = reconparams.curr_proj;
    proj<=(reconparams.curr_proj+reconparams.count_projs);
    ++proj) 
  {
    reconparams.InitLoopParams(proj);

    for (reconparams.curr_col=0; 
         reconparams.curr_col<reconparams.num_cols; 
         ++reconparams.curr_col) 
    {
      ComputeIndexDistance();
      ForwardProjection();
      PartialBackProjection();
    }
  }
}


