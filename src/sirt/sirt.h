#ifndef DISP_APPS_RECONSTRUCTION_SIRT_SIRT_H
#define DISP_APPS_RECONSTRUCTION_SIRT_SIRT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <chrono>
#include "hdf5.h"
#include "string.h"
#include "trace_data.h"
#include "trace_utils.h"
#include "reduction_space_a.h"
#include "data_region_base.h"

class SIRTReconSpace : 
  public AReductionSpaceBase<SIRTReconSpace, float>
{
  private:
      float *coordx = nullptr;
      float *coordy = nullptr;
      float *ax = nullptr;
      float *ay = nullptr;
      float *bx = nullptr;
      float *by = nullptr;
      float *coorx = nullptr;
      float *coory = nullptr;
      float *leng = nullptr;
      float *leng2 = nullptr;
      int *indi = nullptr;

      int num_grids;

  protected:
    // Forward projection
    float CalculateSimdata(
        float *recon,
        int len,
        int *indi,
        float *leng)
    {
      float simdata = 0.;
      for(int i=0; i<len-1; ++i)
        simdata += recon[indi[i]]*leng[i];
      return simdata;
    }

    void UpdateReconReplica(
        float simdata,
        float ray,
        int curr_slice,
        int *indi,
        float *leng2,
        float *leng, 
        int len,
        int suma_beg_offset)
    {
      float upd, a2=0.;

      auto &slice = reduction_objects()[curr_slice];

      start_replica = std::chrono::system_clock::now();
      for (int i=0; i<len-1; ++i) {
        if (indi[i] >= suma_beg_offset) continue;
        a2 += leng2[i];
        slice[suma_beg_offset + indi[i]] += leng[i];
      }
      end_replica = std::chrono::system_clock::now();
      timer_update_replica_loop1 += end_replica-start_replica;

      upd = (ray-simdata) / a2;
      start_replica = std::chrono::system_clock::now();
      for (int i=0; i <len-1; ++i) {
        if (indi[i] >= suma_beg_offset) continue;
        slice[indi[i]] += leng[i]*upd;
      }
      end_replica = std::chrono::system_clock::now();
      timer_update_replica_loop2 += end_replica-start_replica;
    }

  public:
    SIRTReconSpace(int rows, int cols) : 
      AReductionSpaceBase(rows, cols) {}

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
      leng = new float[2*num_grids];
      leng2 = new float[2*num_grids];
      indi = new int[2*num_grids];
    }

    virtual void CopyTo(SIRTReconSpace &target) {
      target.Initialize(num_grids);
    }
    
    void Finalize(){
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
    }

    virtual ~SIRTReconSpace(){
      PrintProfileInfo();
      Finalize();
    }

    /**********************/
    /* Execution Profiler */
    void PrintProfileInfo(){
      std::cout << 
        "Total time: " << timer_all.count() << std::endl <<
        "--Coordinate calculation: " << timer_coordinates.count() << std::endl <<
        "--Sort intersections: " << timer_sort_int.count() << std::endl <<
        "--Merge trim coordinates: " << timer_merge_trim.count() << std::endl <<
        "--Distance calculation: " << timer_dist.count() << std::endl <<
        "--Simdata calculation: " << timer_simdata.count() << std::endl <<
        "--Updating replica: " << timer_update_replica.count() << std::endl <<
        "---Updating replica loop1: " << timer_update_replica_loop1.count() <<  std::endl <<
        "---Updating replica loop2: " << timer_update_replica_loop2.count() << std::endl;
    }

    std::chrono::time_point<std::chrono::system_clock> 
      start, end, start_all, end_all,
      start_replica, end_replica;
    std::chrono::duration<double> 
      timer_coordinates, timer_sort_int,
      timer_dist, timer_simdata, timer_update_replica,
      timer_update_replica_loop1, timer_update_replica_loop2,
      timer_merge_trim, timer_all;
    /**********************/

    void Reduce(MirroredRegionBareBase<float> &input)
    {
      start_all = std::chrono::system_clock::now();

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

      /*
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
         if(rank==0){
         std::cout << "Ray index=" << rays.index() << std::endl;
         std::cout << "Ray count=" << rays.count() << std::endl;
         std::cout << "Beginning projection=" << curr_proj << std::endl;
         std::cout << "Number of projection=" << count_projs << std::endl;
         }
         */

      /* Reconstruction start */
      for (int proj = curr_proj; proj<=(curr_proj+count_projs); ++proj) {
        float theta_q = theta[proj];
        int quadrant = trace_utils::CalculateQuadrant(theta_q);
        float sinq = sinf(theta_q);
        float cosq = cosf(theta_q);

        int curr_slice = metadata.RaySlice(rays.index());
        int curr_slice_offset = curr_slice*num_grids*num_grids;
        float *recon = (&(metadata.recon()[0])+curr_slice_offset);

        /*
           if(rank==0){
           std::cout << "Current projection=" << proj << std::endl;
           std::cout << "Current projection slice=" << curr_slice << std::endl;
           std::cout << "Current slice offset=" << curr_slice_offset << std::endl;
           std::cout << std::endl;
           }
           */

        for (int curr_col=0; curr_col<num_cols; ++curr_col) {
          /// Calculate coordinates
          float xi = -1e6;
          float yi = (1-num_cols)/2. + curr_col+mov;
          start = std::chrono::system_clock::now();
          trace_utils::CalculateCoordinates(
              num_grids, 
              xi, yi, sinq, cosq, 
              gridx, gridy, 
              coordx, coordy);  /// Outputs coordx and coordy
          end = std::chrono::system_clock::now();
          timer_coordinates += end-start;

          /// Merge the (coordx, gridy) and (gridx, coordy)
          /// Output alen and after
          int alen, blen;
          start = std::chrono::system_clock::now();
          trace_utils::MergeTrimCoordinates(
              num_grids, 
              coordx, coordy, 
              gridx, gridy, 
              &alen, &blen, 
              ax, ay, bx, by);
          end = std::chrono::system_clock::now();
          timer_merge_trim += end-start;

          /// Sort the array of intersection points (ax, ay)
          /// The new sorted intersection points are
          /// stored in (coorx, coory).
          /// if quadrant=1 then a_ind = i; if 0 then a_ind = (alen-1-i)
          start = std::chrono::system_clock::now();
          trace_utils::SortIntersectionPoints(
              quadrant, 
              alen, blen, 
              ax, ay, bx, by, 
              coorx, coory);
          end = std::chrono::system_clock::now();
          timer_sort_int += end-start;

          /// Calculate the distances (leng) between the
          /// intersection points (coorx, coory). Find
          /// the indices of the pixels on the
          /// reconstruction grid (ind_recon).
          int len = alen + blen;
          start = std::chrono::system_clock::now();
          trace_utils::CalculateDistanceLengths(
              len, 
              num_grids, 
              coorx, coory, 
              leng, leng2, 
              indi);
          end = std::chrono::system_clock::now();
          timer_dist += end-start;

          /*******************************************************/
          /* Below is for updating the reconstruction grid and
           * is algorithm specific part.
           */
          /// Forward projection
          start = std::chrono::system_clock::now();
          float simdata = CalculateSimdata(recon, len, indi, leng);
          end = std::chrono::system_clock::now();
          timer_simdata += end-start;

          /// Update recon 
          start = std::chrono::system_clock::now();
          int suma_beg_offset = num_grids*num_grids;
          UpdateReconReplica(
              simdata, 
              rays[curr_col], 
              curr_slice, 
              indi, 
              leng2, leng,
              len, 
              suma_beg_offset);
          end = std::chrono::system_clock::now();
          timer_update_replica += end-start;
          /*******************************************************/
        }
      }
      end_all = std::chrono::system_clock::now();
      timer_all += end_all-start_all;
    }
};

#endif    // DISP_APPS_RECONSTRUCTION_SIRT_SIRT_H
