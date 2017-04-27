#include "trace_ph5io.h"
#include "trace_utils.h"
#include "data_region_bare_base.h"
#include "data_region_base.h"
#include <chrono>

TracePH5IO::TracePH5IO(
    DISPCommBase<float> &kcomm,
    std::string &kOutputFilePath,
    std::string &kOutputDatasetPath,
    std::string &kProjectionFilePath,
    std::string &kProjectionDatasetPath,
    std::string &kThetaFilePath,
    std::string &kThetaDatasetPath,
    float kCenter,
    bool kdegree_to_radian)
    : comm {kcomm},
      outputFilePath {kOutputFilePath},
      outputDatasetPath {kOutputDatasetPath},
      projectionFilePath {kProjectionFilePath},
      projectionDatasetPath {kProjectionDatasetPath},
      thetaFilePath {kThetaFilePath},
      thetaDatasetPath {kThetaDatasetPath},
      center {kCenter},
      degree_to_radian {kdegree_to_radian}
{}


TracePH5IO::TracePH5IO(
  DISPCommBase<float> &kcomm,
  TraceRuntimeConfig &config)
    : TracePH5IO(
        kcomm,
        config.kOutputFilePath,
        config.kOutputDatasetPath,
        config.kProjectionFilePath,
        config.kProjectionDatasetPath,
        config.kThetaFilePath,
        config.kThetaDatasetPath,
        config.center,
        config.degree_to_radian)
{}


TraceData TracePH5IO::Read(){
  TraceData trace_data;
  /* Read slice data and setup job information */
#ifdef TIMERON
  std::chrono::duration<double> read_tot(0.);
  auto read_beg = std::chrono::system_clock::now();
#endif
  auto d_metadata = trace_io::ReadMetadata(
      projectionFilePath.c_str(), 
      projectionDatasetPath.c_str());
  int beg_index, n_blocks;
  trace_io::DistributeSlices(
      comm.rank(), comm.size(), 
      d_metadata->dims[1], beg_index, n_blocks);
  auto input_slice = 
    trace_io::ReadSlices(d_metadata, beg_index, n_blocks, 0);

  /* Read theta data */
  auto t_metadata = trace_io::ReadMetadata(
      thetaFilePath.c_str(), 
      thetaDatasetPath.c_str());
  auto theta = trace_io::ReadTheta(t_metadata);
#ifdef TIMERON
  read_tot += (std::chrono::system_clock::now()-read_beg);
#endif
  /* Convert degree values to radian */
  if(degree_to_radian) trace_utils::DegreeToRadian(*theta);

  /* Setup metadata data structure */
  // INFO: TraceMetadata destructor frees theta->data!
  // TraceMetadata internally creates reconstruction object
  trace_data.metadata( new TraceMetadata(
        static_cast<float *>(theta->data),  /// float const *theta,
        0,                                  /// int const proj_id,
        beg_index,                          /// int const slice_id,
        0,                                  /// int const col_id,
        input_slice->metadata->dims[1],     /// int const num_tot_slices,
        input_slice->metadata->dims[0],     /// int const num_projs,
        n_blocks,                           /// int const num_slices,
        input_slice->metadata->dims[2],     /// int const num_cols,
        input_slice->metadata->dims[2],     /// int const num_grids,
        center));         /// float const center

  // INFO: DataRegionBase destructor deletes input_slice.data pointer
  trace_data.sinograms( 
      new DataRegionBase<float, TraceMetadata>(
        static_cast<float *>(input_slice->data),
        trace_data.metadata().count(),
        &trace_data.metadata()));

  return trace_data;
}

void TracePH5IO::Write(TraceData &trace_data){
  hsize_t ndims = 3;
  hsize_t rank_dims[3] = {
    static_cast<hsize_t>(trace_data.metadata().num_slices()),
    static_cast<hsize_t>(trace_data.metadata().num_cols()),
    static_cast<hsize_t>(trace_data.metadata().num_cols())
  };
  hsize_t app_dims[3] = {
    static_cast<hsize_t>(trace_data.metadata().num_total_slices()),
    static_cast<hsize_t>(trace_data.metadata().num_cols()),
    static_cast<hsize_t>(trace_data.metadata().num_cols())
  };

  ADataRegion<float> &recon = trace_data.metadata().recon();

  int recon_slice_data_index =
    trace_data.metadata().num_neighbor_recon_slices()*
    trace_data.metadata().num_grids() * trace_data.metadata().num_grids();

  trace_io::WriteData(
    &recon[recon_slice_data_index],
    ndims, rank_dims, trace_data.metadata().slice_id(),
    ndims, app_dims, 
    0,
    outputFilePath.c_str(), outputDatasetPath.c_str(),
    MPI_COMM_WORLD, MPI_INFO_NULL, H5FD_MPIO_COLLECTIVE);
}


