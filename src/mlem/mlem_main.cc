#include <iomanip>
#include "mpi.h"
#include "trace_h5io.h"
#include "tclap/CmdLine.h"
#include "disp_comm_mpi.h"
#include "data_region_base.h"
#include "disp_engine_reduction.h"
#include "mlem.h"
#include "ordered_subset.h"

class TraceRuntimeConfig {
  public:
    std::string kProjectionFilePath;
    std::string kProjectionDatasetPath;
    std::string kThetaFilePath;
    std::string kThetaDatasetPath;
    std::string kReconOutputPath;
    std::string kReconDatasetPath;
    int iteration;
    float center;
    int thread_count;

    TraceRuntimeConfig(int argc, char **argv, int rank, int size){
      try
      {
        TCLAP::CmdLine cmd("MLEM Iterative Image Reconstruction", ' ', "0.01");
        TCLAP::ValueArg<std::string> argProjectionFilePath(
          "f", "projectionFilePath", "Projection file path", true, "", "string");
        TCLAP::ValueArg<std::string> argProjectionDatasetPath(
          "d", "projectionDatasetPath", "Projection dataset path in hdf5", false,
          "/exchange/data", "string");
        TCLAP::ValueArg<std::string> argThetaFilePath(
          "e", "thetaFilePath", "Theta file path", true, "", "string");
        TCLAP::ValueArg<std::string> argThetaDatasetPath(
          "q", "thetaDatasetPath", "Theta dataset path", false, "/exchange/theta",
          "string");
        TCLAP::ValueArg<std::string> argReconOutputPath(
          "o", "reconOutputPath", "Output file path for reconstructed image (hdf5)",
          false, "./output.h5", "string");
        TCLAP::ValueArg<std::string> argReconDatasetPath(
          "r", "reconDatasetPath", "Reconstruction dataset path in hdf5 file",
          false, "/data", "string");
        TCLAP::ValueArg<int> argIteration(
          "i", "iteration", "Number of iterations", true, 0, "int");
        TCLAP::ValueArg<float> argCenter(
          "c", "center", "Center value", false, 0., "float");
        TCLAP::ValueArg<int> argThreadCount(
          "t", "thread", "Number of threads per process", false, 1, "int");

        cmd.add(argProjectionFilePath);
        cmd.add(argProjectionDatasetPath);
        cmd.add(argThetaFilePath);
        cmd.add(argThetaDatasetPath);
        cmd.add(argReconOutputPath);
        cmd.add(argReconDatasetPath);
        cmd.add(argIteration);
        cmd.add(argCenter);
        cmd.add(argThreadCount);

        cmd.parse(argc, argv);
        kProjectionFilePath = argProjectionFilePath.getValue();
        kProjectionDatasetPath = argProjectionDatasetPath.getValue();
        kThetaFilePath = argThetaFilePath.getValue();
        kThetaDatasetPath = argThetaDatasetPath.getValue();
        kReconOutputPath = argReconOutputPath.getValue();
        kReconDatasetPath = argReconDatasetPath.getValue();
        iteration = argIteration.getValue();
        center = argCenter.getValue();
        thread_count = argThreadCount.getValue();

        std::cout << "MPI rank:"<< rank << "; MPI size:" << size << std::endl;
        if(rank==0)
        {
          std::cout << "Projection file path=" << kProjectionFilePath << std::endl;
          std::cout << "Projection dataset path in hdf5=" << kProjectionDatasetPath << std::endl;
          std::cout << "Theta file path=" << kThetaFilePath << std::endl;
          std::cout << "Theta dataset path=" << kThetaDatasetPath << std::endl;
          std::cout << "Output file path=" << kReconOutputPath << std::endl;
          std::cout << "Recon. dataset path=" << kReconDatasetPath << std::endl;
          std::cout << "Number of iterations=" << iteration << std::endl;
          std::cout << "Center value=" << center << std::endl;
          std::cout << "Number of threads per process=" << thread_count << std::endl;
        }
      }
      catch (TCLAP::ArgException &e)
      {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
      }
    }
};

int main(int argc, char **argv)
{
  /* Initiate middleware's communication layer */
  DISPCommBase<float> *comm =
        new DISPCommMPI<float>(&argc, &argv);
  TraceRuntimeConfig config(argc, argv, comm->rank(), comm->size());

  /* Read slice data and setup job information */
  auto d_metadata = trace_io::ReadMetadata(
        config.kProjectionFilePath.c_str(), 
        config.kProjectionDatasetPath.c_str());
  int beg_index, n_blocks;
  trace_io::DistributeSlices(
      comm->rank(), comm->size(), 
      d_metadata->dims[1], beg_index, n_blocks);
  auto input_slice = 
    trace_io::ReadSlices(d_metadata, beg_index, n_blocks, 0);

  /* Read theta data */
  auto t_metadata = trace_io::ReadMetadata(
        config.kThetaFilePath.c_str(), 
        config.kThetaDatasetPath.c_str());
  auto theta = trace_io::ReadTheta(t_metadata);
  /* Convert degree values to radian */
  trace_utils::DegreeToRadian(*theta);
  size_t ray_count = 
    input_slice->metadata->dims[0]*input_slice->count*input_slice->metadata->dims[2]; 
  trace_utils::RemoveNegatives(
      static_cast<float *>(input_slice->data),
      ray_count);
  trace_utils::RemoveAbnormals(
      static_cast<float *>(input_slice->data),
      ray_count);

  int ddims[3];
  for(int i=0; i<3; ++i) {
    ddims[i] = input_slice->metadata->dims[i];
    std::cout << "Dim[" << i <<"]=" << ddims[i] << std::endl;
  }
  DataSubset data_subset(
    static_cast<float*>(input_slice->data), 
    static_cast<float*>(theta->data), 
    ddims,
    n_blocks, beg_index,
    100, config.center, 0, 1., SubsetType::Ordered);            /// Number of subsets

  ///* Setup metadata data structure */
  //// INFO: TraceMetadata destructor frees theta->data!
  //// TraceMetadata internally creates reconstruction object
  //TraceMetadata trace_metadata(
  //    static_cast<float *>(theta->data),  /// float const *theta,
  //    0,                                  /// int const proj_id,
  //    beg_index,                          /// int const slice_id,
  //    0,                                  /// int const col_id,
  //    input_slice->metadata->dims[1],     /// int const num_tot_slices,
  //    input_slice->metadata->dims[0],     /// int const num_projs,
  //    n_blocks,                           /// int const num_slices,
  //    input_slice->metadata->dims[2],     /// int const num_cols,
  //    input_slice->metadata->dims[2],     /// int const num_grids,
  //    config.center,                      /// float const center,
  //    0,                                  /// int const num_neighbor_recon_slices,
  //    1.);                                /// float const recon_init_val

  //// INFO: DataRegionBase destructor deletes input_slice.data pointer
  //ADataRegion<float> *slices = 
  //  new DataRegionBase<float, TraceMetadata>(
  //      static_cast<float *>(input_slice->data),
  //      trace_metadata.count(),
  //      &trace_metadata);

  /***********************/
  /* Initiate middleware */
  /* Prepare main reduction space and its objects */
  /* The size of the reconstruction object (in reconstruction space) is
   * twice the reconstruction object size, because of the length storage
   */
  int ncols = input_slice->metadata->dims[2]; 
  int nprojs = input_slice->metadata->dims[0];
  int nsinogs = n_blocks;
  int ntotsinogs = input_slice->metadata->dims[1];
  int beg_sinog_id = beg_index;

  auto main_recon_space = new MLEMReconSpace(nsinogs, 2*ncols*ncols);
  main_recon_space->Initialize(ncols*ncols);
  auto &main_recon_replica = main_recon_space->reduction_objects();
  float init_val=0.;
  main_recon_replica.ResetAllItems(init_val);

  /* Prepare processing engine and main reduction space for other threads */
  DISPEngineBase<MLEMReconSpace, float> *engine =
    new DISPEngineReduction<MLEMReconSpace, float>(
        comm, main_recon_space, config.thread_count); 
  /**********************/

  /**************************/
  /* Perform reconstruction */
  /* Define job size per thread request */
  int64_t req_number = ncols;

  DataRegionBase<float, TraceMetadata> *slices = nullptr;
  for(int i=0; i<config.iteration; ++i){
    std::cout << "Iteration: " << i << std::endl;

    for(int j=0; j<data_subset.NSubsets(); ++j){
      std::cout << "Subset: " << j << "/" << data_subset.NSubsets() << std::endl;
      slices = data_subset.DataRegionSubset(j);
      std::cout << "Data region was read" << std::endl;
      engine->RunParallelReduction(*slices, req_number);  /// Reconstruction
      std::cout << "Parallel reconstruction was done" << std::endl;
      engine->ParInPlaceLocalSynchWrapper();              /// Local combination
      std::cout << "Inplace local synch was done" << std::endl;
      
      /// Update reconstruction object
      main_recon_space->UpdateRecon(slices->metadata().recon(), main_recon_replica);
      std::cout << "Recon space was updated" << std::endl;

      /* Periodically write to disk */
      std::stringstream iteration_stream;
      iteration_stream << std::setfill('0') << std::setw(6) << i << j;
      std::string outputpath = iteration_stream.str() + "-recon.h5";
      trace_io::WriteRecon(
          slices->metadata(), *d_metadata,
          outputpath,
          config.kReconDatasetPath);
      /******************************/
      std::cout << "Image was written" << std::endl;

      // Reset iteration
      engine->ResetReductionSpaces(init_val);
      std::cout << "Reduction space was cleaned" << std::endl;
      // slices->ResetMirroredRegionIter();
    }
  }
  /**************************/

  /* Write reconstructed data to disk */
  trace_io::WriteRecon(
      slices->metadata(), *d_metadata, 
      config.kReconOutputPath, 
      config.kReconDatasetPath);

  /* Clean-up the resources */
  std::cout << "1" << std::endl;
  delete d_metadata->dims;
  std::cout << "2" << std::endl;
  delete d_metadata;
  std::cout << "3" << std::endl;
  delete slices;
  std::cout << "4" << std::endl;
  delete t_metadata->dims;
  std::cout << "5" << std::endl;
  delete t_metadata;
  std::cout << "6" << std::endl;
  delete theta;
  std::cout << "7" << std::endl;
  delete engine;
  std::cout << "8" << std::endl;
  delete input_slice;
  std::cout << "9" << std::endl;
}

