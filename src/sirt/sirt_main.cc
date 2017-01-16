#include "mpi.h"
#include "trace_h5io.h"
#include "tclap/CmdLine.h"
#include "disp_comm_mpi.h"
#include "data_region_base.h"
#include "disp_engine_reduction.h"
#include "sirt.h"

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
        TCLAP::CmdLine cmd("SIRT Iterative Image Reconstruction", ' ', "0.01");
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

        if(rank==0)
        {
          std::cout << "MPI rank:"<< rank << "; MPI size:" << size << std::endl;
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
  #ifdef TIMERON
  std::chrono::duration<double> read_tot(0.);
  auto read_beg = std::chrono::system_clock::now();
  #endif
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
  #ifdef TIMERON
  read_tot += (std::chrono::system_clock::now()-read_beg);
  #endif
  /* Convert degree values to radian */
  trace_utils::DegreeToRadian(*theta);

  /* Setup metadata data structure */
  // INFO: TraceMetadata destructor frees theta->data!
  // TraceMetadata internally creates reconstruction object
  TraceMetadata trace_metadata(
      static_cast<float *>(theta->data),  /// float const *theta,
      0,                                  /// int const proj_id,
      beg_index,                          /// int const slice_id,
      0,                                  /// int const col_id,
      input_slice->metadata->dims[1],     /// int const num_tot_slices,
      input_slice->metadata->dims[0],     /// int const num_projs,
      n_blocks,                           /// int const num_slices,
      input_slice->metadata->dims[2],     /// int const num_cols,
      input_slice->metadata->dims[2],     /// int const num_grids,
      config.center);         /// float const center

  // INFO: DataRegionBase destructor deletes input_slice.data pointer
  ADataRegion<float> *slices = 
    new DataRegionBase<float, TraceMetadata>(
        static_cast<float *>(input_slice->data),
        trace_metadata.count(),
        &trace_metadata);

  /***********************/
  /* Initiate middleware */
  /* Prepare main reduction space and its objects */
  /* The size of the reconstruction object (in reconstruction space) is
   * twice the reconstruction object size, because of the length storage
   */
  auto main_recon_space = new SIRTReconSpace(
      n_blocks, 
      2*trace_metadata.num_cols()*trace_metadata.num_cols());
  main_recon_space->Initialize(trace_metadata.num_grids());
  auto &main_recon_replica = main_recon_space->reduction_objects();
  float init_val=0.;
  main_recon_replica.ResetAllItems(init_val);

  /* Prepare processing engine and main reduction space for other threads */
  DISPEngineBase<SIRTReconSpace, float> *engine =
    new DISPEngineReduction<SIRTReconSpace, float>(
        comm,
        main_recon_space,
        config.thread_count);
        /// # threads (0 for auto assign the number of threads)

  /**********************/

  /**************************/
  /* Perform reconstruction */
  /* Define job size per thread request */
  int64_t req_number = trace_metadata.num_cols();

  #ifdef TIMERON
  std::chrono::duration<double> recon_tot(0.), inplace_tot(0.), update_tot(0.);
  #endif
  for(int i=0; i<config.iteration; ++i){
    /// std::cout << "Iteration: " << i << std::endl;
    #ifdef TIMERON
    auto recon_beg = std::chrono::system_clock::now();
    #endif
    engine->RunParallelReduction(*slices, req_number);  /// Reconstruction
    #ifdef TIMERON
    recon_tot += (std::chrono::system_clock::now()-recon_beg);
    auto inplace_beg = std::chrono::system_clock::now();
    #endif
    engine->ParInPlaceLocalSynchWrapper();              /// Local combination
    #ifdef TIMERON
    inplace_tot += (std::chrono::system_clock::now()-inplace_beg);

    /// Update reconstruction object
    auto update_beg = std::chrono::system_clock::now();
    #endif
    main_recon_space->UpdateRecon(trace_metadata.recon(), main_recon_replica);
    #ifdef TIMERON
    update_tot += (std::chrono::system_clock::now()-update_beg);

    if(comm->rank()==0){
      std::cout << "Reconstruction time=" << recon_tot.count() <<
      std::cout << "; Local combination time=" << inplace_tot.count() << 
      std::cout << "; Update time=" << update_tot.count() <<
      std::cout << "; Read time=" << read_tot.count() << std::endl;
    }
    #endif
    
    /// Reset iteration
    engine->ResetReductionSpaces(init_val);
    slices->ResetMirroredRegionIter();
  }
  /**************************/

  /* Write reconstructed data to disk */
  #ifdef TIMERON
  std::chrono::duration<double> write_tot(0.);
  auto write_beg = std::chrono::system_clock::now();
  #endif
  trace_io::WriteRecon(
      trace_metadata, *d_metadata, 
      config.kReconOutputPath, 
      config.kReconDatasetPath);
  #ifdef TIMERON
  write_tot += (std::chrono::system_clock::now()-write_beg);

  if(comm->rank()==0){
    std::cout << "Reconstruction time=" << recon_tot.count() << std::endl;
    std::cout << "Local combination time=" << inplace_tot.count() << std::endl;
    std::cout << "Update time=" << update_tot.count() << std::endl;
    std::cout << "Read time=" << read_tot.count() << std::endl;
    std::cout << "Write time=" << write_tot.count() << std::endl;
  }
  #endif

  /* Clean-up the resources */
  delete d_metadata->dims;
  delete d_metadata;
  delete slices;
  delete t_metadata->dims;
  delete t_metadata;
  delete theta;
  delete engine;
  delete input_slice;
  delete comm;
}

