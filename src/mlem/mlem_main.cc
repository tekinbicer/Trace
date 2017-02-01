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
    int block_remove_iteration;
    float center;
    int thread_count;
    int subsets = 1;
    int block_remove_subsets = 1;
    int write_freq = 0;
    int write_block_freq = 0;

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
        TCLAP::ValueArg<int> argBlockRemoveIteration(
          "", "block-remove-iteration", "Number of iterations for removing subset blocks on 3D image", false, 0, "int");
        TCLAP::ValueArg<float> argCenter(
          "c", "center", "Center value", false, 0., "float");
        TCLAP::ValueArg<int> argThreadCount(
          "t", "thread", "Number of threads per process", false, 1, "int");
        TCLAP::ValueArg<float> argSubsets(
          "", "subsets", "Ordered subsets", false, 1, "int");
        TCLAP::ValueArg<float> argBlockRemoveSubsets(
          "", "block-remove-subsets", "Number of subsets for removing block lines", false, 1, "int");
        TCLAP::ValueArg<float> argWriteFreq(
          "", "write-frequency", "Write frequency", false, 0, "int");
        TCLAP::ValueArg<float> argWriteBlockFreq(
          "", "write-block-frequency", "Write frequency during subset block cleaning", false, 0, "int");

        cmd.add(argProjectionFilePath);
        cmd.add(argProjectionDatasetPath);
        cmd.add(argThetaFilePath);
        cmd.add(argThetaDatasetPath);
        cmd.add(argReconOutputPath);
        cmd.add(argReconDatasetPath);
        cmd.add(argIteration);
        cmd.add(argBlockRemoveIteration);
        cmd.add(argCenter);
        cmd.add(argThreadCount);
        cmd.add(argSubsets);
        cmd.add(argBlockRemoveSubsets);
        cmd.add(argWriteFreq);
        cmd.add(argWriteBlockFreq);

        cmd.parse(argc, argv);
        kProjectionFilePath = argProjectionFilePath.getValue();
        kProjectionDatasetPath = argProjectionDatasetPath.getValue();
        kThetaFilePath = argThetaFilePath.getValue();
        kThetaDatasetPath = argThetaDatasetPath.getValue();
        kReconOutputPath = argReconOutputPath.getValue();
        kReconDatasetPath = argReconDatasetPath.getValue();
        iteration = argIteration.getValue();
        block_remove_iteration = argBlockRemoveIteration.getValue();
        center = argCenter.getValue();
        thread_count = argThreadCount.getValue();
        subsets = argSubsets.getValue();
        block_remove_subsets = argBlockRemoveSubsets.getValue();
        write_freq= argWriteFreq.getValue();
        write_block_freq= argWriteBlockFreq.getValue();

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
          std::cout << "Number of block remove iterations=" << block_remove_iteration << std::endl;
          std::cout << "Center value=" << center << std::endl;
          std::cout << "Number of threads per process=" << thread_count << std::endl;
          std::cout << "Subsets=" << subsets  << std::endl;
          std::cout << "Block remove subsets=" << block_remove_subsets  << std::endl;
          std::cout << "Write frequency=" << write_freq << std::endl;
          std::cout << "Write block frequency=" << write_block_freq << std::endl;
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
    config.subsets, config.center, 0, 1., SubsetType::Ordered);/// Number of subsets

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


  /**************************************/
  /* Perform reconstruction             */
  /* Define job size per thread request */
  int64_t req_number = ncols;

  DataRegionBase<float, TraceMetadata> *slices = nullptr;
  int counter=0;
  for(int i=0; i<config.iteration; ++i){
    std::cout << "Iteration: " << i << std::endl;

    for(int j=0; j<data_subset.NSubsets(); ++j){
      std::cout << "Subset: " << j << "/" << data_subset.NSubsets() << std::endl;
      slices = data_subset.DataRegionSubset(j);
      std::cout << "Data region was read" << std::endl;
      std::cout << "Parallel reconstruction is about to start: recon=" << &(slices->metadata().recon()[0]) << std::endl;
      engine->RunParallelReduction(*slices, req_number);  /// Reconstruction
      std::cout << "Parallel reconstruction was done" << std::endl;
      engine->ParInPlaceLocalSynchWrapper();              /// Local combination
      std::cout << "Inplace local synch was done" << std::endl;
      
      /// Update reconstruction object
      main_recon_space->UpdateRecon(slices->metadata().recon(), main_recon_replica);
      std::cout << "Recon space was updated" << std::endl;

      /* Periodically write to disk */
      if(config.write_freq>0 && counter%config.write_freq==0){
      std::stringstream iteration_stream;
      iteration_stream << std::setfill('0') << std::setw(6) << i << j;
      std::string outputpath = iteration_stream.str() + "-recon.h5";
      trace_io::WriteRecon(
          slices->metadata(), *d_metadata,
          outputpath,
          config.kReconDatasetPath);
      std::cout << "Image was written" << std::endl;
      }
      ++counter;
      /******************************/

      // Reset iteration
      engine->ResetReductionSpaces(init_val);
      std::cout << "Reduction space was cleaned" << std::endl;
      // slices->ResetMirroredRegionIter();
    }
  }
  /**************************/

  /* Write reconstructed data to disk */
  std::string outputpath = "output-recon.h5";
  trace_io::WriteRecon(
      slices->metadata(), *d_metadata, 
      outputpath, 
      config.kReconDatasetPath);




  /***************************/
  /* Remove block signs      */
  DataSubset data_br_subset(
    static_cast<float*>(input_slice->data), 
    static_cast<float*>(theta->data), 
    ddims,
    n_blocks, beg_index,
    config.block_remove_subsets, config.center, 0, 1., SubsetType::Ordered);/// Number of subsets

  counter=0;
  for(int i=0; i<config.block_remove_iteration; ++i){
    std::cout << "Block remove iteration: " << i << std::endl;

    for(int j=0; j<data_br_subset.NSubsets(); ++j){
      std::cout << "Subset: " << j << "/" << data_subset.NSubsets() << std::endl;
      slices = data_br_subset.DataRegionSubset(j);
      std::cout << "Data region was read" << std::endl;
      engine->RunParallelReduction(*slices, req_number);  /// Reconstruction
      std::cout << "Parallel reconstruction was done" << std::endl;
      engine->ParInPlaceLocalSynchWrapper();              /// Local combination
      std::cout << "Inplace local synch was done" << std::endl;
      
      /// Update reconstruction object
      main_recon_space->UpdateRecon(slices->metadata().recon(), main_recon_replica);
      std::cout << "Recon space was updated" << std::endl;

      /* Periodically write to disk */
      if(config.write_block_freq>0 && counter%config.write_block_freq==0){
      std::stringstream iteration_stream;
      iteration_stream << std::setfill('0') << std::setw(6) << i << j;
      std::string outputpath = iteration_stream.str() + "-block-recon.h5";
      trace_io::WriteRecon(
          slices->metadata(), *d_metadata,
          outputpath,
          config.kReconDatasetPath);
      std::cout << "Image was written" << std::endl;
      }
      ++counter;
      /******************************/

      // Reset iteration
      engine->ResetReductionSpaces(init_val);
      std::cout << "Reduction space was cleaned" << std::endl;
      // slices->ResetMirroredRegionIter();
    }
  }
  /**************************/

  /* Write reconstructed data to disk */
  std::string outputpath_cleaned = "output-cleaned-recon.h5";
  trace_io::WriteRecon(
      slices->metadata(), *d_metadata, 
      outputpath_cleaned, 
      config.kReconDatasetPath);


  /* Clean-up the resources */
  delete d_metadata->dims;
  delete d_metadata;
  delete slices;
  delete t_metadata->dims;
  delete t_metadata;
  delete theta;
  delete engine;
  delete input_slice;
}

