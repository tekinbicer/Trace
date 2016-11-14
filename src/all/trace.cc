#include "mpi.h"
#include "tclap/CmdLine.h"
#include "data_region_base.h"
#include "disp_comm_mpi.h"
#include "disp_engine_reduction.h"
#include "sirt.h"

class TraceRuntimeConfig {
  public:
    std::string kReconstructionAlg;

    std::string kInputFilePath;
    std::string kOutputFilePath;
    std::string kOutputDatasetPath;

    std::string kProjectionFilePath;
    std::string kProjectionDatasetPath;
    std::string kThetaFilePath;
    std::string kThetaDatasetPath;

    int iteration_count;
    int thread_count;
    float center;


    /// Setup these according to application
    float b0=10., b1=1., d0=1., d1=1., regw=1.;
    bool degree_to_radian=false;

    std::shared_ptr<DISPCommBase<float>> comm;

    TraceRuntimeConfig(int argc, char **argv){
      try
      {
        TCLAP::CmdLine 
          cmd("Trace: Iterative Image Reconstruction Engine", '=', "0.01");

        TCLAP::SwitchArg argDegreeToRadian(
          "", "degree-to-radian", "Transform theta values from degree to "
          "radian", cmd, false);
        TCLAP::ValueArg<float> argRegw(
          "", "regw", "Regw 1 value. Default value is 1.", 
          false, 1., "float");
        cmd.add(argRegw);
        TCLAP::ValueArg<float> argBetaOne(
          "", "beta-one", "Beta 1 value. Default value is 1.", 
          false, 1., "float");
        cmd.add(argBetaOne);
        TCLAP::ValueArg<float> argBetaZero(
          "", "beta-zero", "Beta 0 value. Default value is 10.", 
          false, 10., "float");
        cmd.add(argBetaZero);
        TCLAP::ValueArg<float> argDeltaOne(
          "", "delta-one", "Delta 1 value. Default value is 1.", 
          false, 1., "float");
        cmd.add(argDeltaOne);
        TCLAP::ValueArg<float> argDeltaZero(
          "", "delta-zero", "Delta 0 value. Default value is 1.", 
          false, 1., "float");
        cmd.add(argDeltaZero);

        TCLAP::ValueArg<std::string> argProjectionFilePath(
          "p", "projection-file-path", 
          "Input projection file path.", false, "", "string");
        cmd.add(argProjectionFilePath);
        TCLAP::ValueArg<std::string> argProjectionDatasetPath(
          "d", "projection-dataset-path", "Projection dataset path in HDF5"
          "file. Default value is \"/exchange/data\".", 
          false, "/exchange/data", "string");
        cmd.add(argProjectionDatasetPath);
        TCLAP::ValueArg<std::string> argThetaFilePath(
          "e", "theta-file-path", 
          "Input theta file path.", false, "", "string");
        cmd.add(argThetaFilePath);
        TCLAP::ValueArg<std::string> argThetaDatasetPath(
          "q", "theta-dataset-path", "Theta dataset path in HDF5 file. "
          "Default value is \"/exchange/theta\".", 
          false, "/exchange/theta", "string");
        cmd.add(argThetaDatasetPath);
        
        TCLAP::ValueArg<std::string> argOutputDatasetPath(
          "m", "output-dataset-path", 
          "Reconstructed 3D image dataset path in hdf5 file.",
          false, "/image", "string");
        cmd.add(argOutputDatasetPath);

        TCLAP::ValueArg<float> argCenter(
          "c", "center", "Center value for input dataset.", false, 0., "float");
        cmd.add(argCenter);
        TCLAP::ValueArg<int> argThreadCount(
          "t", "thread-count", 
          "Number of threads per process. Defaule value is 1", false, 1, "int");
        cmd.add(argThreadCount);
        TCLAP::ValueArg<int> argIterationCount(
          "i", "iteration-count", "Number of iterations.", true, 0, "int");
        cmd.add(argIterationCount);

        TCLAP::ValueArg<std::string> argOutputFilePath(
          "o", "output-file-path", 
          "Output HDF5 file path. The reconstructed 3D image is stored "
          "in \"/image\" dataset. Default file is ./output.h5.", 
          false, "./output.h5", "string");
        cmd.add(argOutputFilePath);
        TCLAP::ValueArg<std::string> argInputFilePath(
          "f", "input-file-path", 
          "Input HDF5 file path. This file should include both projection "
          "and theta datasets.", 
          false, "", "string");
        cmd.add(argInputFilePath);

        std::vector<std::string> allowedAlgs{"sirt", "mlem", "pml", "apmlr"};
        TCLAP::ValuesConstraint<std::string> allowedVals(allowedAlgs);
        TCLAP::ValueArg<std::string> argReconstructionAlg(
          "r", "reconstruction-algorithm", 
          "Reconstruction algorithm to be applied.", 
          true, "", &allowedVals);
        cmd.add(argReconstructionAlg);

        cmd.parse(argc, argv);

        bool mxor = argProjectionFilePath.isSet() | argThetaFilePath.isSet();
        if(!(mxor xor argInputFilePath.isSet())){
          std::cerr << "Both input and projection (or theta) files are set "
            "(or unset). Only one of them should be given." 
            << std::endl;
          std::exit(0);
        }
        if(mxor){
          if(!(argProjectionFilePath.isSet() and argThetaFilePath.isSet())){
            std::cerr << "Both projection and theta files need to be set." 
              << std::endl;
            std::exit(0);
          }
          kProjectionFilePath = argProjectionFilePath.getValue();
          kThetaFilePath = argThetaFilePath.getValue();
        } else {
          kProjectionFilePath = argInputFilePath.getValue();
          kThetaFilePath = argInputFilePath.getValue();
        }

        kReconstructionAlg = argReconstructionAlg.getValue();
        kInputFilePath = argInputFilePath.getValue();
        kOutputFilePath = argOutputFilePath.getValue();
        kProjectionDatasetPath = argProjectionDatasetPath.getValue();
        kThetaDatasetPath = argThetaDatasetPath.getValue();
        kOutputDatasetPath = argOutputDatasetPath.getValue();
        iteration_count = argIterationCount.getValue();
        center = argCenter.getValue();
        thread_count = argThreadCount.getValue();
        b0 = argBetaZero.getValue();
        b1 = argBetaOne.getValue();
        d0 = argDeltaZero.getValue();
        d1 = argDeltaOne.getValue();
        regw = argRegw.getValue();
        degree_to_radian = argDegreeToRadian.getValue();

        std::cout << "Reconstruction algorithm=" << kReconstructionAlg << 
          std::endl;
        std::cout << "Projection file path=" << kProjectionFilePath << 
          std::endl;
        std::cout << "Projection dataset path in hdf5=" << 
          kProjectionDatasetPath << std::endl;
        std::cout << "Theta file path=" << kThetaFilePath << std::endl;
        std::cout << "Theta dataset path=" << kThetaDatasetPath << std::endl;
        std::cout << "Output file path=" << kOutputFilePath << std::endl;
        std::cout << "Recon. dataset path=" << kOutputDatasetPath << std::endl;
        std::cout << "Number of iterations=" << iteration_count << std::endl;
        std::cout << "Center value=" << center << std::endl;
        std::cout << "Number of threads per process=" << thread_count << 
          std::endl;
        std::cout << "Degree to radian=" << degree_to_radian  << 
          std::endl;
        if(argReconstructionAlg.getValue()=="pml" ||
            argReconstructionAlg.getValue()=="apmlr"){
          std::cout << "Beta 0=" << b0 << std::endl;
        }
        if(argReconstructionAlg.getValue()=="apmlr"){
          std::cout << "Beta 1=" << b1 << std::endl;
          std::cout << "Delta 0=" << d0 << std::endl;
          std::cout << "Delta 1=" << d1 << std::endl;
          std::cout << "Regw=" << regw << std::endl;
        }
      }
      catch (TCLAP::ArgException &e)
      {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << 
          std::endl;
      }

      /* Initiate middleware's communication layer */
      comm = std::make_shared<DISPCommMPI<float>>(&argc, &argv);

    }
};

struct TraceData{
  public:
    ADataRegion<float> *sinograms_=nullptr;
    TraceMetadata *metadata_=nullptr;

    TraceMetadata& metadata() const { return *metadata_; };
    void metadata(TraceMetadata *metadata__) { metadata_ = metadata__; };
    ADataRegion<float>& sinograms() const { return *sinograms_; };
    void sinograms(ADataRegion<float> *sinograms__) { sinograms_ = sinograms__; };

    ~TraceData(){
      delete sinograms_;
      delete metadata_;
    }
};

class TracePH5IO {
  private:
    DISPCommBase<float> const &comm;

    std::string const &inputFilePath;
    std::string const &outputFilePath;
    std::string const &outputDatasetPath;

    std::string const &projectionFilePath;
    std::string const &projectionDatasetPath;
    std::string const &thetaFilePath;
    std::string const &thetaDatasetPath;

    float center;

    bool const degree_to_radian;

    struct TraceData &trace_data;

  public:
    TracePH5IO(
      DISPCommBase<float> &kcomm,
      std::string &kInputFilePath,
      std::string &kOutputFilePath,
      std::string &kOutputDatasetPath,
      std::string &kProjectionFilePath,
      std::string &kProjectionDatasetPath,
      std::string &kThetaFilePath,
      std::string &kThetaDatasetPath,
      float kCenter,
      bool kdegree_to_radian,
      struct TraceData &kTraceData):
        comm {kcomm},
        inputFilePath {kInputFilePath},
        outputFilePath {kOutputFilePath},
        outputDatasetPath {kOutputDatasetPath},
        projectionFilePath {kProjectionFilePath},
        projectionDatasetPath {kProjectionDatasetPath},
        thetaFilePath {kThetaFilePath},
        thetaDatasetPath {kThetaDatasetPath},
        center {kCenter},
        degree_to_radian {kdegree_to_radian},
        trace_data {kTraceData}
    {}

    TraceData& read(){
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
};

class TraceEngine {
  private:
    TraceData &trace_data;
    int thread_count;
    int iteration;

    float recon_space_init_val;

  public:
    TraceEngine(
      std::string recon_alg, 
      TraceData &trace_data,
      DISPCommBase<float> &comm,
      int thread_count,
      int iteration
    ):
      trace_data {trace_data},
      thread_count {thread_count},
      iteration {iteration}
    {
      if(recon_alg=="sirt"){
        recon_space_init_val=0.;
        auto main_recon_space = new SIRTReconSpace(
            trace_data.metadata().num_slices(), 
            2*trace_data.metadata().num_cols()*trace_data.metadata().num_cols());
        main_recon_space->Initialize(trace_data.metadata().num_grids());
        main_recon_space->reduction_objects().ResetAllItems(recon_space_init_val);

        /* Prepare processing engine and main reduction space for other threads */
        auto engine =
          new DISPEngineReduction<SIRTReconSpace, float>(
              &comm,
              main_recon_space,
              thread_count);
              /// thread_count=0 for automatically assign the number of threads
      }
      else if(recon_alg=="mlem"){
        std::cerr << "Algorithm is not ready: " << recon_alg << std::endl;
        exit(0);
      }
      else if(recon_alg=="pml"){
        std::cerr << "Algorithm is not ready: " << recon_alg << std::endl;
        exit(0);
      }
      else if(recon_alg=="apmlr"){
        std::cerr << "Algorithm is not ready: " << recon_alg << std::endl;
        exit(0);
      }
      else{
        std::cerr << "Unknown algorithm: " << recon_alg << std::endl;
        exit(0);
      }
    }

    void process(TraceData &trace_data, int thread_count, int iteration){

    }

};

int main(int argc, char **argv)
{
  TraceRuntimeConfig config(argc, argv);

  //TracePH5IO trace_io(
  //      config.comm, 
  //      config.kProjectionFilePath, config.kProjectionDatasetPath,
  //      config.kThetaFilePath, config.kThetaDatasetPath,
  //      config.kOutputFilePath, config.kOutputDatasetPath,
  //      config.degree_to_radian);

  //TraceData &trace_data = trace_io.read();

  //TraceEngine trace_engine(config.kReconstructionAlg, trace_data.metadata());
  //engine.process(
  //  trace_data, 
  //  config.thread_count, config.iteration);

  //trace_io.write(recon_engine.image());
}

