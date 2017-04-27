#include "trace_runtime_config.h"

void TraceRuntimeConfig::Print() {
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
    if(kReconstructionAlg=="pml" ||
        kReconstructionAlg=="apmlr"){
      std::cout << "Beta 0=" << b0 << std::endl;
    }
    if(kReconstructionAlg=="apmlr"){
      std::cout << "Beta 1=" << b1 << std::endl;
      std::cout << "Delta 0=" << d0 << std::endl;
      std::cout << "Delta 1=" << d1 << std::endl;
      std::cout << "Regw=" << regw << std::endl;
    }
}

TraceRuntimeConfig::TraceRuntimeConfig(int argc, char **argv)
{
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
        "", "projection-file-path", 
        "Input projection file path.", false, "", "string");
    cmd.add(argProjectionFilePath);
    TCLAP::ValueArg<std::string> argProjectionDatasetPath(
        "", "projection-dataset-path", "Projection dataset path in HDF5"
        "file. Default value is \"/exchange/data\".", 
        false, "/exchange/data", "string");
    cmd.add(argProjectionDatasetPath);
    TCLAP::ValueArg<std::string> argThetaFilePath(
        "", "theta-file-path", 
        "Input theta file path.", false, "", "string");
    cmd.add(argThetaFilePath);
    TCLAP::ValueArg<std::string> argThetaDatasetPath(
        "", "theta-dataset-path", "Theta dataset path in HDF5 file. "
        "Default value is \"/exchange/theta\".", 
        false, "/exchange/theta", "string");
    cmd.add(argThetaDatasetPath);

    TCLAP::ValueArg<std::string> argOutputDatasetPath(
        "", "output-dataset-path", 
        "Reconstructed 3D image dataset path in hdf5 file.",
        false, "/image", "string");
    cmd.add(argOutputDatasetPath);

    TCLAP::ValueArg<float> argCenter(
        "", "center", "Center value for input dataset.", false, 0., "float");
    cmd.add(argCenter);
    TCLAP::ValueArg<int> argThreadCount(
        "", "thread-count", 
        "Number of threads per process. Defaule value is 1", false, 1, "int");
    cmd.add(argThreadCount);
    TCLAP::ValueArg<int> argIterationCount(
        "", "iteration-count", "Number of iterations.", true, 0, "int");
    cmd.add(argIterationCount);

    TCLAP::ValueArg<std::string> argOutputFilePath(
        "", "output-file-path", 
        "Output HDF5 file path. The reconstructed 3D image is stored "
        "in \"/image\" dataset. Default file is ./output.h5.", 
        false, "./output.h5", "string");
    cmd.add(argOutputFilePath);
    TCLAP::ValueArg<std::string> argInputFilePath(
        "", "input-file-path", 
        "Input HDF5 file path. This file should include both projection "
        "and theta datasets.", 
        false, "", "string");
    cmd.add(argInputFilePath);

    std::vector<std::string> allowedAlgs{"sirt", "mlem", "pml", "apmlr"};
    TCLAP::ValuesConstraint<std::string> allowedVals(allowedAlgs);
    TCLAP::ValueArg<std::string> argReconstructionAlg(
        "", "reconstruction-algorithm", 
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

  }
  catch (TCLAP::ArgException &e)
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << 
      std::endl;
  }
}
