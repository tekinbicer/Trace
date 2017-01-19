#ifndef TRACE_COMMON_TRACERUNTIMECONFIG_H
#define TRACE_COMMON_TRACERUNTIMECONFIG_H

class TraceRuntimeConfig {
  private:

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


    DISPCommBase<float>& comm();

    TraceRuntimeConfig(int argc, char **argv);
};

#endif /// TRACE_COMMON_TRACERUNTIMECONFIG_H
