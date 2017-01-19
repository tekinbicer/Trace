#ifndef TRACE_COMMON_TRACEPH5IO_H
#define TRACE_COMMON_TRACEPH5IO_H

class TracePH5IO {
  private:
    DISPCommBase<float> const &comm;

    std::string const &outputFilePath;
    std::string const &outputDatasetPath;

    std::string const &projectionFilePath;
    std::string const &projectionDatasetPath;
    std::string const &thetaFilePath;
    std::string const &thetaDatasetPath;

    float center;
    bool const degree_to_radian;

  public:
    TracePH5IO(
      DISPCommBase<float> &kcomm,
      std::string &kOutputFilePath,
      std::string &kOutputDatasetPath,
      std::string &kProjectionFilePath,
      std::string &kProjectionDatasetPath,
      std::string &kThetaFilePath,
      std::string &kThetaDatasetPath,
      float kCenter,
      bool kdegree_to_radian);

    TracePH5IO(TraceRuntimeConfig &config);
    

    TraceData Read();

    void Write();
};

#endif // TRACE_COMMON_TRACEPH5IO_H
