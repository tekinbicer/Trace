#ifndef TRACE_COMMON_TRACEDATA_H
#define TRACE_COMMON_TRACEDATA_H

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

#endif
