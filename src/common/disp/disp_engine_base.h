#ifndef DISP_SRC_DISP_ENGINE_BASE_H_
#define DISP_SRC_DISP_ENGINE_BASE_H_

#include <sched.h>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <chrono>
#include "data_region_a.h"
#include "disp_comm_base.h"
#include "reduction_space_a.h"

enum ReplicationTypes{
  FULL_REPLICATION,
  REPLICATED,
  SINGLE
};

template <typename RST, typename DT>
class DISPEngineBase{
  protected:
    std::vector<AReductionSpaceBase<RST, DT> *> reduction_spaces_;

    int num_reduction_threads_;
    int num_procs_;
    std::mutex partitioner_mutex_;
    ReplicationTypes replication_type_;

    DISPCommBase<DT> *comm_;

    virtual void ReductionWrapper(AReductionSpaceBase<RST, DT> &reduction_space,
        ADataRegion<DT> &input_data, int &req_units)=0;

    virtual MirroredRegionBareBase<DT>* PartitionWrapper(ADataRegion<DT> &input_data, 
        int req_units)=0;
    virtual MirroredRegionBareBase<DT>* Partitioner(ADataRegion<DT> &input_data, 
        int req_units)=0;

    virtual void SeqInPlaceLocalSynch(
        std::vector<AReductionSpaceBase<RST, DT>*> &reduction_spaces)=0;

    void DeleteReductionSpaces();
    virtual void ResetAllReductionObjects(
        AReductionSpaceBase<RST, DT> &reduction_space,
        DT &val)=0;

    int NumProcessors();

    void InitReductionSpaces(int num_replicas, 
                             AReductionSpaceBase<RST, DT> *conf_space);

  public:
    DISPEngineBase(
        DISPCommBase<DT> *comm,
        AReductionSpaceBase<RST, DT> *conf_reduction_space, 
        int num_reduction_threads);
    virtual ~DISPEngineBase();

    virtual void RunParallelReduction(ADataRegion<DT> &input_data, 
                                      int req_units)=0;

    virtual void SeqInPlaceLocalSynchWrapper() = 0;
    virtual void ParInPlaceLocalSynchWrapper() = 0;
    virtual void DistInPlaceGlobalSynchWrapper() = 0;
    virtual void ResetReductionSpaces(DT &val) = 0;

    int num_procs() const;
    int num_reduction_threads() const;

    void Print();
};

#include "disp_engine_base.inl"

#endif    // DISP_SRC_DISP_ENGINE_BASE_H_
