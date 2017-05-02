#ifndef DISP_SRC_DISP_ENGINE_REDUCTION_H_
#define DISP_SRC_DISP_ENGINE_REDUCTION_H_

#include "disp_engine_base.h"
#include "mirrored_region_bare_base.h"
#include "reduction_space_a.h"
#include <deque>

template <typename RST, typename DT>
class DISPEngineReduction : public DISPEngineBase<RST, DT>{
  protected:
    std::mutex work_queue_mutex;

    virtual MirroredRegionBareBase<DT>* Partitioner(ADataRegion<DT> &input_data,
        int req_units);

    virtual MirroredRegionBareBase<DT>* PartitionWrapper(
        ADataRegion<DT> &input_data, int req_units);

    /**
     * \brief Wrapper for reduction thread.
     *
     * \param[in] reduction_space This thread's reduction space.
     * \param[in] input_data The whole input data.
     * \param[in] req_units The number of items that will be assigned after each
     * request. this thread will have from
     */
    virtual void ReductionWrapper(AReductionSpaceBase<RST, DT> &reduction_space,
                                  ADataRegion<DT> &input_data, int &req_units);

    virtual void SeqInPlaceLocalSynch(
        std::vector<AReductionSpaceBase<RST, DT> *> &reduction_spaces);

    virtual void ResetAllReductionObjects(AReductionSpaceBase<RST, DT> 
                                          &reduction_space, DT &val);


  public:
    virtual void GlobalInPlaceSynch(
        DataRegion2DBareBase<DT> &dr, 
        DISPCommBase<DT> &comm);

    virtual void DistInPlaceGlobalSynchWrapper();

    virtual void SeqInPlaceLocalSynchWrapper();

    void ParInPlaceLocalSynchHelper(
        std::deque<std::vector<AReductionSpaceBase<RST, DT> *>*> &work_queue);

    void PartitionReductionSpaces(
        std::vector<AReductionSpaceBase<RST, DT> *> &input_spaces,
        int partition_size,
        std::deque<std::vector<AReductionSpaceBase<RST, DT> *>*> &work_queue,
        std::vector<AReductionSpaceBase<RST, DT> *> &target_spaces);

    virtual void ParInPlaceLocalSynch(
        std::vector<AReductionSpaceBase<RST, DT> *> &reduction_spaces,
        int partition_size,
        int num_threads);

    virtual void ParInPlaceLocalSynchWrapper();

    virtual void RunParallelReduction(ADataRegion<DT> &input_data, int req_units);


    virtual void ResetReductionSpaces(DT &val);

    DISPEngineReduction(
        DISPCommBase<DT> *comm,
        AReductionSpaceBase<RST, DT> *conf_reduction_space_i,
        int num_reduction_threads);
};

#include "disp_engine_reduction.inl"

#endif    // DISP_SRC_DISP_ENGINE_REDUCTION_H_
