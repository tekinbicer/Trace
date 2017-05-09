#ifndef _REDUCTION_SPACE_BASE_A_
#define _REDUCTION_SPACE_BASE_A_

#include "data_region_2d_bare_base.h"
#include "mirrored_region_bare_base.h"

/// CT: Derived class type
/// DT: Data type on which Reduce function operate
template <typename CT, typename DT>
class AReductionSpaceBase{
  private:
    DataRegion2DBareBase<DT> *reduction_objects_ = nullptr;

  public:
    void Process(MirroredRegionBareBase<DT> &input);

    // Default operation is sum
    virtual void LocalSynchWith(CT &input_reduction_space);


    // Derived class can use this function to perform
    // deep copies
    virtual void CopyTo(CT &target)=0;

    virtual CT *Clone()=0;

    DataRegion2DBareBase<DT>& reduction_objects();

    DataRegionBareBase<DT>& operator[](int row);

    AReductionSpaceBase(DataRegion2DBareBase<DT> *reduction_objects);

    AReductionSpaceBase(size_t rows, size_t cols);

    /// Make sure to call derived class (CT) destructor 
    /// instead of base class (this one) destructor 
    virtual ~AReductionSpaceBase();
};

#include "reduction_space_a.inl"

#endif
