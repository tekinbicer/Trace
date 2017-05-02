/** \file data_region_base.h 
 */
#ifndef DISP_SRC_DISP_DATA_REGION_BASE_H
#define DISP_SRC_DISP_DATA_REGION_BASE_H

/**
 * \class DataRegionBase
 * \brief Implementation of abstract DataRegion class with metadata data structure.
 *
 * Contact: bicer@anl.gov
 */

#include <iostream>
#include <vector>
#include "data_region_a.h"
#include "mirrored_region_base.h"

template <typename T, typename I> 
class DataRegionBase : public ADataRegion<T> {
  protected:
    I *metadata_ = nullptr;

    virtual MirroredRegionBase<T, I>* MirrorRegion(const size_t index, const size_t count);

  public:

    /** Constructors */
    explicit DataRegionBase(const size_t count, I * const metadata);
    explicit DataRegionBase(T * const data, const size_t count, I * const metadata);
    DataRegionBase(const ADataRegion<T> &region);
    DataRegionBase(ADataRegion<T> &&region);

    /** Assignments */
    DataRegionBase<T, I>& operator=(const ADataRegion<T> &region);
    DataRegionBase<T, I>& operator=(const ADataRegion<T> &&region);

    /** Accessors/Mutators */
    I& metadata() const;

    /** Mirrored region functions */
    virtual MirroredRegionBareBase<T>* NextMirroredRegion(const size_t count);

    virtual ADataRegion<T>* Clone();
};

#include "data_region_base.inl"

#endif    // DISP_SRC_DISP_DATA_REGION_BARE_BASE_H
