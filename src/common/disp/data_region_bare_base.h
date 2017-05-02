/** \file data_region_bare_base.h 
 */
#ifndef DISP_SRC_DISP_DATA_REGION_BARE_BASE_H
#define DISP_SRC_DISP_DATA_REGION_BARE_BASE_H

/**
 * \class DataRegionBareBase
 * \brief Basic implementation of DataRegionA class.
 *
 * Contact: bicer@anl.gov
 */

#include <iostream>
#include <vector>
#include "data_region_a.h"
#include "mirrored_region_bare_base.h"

template <typename T> 
class DataRegionBareBase : public ADataRegion<T> {
  protected:
    virtual MirroredRegionBareBase<T>* MirrorRegion(const size_t index, const size_t count);

  public:
    virtual ADataRegion<T>* Clone();

    // Constructors
    explicit DataRegionBareBase(const size_t count);

    explicit DataRegionBareBase(T * const data, const size_t count);

    DataRegionBareBase(const ADataRegion<T> &region);
    DataRegionBareBase(const ADataRegion<T> &&region);

    // Assignments
    DataRegionBareBase<T>& operator=(const ADataRegion<T> &region);
    DataRegionBareBase<T>& operator=(const ADataRegion<T> &&region);

    virtual MirroredRegionBareBase<T>* NextMirroredRegion(const size_t count);
};

#include "data_region_bare_base.inl"

#endif    // DISP_SRC_DISP_DATA_REGION_BARE_BASE_H
