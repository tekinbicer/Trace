/** \file mirrored_region_bare_base.h 
 */
#ifndef DISP_SRC_DISP_MIRRORED_REGION_BARE_BASE_H
#define DISP_SRC_DISP_MIRRORED_REGION_BARE_BASE_H

/**
 * \class MirroredRegionBareBase
 * \brief Basic class that represents the mirrored regions.
 *
 * Contact: bicer@anl.gov
 */

#include <vector>
#include "data_region_a.h"

template <typename T>
class ADataRegion;

template <typename T>
class MirroredRegionBareBase {
  private:
    T * const data_;
    const size_t count_;

    const size_t index_;
    ADataRegion<T> const * const parent_;

  public:
    virtual ~MirroredRegionBareBase();

    explicit MirroredRegionBareBase(MirroredRegionBareBase<T> &region);

    explicit MirroredRegionBareBase(
        ADataRegion<T> const * const parent, 
        T * const data, size_t count, size_t index);

    T& operator[](const size_t index) const;

    size_t count() const;
    size_t index() const;

    virtual MirroredRegionBareBase<T>* Clone();
};

#include "mirrored_region_bare_base.inl"

#endif    // DISP_SRC_DISP_MIRRORED_REGION_BARE_BASE_H
