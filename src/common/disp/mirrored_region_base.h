/** \file mirrored_region_base.h 
 */
#ifndef DISP_SRC_DISP_MIRRORED_REGION_BASE_H
#define DISP_SRC_DISP_MIRRORED_REGION_BASE_H

/**
 * \class MirroredRegionBase
 * \brief Extended mirrored region class with metadata information.
 *
 * Contact: bicer@anl.gov
 */

#include <vector>
#include "mirrored_region_bare_base.h"


template <typename T, typename I>
class MirroredRegionBase : public MirroredRegionBareBase<T> {
  private:
    I * const metadata_;

  public:
    virtual ~MirroredRegionBase();

    explicit MirroredRegionBase(MirroredRegionBase<T, I> &region);

    explicit MirroredRegionBase(
        ADataRegion<T> const * const parent, 
        T * const data, size_t count, size_t index, 
        I * const metadata);

    I& metadata() const;

    virtual MirroredRegionBase<T, I>* Clone();
};

#include "mirrored_region_base.inl"

#endif    // DISP_SRC_DISP_MIRRORED_REGION_BASE_H
