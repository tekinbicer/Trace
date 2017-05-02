/** \file data_region_2d_bare_base.h 
 */
#ifndef DISP_SRC_DISP_DATA_REGION_2D_BARE_BASE_H
#define DISP_SRC_DISP_DATA_REGION_2D_BARE_BASE_H

/**
 * \class DataRegion2DBareBase
 * \brief This class uses DataRegionBareBase in order to provide 
 * replicated output objects.
 *
 * Contact: bicer@anl.gov
 */

#include <iostream>
#include <vector>
#include "data_region_bare_base.h"

template <typename T> 
class DataRegion2DBareBase {
  private:
    std::vector<DataRegionBareBase<T>*> regions_;

    size_t rows_;
    size_t cols_;

  public:
    DataRegion2DBareBase(size_t rows, size_t cols);

    ~DataRegion2DBareBase();

    /* Copy constructor */
    DataRegion2DBareBase(const DataRegion2DBareBase &dr);

    /* Copy assignment */
    DataRegion2DBareBase<T>& operator=(const DataRegion2DBareBase &dr);

    /* TODO: Move constructor & assignment */

    DataRegionBareBase<T>& operator[](size_t row) const; 

    const std::vector<DataRegionBareBase<T>*>& regions() const;
    size_t num_rows() const;
    size_t num_cols() const;
    size_t rows() const;
    size_t cols() const;
    size_t count() const;
    size_t size() const;


    T& item(size_t row, size_t col) const; 

    void item(size_t row, size_t col, T &val) const;

    void ResetAllItems(T &val) const;

    void ResetAllMirroredRegions();
   
    void copy(DataRegion2DBareBase<T> &dr);
};

#include "data_region_2d_bare_base.inl"

#endif    // DISP_SRC_DISP_DATA_REGION_2D_BARE_BASE_H
