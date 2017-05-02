template <typename T>
DataRegion2DBareBase<T>::DataRegion2DBareBase(size_t rows, size_t cols){
  regions_.reserve(rows);

  for(size_t i=0; i<rows; i++)
    regions_.push_back(new DataRegionBareBase<T>(cols));

  rows_ = rows;
  cols_ = cols;
}

template <typename T>
DataRegion2DBareBase<T>::~DataRegion2DBareBase(){
  for(auto &region : regions_)
    delete region;
  rows_ = 0;
  cols_ = 0;
}

template <typename T>
DataRegion2DBareBase<T>::DataRegion2DBareBase(const DataRegion2DBareBase &dr){
  regions_.reserve(dr.rows());

  for(auto t : dr.regions()){
    auto r = new DataRegionBareBase<T>(*t);
    regions_.push_back(r);
  }

  rows_ = dr.rows();
  cols_ = dr.cols();
}

template <typename T>
DataRegion2DBareBase<T>& DataRegion2DBareBase<T>::operator=(const DataRegion2DBareBase &dr){
  if(this==&dr) return *this;

  for(auto rg : regions_)
    delete rg;

  regions_.clear();
  regions_.reserve(dr.rows());

  for(auto rg : dr.regions())
    regions_.push_back(new DataRegionBareBase<T>(*rg));
  rows_ = dr.rows();
  cols_ = dr.cols();

  return *this;
}

template <typename T>
DataRegionBareBase<T>& DataRegion2DBareBase<T>::operator[](size_t row) const { 
  return *(regions_[row]); 
}

template <typename T>
const std::vector<DataRegionBareBase<T>*>& DataRegion2DBareBase<T>::regions() const { return regions_; };
template <typename T>
size_t DataRegion2DBareBase<T>::num_rows() const { return rows_; }
template <typename T>
size_t DataRegion2DBareBase<T>::num_cols() const { return cols_; }
template <typename T>
size_t DataRegion2DBareBase<T>::rows() const { return rows_; }
template <typename T>
size_t DataRegion2DBareBase<T>::cols() const { return cols_; }
template <typename T>
size_t DataRegion2DBareBase<T>::count() const { return rows_*cols_; }
template <typename T>
size_t DataRegion2DBareBase<T>::size() const {return rows_*cols_*sizeof(T); }


template <typename T>
T& DataRegion2DBareBase<T>::item(size_t row, size_t col) const {
  if(row >= rows_ || col >= cols_)
    throw std::out_of_range("Tried to access out of range offset!");
  return (*regions_[row])[col];
}

template <typename T>
void DataRegion2DBareBase<T>::item(size_t row, size_t col, T &val) const {
  if(row >= rows_ || col >= cols_)
    throw std::out_of_range("Tried to access out of range offset!");
  (*regions_[row])[col] = val;
}

template <typename T>
void DataRegion2DBareBase<T>::ResetAllItems(T &val) const {
  for(size_t i=0; i<rows_; i++)
    for(size_t j=0; j<cols_; j++)
      (*regions_[i])[j] = val;
}

template <typename T>
void DataRegion2DBareBase<T>::ResetAllMirroredRegions(){
  for(size_t i=0; i<rows_; i++)
    (*regions_[i]).ResetMirroredRegionIter();
}

template <typename T>
void DataRegion2DBareBase<T>::copy(DataRegion2DBareBase<T> &dr){
  if(dr.rows() != rows_ || dr.cols() != cols_)
    throw std::out_of_range("DataRegions' ranges do not overlap!");

  dr = *this;
}

