template <typename T>
MirroredRegionBareBase<T>::~MirroredRegionBareBase() 
{}

template <typename T>
MirroredRegionBareBase<T>::MirroredRegionBareBase(MirroredRegionBareBase<T> &region)
  : MirroredRegionBareBase<T>(
      region.parent_, 
      region.data_, 
      region.count_, 
      region.index_)
{}

template <typename T>
MirroredRegionBareBase<T>::MirroredRegionBareBase(
        ADataRegion<T> const * const parent, 
        T * const data, size_t count, size_t index)
  : data_{data},
    count_{count},
    index_{index},
    parent_{parent}
{
  if(data == nullptr)
    throw std::invalid_argument("Data pointer cannot be null!");
  if(count<1)
    throw std::invalid_argument("Number of data items cannot be less than 1!");
  if(parent == nullptr)
    throw std::invalid_argument("Parent pointer cannot be null!");
  if(data+count > &((*parent)[parent->count()]))
    throw std::out_of_range("data+count is out of range considering parent data ptr!");
}

template <typename T>
T& MirroredRegionBareBase<T>::operator[](const size_t index) const { return data_[index]; }

template <typename T>
size_t MirroredRegionBareBase<T>::count() const { return count_; }
template <typename T>
size_t MirroredRegionBareBase<T>::index() const { return index_; }

template <typename T>
MirroredRegionBareBase<T>* MirroredRegionBareBase<T>::Clone()
{
  MirroredRegionBareBase<T> *mirror = 
    new MirroredRegionBareBase<T>(parent_, data_, count_, index_);
  return mirror;
}
