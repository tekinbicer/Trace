/** Clone */
template <typename T>
ADataRegion<T>* DataRegionBareBase<T>::Clone(){
  return new DataRegionBareBase<T>(*this);
}

/** Constructors */
template <typename T>
DataRegionBareBase<T>::DataRegionBareBase(const size_t count)
: ADataRegion<T>(count)
{}
template <typename T>
DataRegionBareBase<T>::DataRegionBareBase(T * const data, const size_t count)
: ADataRegion<T>(data, count)
{}
template <typename T>
DataRegionBareBase<T>::DataRegionBareBase(const ADataRegion<T> &region)
: ADataRegion<T>(region)
{}
template <typename T>
DataRegionBareBase<T>::DataRegionBareBase(const ADataRegion<T> &&region)
: ADataRegion<T>(std::move(region))
{}

/** Assignments */
template <typename T>
DataRegionBareBase<T>& DataRegionBareBase<T>::operator=(const ADataRegion<T> &region)
{
  ADataRegion<T>::operator=(region);
}
template <typename T>
DataRegionBareBase<T>& DataRegionBareBase<T>::operator=(const ADataRegion<T> &&region)
{
  ADataRegion<T>::operator=(std::move(region));
}

/** MirroredRegion Functions */
template <typename T>
MirroredRegionBareBase<T>* DataRegionBareBase<T>::NextMirroredRegion(const size_t count)
{
  if(ADataRegion<T>::index_ >= ADataRegion<T>::count()) 
    return nullptr;

  size_t num_units = 
    (ADataRegion<T>::index_+count > ADataRegion<T>::count()) ? 
    (ADataRegion<T>::count() - ADataRegion<T>::index_) : count;

  auto region = MirrorRegion(ADataRegion<T>::index_, num_units);
  ADataRegion<T>::index_ += num_units;

  return region;
}

template <typename T>
MirroredRegionBareBase<T>* DataRegionBareBase<T>::MirrorRegion(const size_t index, const size_t count)
{
  if(index+count > ADataRegion<T>::count()) 
    throw std::out_of_range("Mirrored region is out of range!");

  MirroredRegionBareBase<T> *mirrored_region = 
    new MirroredRegionBareBase<T>(this, &ADataRegion<T>::operator[](0)+index, count, index);
  ADataRegion<T>::mirrored_regions_.push_back(mirrored_region);

  return mirrored_region;
}
