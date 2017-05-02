/** Constructors */
template <typename T, typename I> 
DataRegionBase<T,I>::DataRegionBase(const size_t count, I * const metadata)
: ADataRegion<T>(count)
  , metadata_{metadata}
{}

template <typename T, typename I> 
DataRegionBase<T,I>::DataRegionBase(T * const data, const size_t count, I * const metadata)
: ADataRegion<T>(data, count)
  , metadata_{metadata}
{}

template <typename T, typename I> 
DataRegionBase<T,I>::DataRegionBase(const ADataRegion<T> &region)
: ADataRegion<T>(region)
{
  DataRegionBase<T, I> *dr = dynamic_cast<DataRegionBase<T, I>*>(&region);
  if(dr != nullptr) metadata_ = dr->metadata_;
}

template <typename T, typename I> 
DataRegionBase<T,I>::DataRegionBase(ADataRegion<T> &&region)
: ADataRegion<T>(std::move(region))
{
  DataRegionBase<T, I> *dr = dynamic_cast<DataRegionBase<T, I>*>(&region);
  if(dr != nullptr){
    metadata_ = dr->metadata_;
    dr->metadata_ = nullptr;
  }
  ADataRegion<T>::operator=(std::move(region));
}


/* Assignments */
template <typename T, typename I> 
DataRegionBase<T, I>& DataRegionBase<T,I>::operator=(const ADataRegion<T> &region)
{
  DataRegionBase<T, I> *dr = dynamic_cast<DataRegionBase<T, I>*>(&region);
  if(dr != nullptr) metadata_ = dr->metadata_;
  ADataRegion<T>::operator=(region);

  return *this;
}
template <typename T, typename I> 
DataRegionBase<T, I>& DataRegionBase<T,I>::operator=(const ADataRegion<T> &&region)
{
  DataRegionBase<T, I> *dr = dynamic_cast<DataRegionBase<T, I>*>(&region);
  if(dr != nullptr){
    metadata_ = dr->metadata_;
    dr->metadata_ = nullptr;
  }
  ADataRegion<T>::operator=(std::move(region));

  return *this;
}


/** Accessors/Mutators */
template <typename T, typename I> 
I& DataRegionBase<T,I>::metadata() const { return *metadata_; }


/** Mirrored region functions */
template <typename T, typename I> 
MirroredRegionBareBase<T>* DataRegionBase<T,I>::NextMirroredRegion(const size_t count)
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

template <typename T, typename I> 
MirroredRegionBase<T, I>* DataRegionBase<T,I>::MirrorRegion(const size_t index, const size_t count)
{
  if(index+count > ADataRegion<T>::count()) 
    throw std::out_of_range("Mirrored region is out of range!");

  MirroredRegionBase<T, I> *mirrored_region = 
    new MirroredRegionBase<T, I>(
        this, 
        &ADataRegion<T>::operator[](0)+index, 
        count, 
        index, 
        metadata_);
  ADataRegion<T>::mirrored_regions_.push_back(mirrored_region);

  return mirrored_region;
}

template <typename T, typename I> 
ADataRegion<T>* DataRegionBase<T,I>::Clone()
{
  ADataRegion<T> *region = new DataRegionBase<T, I>(*this);
  return region;
}
