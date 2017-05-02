template <typename T, typename I>
MirroredRegionBase<T,I>::~MirroredRegionBase() {}

template <typename T, typename I>
MirroredRegionBase<T,I>::MirroredRegionBase(MirroredRegionBase<T,I> &region)
      : MirroredRegionBareBase<T>(region),
        metadata_{region.metadata_}
{
  if(metadata_ == nullptr)
    throw std::invalid_argument("Metadata pointer cannot be null!");
}

template <typename T, typename I>
MirroredRegionBase<T,I>::MirroredRegionBase(
        ADataRegion<T> const * const parent, 
        T * const data, size_t count, size_t index, 
        I * const metadata)
      : MirroredRegionBareBase<T>(parent, data, count, index),
        metadata_{metadata}
{
  if(metadata == nullptr)
    throw std::invalid_argument("Metadata pointer cannot be null!");
}

template <typename T, typename I>
I& MirroredRegionBase<T,I>::metadata() const { return *metadata_; }

template <typename T, typename I>
MirroredRegionBase<T,I>* MirroredRegionBase<T,I>::Clone()
{
  MirroredRegionBase<T, I> *mirror = 
    new MirroredRegionBase<T, I>(*this);
  return mirror;
}
