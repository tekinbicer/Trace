/** Constructors & Assignments */
template <typename T>
ADataRegion<T>::ADataRegion(const size_t count)
  : data_{new T[count]}
  , count_{count}
  , index_{0}
{}

template <typename T>
ADataRegion<T>::ADataRegion(T * const data, const size_t count)
  : data_{data}
  , count_{count}
  , index_{0}
{}

template <typename T>
ADataRegion<T>::ADataRegion(const ADataRegion<T> &region)
  : data_{new T[region.count_]}
  , count_{region.count_}
  , index_{0}
{
  std::copy(region.data_, region.data_ + region.count_, data_);
}

template <typename T>
ADataRegion<T>::ADataRegion (ADataRegion<T> &&region)
  : data_{nullptr}
  , count_{0}
  , index_{0}
{
  *this = std::move(region);
  std::cout << "ADataRegion:count=" << count_ << std::endl;
}

template <typename T>
ADataRegion<T>& ADataRegion<T>::operator=(const ADataRegion<T> &region)
{
  if(this != &region) {
    if(region.count_ != count_){
      Clear();
      data_ = new T[region.count_];
    }
    std::copy(region.data_, region.data_ + region.count_, data_);

    count_ = region.count_;
  }
  return *this;
}

template <typename T>
ADataRegion<T>& ADataRegion<T>::operator=(ADataRegion<T> &&region)
{
  if(this != &region){
    Clear();
    data_ = region.data_;
    count_ = region.count_;

    mirrored_regions_ = std::move(region.mirrored_regions_);

    region.data_ = nullptr;
    region.count_ = 0;
  }
  return *this;
}
/** End of constructors and assingments */

template <typename T>
ADataRegion<T>::~ADataRegion(){
  Clear();
}

template <typename T>
void ADataRegion<T>::Clear(){
  //std::cout << "ADataRegion: In the clear" << std::endl;
  if(data_ != nullptr){
    delete [] data_;
    data_ = nullptr;
  }
  count_ = 0;
  for(auto &m_region : mirrored_regions_){
    delete m_region;
    m_region = nullptr;
  }
  index_ = 0;
}

template <typename T>
void ADataRegion<T>::ResetMirroredRegionIter(){
  index_ = 0;
  DeleteMirroredRegions();
}

template <typename T>
int ADataRegion<T>::CompareBoundary(ADataRegion<T> &region){
  return (region.count() != count_) ? -1 : 0;
}

template <typename T>
size_t ADataRegion<T>::size() const {
  return count_*sizeof(T);
}

template <typename T>
T& ADataRegion<T>::item (size_t index) const {
  if(index >= count_) 
    throw std::out_of_range("Tried to access out of range index!");
  return data_[index];
}

template <typename T>
void ADataRegion<T>::item(size_t index, T &&value){
  if(index >= count_) 
    throw std::out_of_range("Tried to update out of range index!");
  data_[index] = value;
}

template <typename T>
void ADataRegion<T>::ResetAllItems(T &value){
  for(size_t i=0; i<count_; i++)
    data_[i] = value;
}

template <typename T>
void ADataRegion<T>::DeleteMirroredRegions(){
  for(auto &m_region : mirrored_regions_){
    delete m_region;
    m_region = nullptr;
  }
  mirrored_regions_.clear();
}

template <typename T>
size_t ADataRegion<T>::index() const { return index_; }

template <typename T>
T& ADataRegion<T>::operator[](size_t index) { return data_[index]; }
template <typename T>
const T& ADataRegion<T>::operator[](size_t index) const { return data_[index]; }

template <typename T>
size_t ADataRegion<T>::count() const { return count_; }
template <typename T>
size_t ADataRegion<T>::num_cols() const { return count_; }
template <typename T>
size_t ADataRegion<T>::cols() const { return count_; }
