template <typename CT, typename DT>
void AReductionSpaceBase<CT,DT>::Process(MirroredRegionBareBase<DT> &input)
{
  static_cast<CT*>(this)->Reduce(input);
}

// Default operation is sum
template <typename CT, typename DT>
void AReductionSpaceBase<CT,DT>::LocalSynchWith(CT &input_reduction_space)
{
  auto &ri = input_reduction_space.reduction_objects();
  auto &ro = *reduction_objects_;

  if(ri.num_rows()!=ro.num_rows() || ri.num_cols()!=ro.num_cols())
    throw std::range_error("Local and destination reduction objects have different dimension sizes!");

  for(size_t i=0; i<ro.num_rows(); i++)
    for(size_t j=0; j<ro.num_cols(); j++)
      ro[i][j] += ri[i][j];
}


template <typename CT, typename DT>
CT* AReductionSpaceBase<CT,DT>::Clone()
{
  if(reduction_objects_ == nullptr)
    throw std::invalid_argument("reduction objects point to nullptr!");

  auto &red_objs = *reduction_objects_;

  CT *cloned_obj = new CT(red_objs.rows(), red_objs.cols());
  (*cloned_obj).reduction_objects() = red_objs;

  static_cast<CT*>(this)->CopyTo(*cloned_obj);

  return cloned_obj;
}

template <typename CT, typename DT>
DataRegion2DBareBase<DT>& AReductionSpaceBase<CT,DT>::reduction_objects()
{
  return *reduction_objects_; 
}

template <typename CT, typename DT>
DataRegionBareBase<DT>& AReductionSpaceBase<CT,DT>::operator[](int row)
{
  return (*reduction_objects_)[row]; 
}

template <typename CT, typename DT>
AReductionSpaceBase<CT,DT>::AReductionSpaceBase(DataRegion2DBareBase<DT> *reduction_objects) 
  : reduction_objects_(reduction_objects) 
{}

template <typename CT, typename DT>
AReductionSpaceBase<CT,DT>::AReductionSpaceBase(size_t rows, size_t cols)
{
  reduction_objects_ = new DataRegion2DBareBase<DT>(rows, cols);
}

template <typename CT, typename DT>
AReductionSpaceBase<CT,DT>::~AReductionSpaceBase()
{
  delete reduction_objects_; 
}
