template <typename DT>
void DISPCommMPI<DT>::MPI_AllreduceInPlaceWithType(
    DataRegion2DBareBase<DT> &dr,
    MPI_Datatype input_type,
    MPI_Op op_type,
    MPI_Comm comm
    )
{
  for(size_t i=0; i<dr.num_rows(); i++)
    MPI_Allreduce(MPI_IN_PLACE, &dr[i][0], dr.num_cols(), input_type,
        op_type, comm);
}

template <typename DT>
DISPCommMPI<DT>::DISPCommMPI(int *argc, char ***argv){
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &(this->rank_));
  MPI_Comm_size(MPI_COMM_WORLD, &(this->size_));
}

template <typename DT>
DISPCommMPI<DT>::~DISPCommMPI(){
  MPI_Finalize(); 
}

template <typename DT>
void DISPCommMPI<DT>::Finalize(){}

template <typename DT>
void DISPCommMPI<DT>::GlobalInPlaceCombination(DataRegion2DBareBase<DT> &dr){
  if(std::is_same<float, DT>::value)
    MPI_AllreduceInPlaceWithType(dr, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  else if(std::is_same<double, DT>::value)
    MPI_AllreduceInPlaceWithType(dr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  else if(std::is_same<int, DT>::value)
    MPI_AllreduceInPlaceWithType(dr, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  else if(std::is_same<int64_t, DT>::value)
    MPI_AllreduceInPlaceWithType(dr, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  else 
    throw std::runtime_error("Unknown data type for MPI collective call");

}
