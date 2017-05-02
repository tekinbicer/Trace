template <typename RST, typename DT>
void DISPEngineBase<RST,DT>::InitReductionSpaces(int num_replicas, AReductionSpaceBase<RST,DT> *conf_space)
{
  reduction_spaces_.push_back(conf_space);
  for(int i=1; i<num_replicas; i++)
    reduction_spaces_.push_back(conf_space->Clone());
}

template <typename RST, typename DT>
void DISPEngineBase<RST,DT>::DeleteReductionSpaces()
{
  for(auto &obj : reduction_spaces_)
    delete obj;
}


template <typename RST, typename DT>
int DISPEngineBase<RST,DT>::NumProcessors(){
  int num_procs = 0;

  // With C++11 
  num_procs = std::thread::hardware_concurrency();
  // Below code should works with Linux, Solaris, AIX and Mac OS X(>= 10.4)
  //num_procs = sysconf(_SC_NPROCESSORS_ONLN);

  // If unable to compute hardware threads, set default to 1
  if(num_procs<1) num_procs = 1;
    // throw std::length_error("Number of available processors is <1!");

  return num_procs;
}

template <typename RST, typename DT>
DISPEngineBase<RST,DT>::DISPEngineBase(
    DISPCommBase<DT> *comm,
    AReductionSpaceBase<RST, DT> *conf_reduction_space_i, 
    int num_reduction_threads):
  replication_type_(FULL_REPLICATION)
{
  num_procs_ = NumProcessors();
  num_reduction_threads_ = 
    (num_reduction_threads<1) ? num_procs_ : num_reduction_threads;

  comm_ = comm;

  InitReductionSpaces(num_reduction_threads_, conf_reduction_space_i);
}

template <typename RST, typename DT>
DISPEngineBase<RST,DT>::~DISPEngineBase() {
  DeleteReductionSpaces();
}

template <typename RST, typename DT>
void DISPEngineBase<RST,DT>::Print(){
  DT total=0;
  std::cout << "DISPEngineBase::Print function is called.." << std::endl;
  for(auto &reduction_space : reduction_spaces_){
    std::cout << "reduction_space: " << 
      reduction_space->reduction_objects()[0][0] << std::endl;
    total+=reduction_space->reduction_objects()[0][0];
  }

  std::cout << "Total=" << total << std::endl;
}

template <typename RST, typename DT>
int DISPEngineBase<RST,DT>::num_procs() const {return num_procs_;};
template <typename RST, typename DT>
int DISPEngineBase<RST,DT>::num_reduction_threads() const {return num_reduction_threads_;};
