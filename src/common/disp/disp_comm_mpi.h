#ifndef DISP_SRC_DISP_COMM_MPI_H
#define DISP_SRC_DISP_COMM_MPI_H

#include <iostream>
#include "disp_comm_base.h"
#include "mpi.h"

template <typename DT>
class DISPCommMPI : public DISPCommBase<DT> {
  private:
    void MPI_AllreduceInPlaceWithType(
        DataRegion2DBareBase<DT> &dr,
        MPI_Datatype input_type,
        MPI_Op op_type,
        MPI_Comm comm);

  public:
    DISPCommMPI(int *argc, char ***argv);
    ~DISPCommMPI();

    void Finalize();

    void GlobalInPlaceCombination(DataRegion2DBareBase<DT> &dr);
};

#include "disp_comm_mpi.inl"

#endif    // DISP_SRC_DISP_COMM_MPI_H
