#include "trace_comm.h"

void trace_comm::GlobalNeighborUpdate(
    ADataRegion<float> &recon_a,
    int num_slices,
    size_t slice_size)  // Reconstruction object
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  int mpi_size; 
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  float *recon = &recon_a[0];

  float *send_ptr, *recv_ptr;

  MPI_Status statuses[4];
  MPI_Request requests[4];


  /* No neighboring rank */
  if(mpi_size==1) return;

  if(mpi_rank == 0){
    /* With bottom */
    send_ptr = ((float*)recon) + slice_size*(num_slices-2);
    recv_ptr = ((float*)recon) + slice_size*(num_slices-1);

    MPI_Isend(send_ptr, slice_size, MPI_FLOAT, mpi_rank+1, 0,
      MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recv_ptr, slice_size, MPI_FLOAT, mpi_rank+1, 0,
      MPI_COMM_WORLD, &requests[1]);

    MPI_Waitall(2, requests, statuses);

  } else if (0<mpi_rank && mpi_rank<mpi_size-1){
    /* With up */
    send_ptr = ((float*)recon) + slice_size;
    recv_ptr = ((float*)recon);
    MPI_Isend(send_ptr, slice_size, MPI_FLOAT, mpi_rank-1, 0, 
        MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recv_ptr, slice_size, MPI_FLOAT, mpi_rank-1, 0,
        MPI_COMM_WORLD, &requests[1]);

    /* With bottom */
    float *send_ptr_b = ((float*)recon) + slice_size*(num_slices-2);
    float *recv_ptr_b = ((float*)recon) + slice_size*(num_slices-1);
    MPI_Isend(send_ptr_b, slice_size, MPI_FLOAT, mpi_rank+1, 0,
      MPI_COMM_WORLD, &requests[2]);
    MPI_Irecv(recv_ptr_b, slice_size, MPI_FLOAT, mpi_rank+1, 0,
      MPI_COMM_WORLD, &requests[3]);

    MPI_Waitall(4, requests, statuses);

  } else{ /* mpi_rank == mpi_size-1 */
    /* With up */
    send_ptr = ((float*)recon) + slice_size;
    recv_ptr = ((float*)recon);
    MPI_Isend(send_ptr, slice_size, MPI_FLOAT, mpi_rank-1, 0, 
        MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recv_ptr, slice_size, MPI_FLOAT, mpi_rank-1, 0,
        MPI_COMM_WORLD, &requests[1]);

    MPI_Waitall(2, requests, statuses);
  }
}

void trace_comm::SetupSliceCommGroups(
  int slice_id, 
  MPI_Comm &newcomm)
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_split(MPI_COMM_WORLD, slice_id, mpi_rank, &newcomm);
}


void trace_comm::UpdateReconReplicasCommGroups(
  DataRegion2DBareBase<float> &dr,
  MPI_Comm &group_comm)
{
  for(size_t i=0; i<dr.num_rows(); i++)
    MPI_Allreduce(MPI_IN_PLACE, &dr[i][0], dr.num_cols(), MPI_FLOAT,
        MPI_SUM, group_comm);
}


void trace_comm::UpdateSliceCommGroups(
  ADataRegion<float> &recon_a,
  size_t slice_size,
  MPI_Comm &slice_comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &recon_a[0], slice_size, MPI_FLOAT, MPI_SUM, slice_comm);
}


