#if defined ILLUM && defined SOLVERS && defined USE_MPI
#include "illum_mpi_struct.h"
#include "mpi.h"

MPI_Datatype MPI_RECV_ILLUM_T;
MPI_Datatype MPI_SUB_ARRAY_RECV_T;

void illum::MpiInitStruct(const grid_directions_t &grid) {
  // MPI_Type_create_resized(MPI_DOUBLE, 0, sizeof(double) * grid.size, &MPI_RECV_ILLUM_T);
  int array_size[] = {static_cast<int>(grid.size)};
  int array_subsize[] = {static_cast<int>(1)};
  int array_start[] = {0};

  /* Create a subarray datatype */
  MPI_Type_create_subarray(1, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_DOUBLE, &MPI_RECV_ILLUM_T);
  MPI_Type_commit(&MPI_RECV_ILLUM_T);
}

void illum::spectrum::MpiInitStruct(const grid_t &grid) {

  int array_size[] = {static_cast<int>(grid.size_frq * grid.size_dir)};
  int array_subsize[] = {static_cast<int>(grid.size_frq)};
  int array_start[] = {0};

  /* Create a subarray datatype */
  MPI_Type_create_subarray(1, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_DOUBLE, &MPI_SUB_ARRAY_RECV_T);
  MPI_Type_commit(&MPI_SUB_ARRAY_RECV_T);

  // MPI_Send(Illum_local.data(), N*F,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  // MPI_Recv(Illum.data(), N,subarray2,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

#endif