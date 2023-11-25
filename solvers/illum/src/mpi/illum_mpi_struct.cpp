#include "illum_mpi_struct.h"
#include "mpi.h"

MPI_Datatype MPI_RECV_ILLUM_T;

void illum::MpiInitStruct(const grid_directions_t &grid) {
  MPI_Type_create_resized(MPI_DOUBLE, 0, sizeof(double) * grid.size, &MPI_RECV_ILLUM_T);
}