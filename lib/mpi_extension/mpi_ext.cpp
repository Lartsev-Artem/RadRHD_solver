#include <cstdint>

#ifdef USE_MPI

#include "dbgdef.h"
#include "mpi.h"

static int8_t id = -1;
static int8_t np = -1;

static void init_value(const MPI_Comm &comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  id = (int8_t)rank;
  np = (int8_t)size;

  WRITE_LOG("MPI Claster config: %d %d\n", id, np);
}

int8_t get_mpi_id() {
  if (id < 0) {
    init_value(MPI_COMM_WORLD);
  }
  return id;
}

int8_t get_mpi_np() {
  if (np < 0) {
    init_value(MPI_COMM_WORLD);
  }
  return np;
}

int8_t get_mpi_id(const MPI_Comm &comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  return (int8_t)rank;
}

int8_t get_mpi_np(const MPI_Comm &comm) {
  int size;
  MPI_Comm_size(comm, &size);
  return (int8_t)size;
}

MPI_Datatype MPI_flux_t;
MPI_Datatype MPI_flux_illum_elem_t;
MPI_Datatype MPI_hllc_value_t;
MPI_Datatype MPI_flux_all_t;
MPI_Datatype MPI_flux_elem_t;

#else
int8_t get_mpi_id() { return 0; }
int8_t get_mpi_np() { return 1; }
#endif // USE_MPI
