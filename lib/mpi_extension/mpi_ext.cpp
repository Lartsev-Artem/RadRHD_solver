#include <cstdint>

#ifdef USE_MPI

#include "mpi.h"

static int8_t id = -1;
static int8_t np = -1;

static void init_value(const MPI_Comm &comm) {
  MPI_Comm_rank(comm, (int *)&id);
  MPI_Comm_size(comm, (int *)&np);
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

MPI_Datatype MPI_flux_t;
MPI_Datatype MPI_flux_illum_elem_t;
MPI_Datatype MPI_hllc_value_t;
MPI_Datatype MPI_flux_all_t;
MPI_Datatype MPI_flux_elem_t;

#else
int8_t get_mpi_id() { return 0; }
int8_t get_mpi_np() { return 1; }
#endif // USE_MPI