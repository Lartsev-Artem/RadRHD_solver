#include <cstdint>

#ifdef USE_MPI

#include "dbgdef.h"
#include "mpi.h"

static int8_t id = -1;
static int8_t np = -1;

namespace mpi_private {
static void init_value(const MPI_Comm &comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  id = (int8_t)rank;
  np = (int8_t)size;

  WRITE_LOG("MPI Claster config: %d %d\n", id, np);
}
}; // namespace mpi_private

int8_t get_mpi_id() {
  if (id < 0) {
    mpi_private::init_value(MPI_COMM_WORLD);
  }
  return id;
}

int8_t get_mpi_np() {
  if (np < 0) {
    mpi_private::init_value(MPI_COMM_WORLD);
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

#include "solvers_struct.h"
MPI_Datatype MPI_flux_t;
MPI_Datatype MPI_phys_val_t;
MPI_Datatype MPI_flux_elem_t;

MPI_Datatype MPI_hllc_value_t;

#if 1
void InitMPiStruct() {

  // структура потоков
  {
    int len[3] = {1, 3, 1};
    MPI_Aint pos[3] = {offsetof(flux_t, d), offsetof(flux_t, v), offsetof(flux_t, p)};
    MPI_Datatype typ[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(3, len, pos, typ, &MPI_flux_t);
    MPI_Type_commit(&MPI_flux_t);
  }

  MPI_Type_create_resized(MPI_flux_t, 0, sizeof(elem_t), &MPI_phys_val_t); // структура физических переменных
  WRITE_LOG("sizeof(elem_t)= %lu\n", sizeof(elem_t));

  // перессылка потоков из ячейки
  {
    int len[2] = {1, 1};
    MPI_Aint pos[2] = {offsetof(elem_t, phys_val), offsetof(elem_t, conv_val)};
    MPI_Datatype typ[2] = {MPI_flux_t, MPI_flux_t};
    MPI_Type_create_struct(2, len, pos, typ, &MPI_flux_elem_t);
    MPI_Type_commit(&MPI_flux_elem_t);
    MPI_Type_create_resized(MPI_flux_elem_t, 0, sizeof(elem_t), &MPI_flux_elem_t); // структура физических переменных
  }

  //настройки динамического расчета
  {
    int len[5] = {1, 1, 1, 1, 1};
    MPI_Aint pos[5] = {offsetof(hllc_value_t, T), offsetof(hllc_value_t, CFL), offsetof(hllc_value_t, h_min), offsetof(hllc_value_t, save_timer), offsetof(hllc_value_t, tau)};
    MPI_Datatype typ[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(5, len, pos, typ, &MPI_hllc_value_t);
    MPI_Type_commit(&MPI_hllc_value_t);
  }

  WRITE_LOG("end InitMPiStruct\n");

  return;
}
#else //устаревший формат?
void InitMPiStruct() {

  // структура потоков
  {
    int len[3 + 1] = {1, 3, 1, 1};
    MPI_Aint pos[4] = {offsetof(flux_t, d), offsetof(flux_t, v), offsetof(flux_t, p), sizeof(flux_t)};
    MPI_Datatype typ[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UB};
    MPI_Type_create_struct(4, len, pos, typ, &MPI_flux_t);
    MPI_Type_commit(&MPI_flux_t);
  }

  // структура физических переменных
  {
    int len[1 + 1] = {1, 1};
    MPI_Aint pos[2] = {offsetof(elem_t, phys_val), sizeof(elem_t)};
    MPI_Datatype typ[2] = {MPI_flux_t, MPI_UB};
    MPI_Type_create_struct(2, len, pos, typ, &MPI_phys_val_t);
    MPI_Type_commit(&MPI_phys_val_t);
  }

  // перессылка потоков из ячейки
  {
    int len[2 + 1] = {1, 1, 1};
    MPI_Aint pos[3] = {offsetof(elem_t, phys_val), offsetof(elem_t, conv_val), sizeof(elem_t)};
    MPI_Datatype typ[3] = {MPI_flux_t, MPI_flux_t, MPI_UB};
    MPI_Type_create_struct(3, len, pos, typ, &MPI_flux_elem_t);
    MPI_Type_commit(&MPI_flux_elem_t);
  }

  //настройки динамического расчета
  {
    int len[5 + 1] = {1, 1, 1, 1, 1, 1};
    MPI_Aint pos[6] = {offsetof(hllc_value_t, T), offsetof(hllc_value_t, CFL), offsetof(hllc_value_t, h_min), offsetof(hllc_value_t, save_timer), offsetof(hllc_value_t, tau), sizeof(hllc_value_t)};
    MPI_Datatype typ[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UB};
    MPI_Type_create_struct(6, len, pos, typ, &MPI_hllc_value_t);
    MPI_Type_commit(&MPI_hllc_value_t);
  }

  return;
}
#endif

#else
int8_t get_mpi_id() { return 0; }
int8_t get_mpi_np() { return 1; }
#endif // USE_MPI
