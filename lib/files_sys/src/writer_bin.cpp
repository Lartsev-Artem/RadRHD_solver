#include "global_def.h"

#include "writer_bin.h"

#include "mpi_ext.h"
#include "solvers_config.h"
#include "solvers_types.h"

/*! \addtogroup file_sys Файловый модуль
    @{
*/

namespace files_sys {
namespace bin {

int WriteNormals(const std::string &name_file_normals, std::vector<Normals> &normals) {
  FILE *f;
  OPEN_FILE(f, name_file_normals.c_str(), "wb");

  int n = normals.size();
  fwrite(&n, sizeof(int), 1, f);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < CELL_SIZE; j++) {
      fwrite(&normals[i].n[j], sizeof(Vector3), 1, f);
    }
  }
  fclose(f);

  return e_completion_success;
}

// поэлементная запись данных в файл
#define WRITE_FILE_ELEM(name_file, data, value)  \
  {                                              \
    FILE *f;                                     \
    OPEN_FILE(f, name_file, "wb");               \
    int n = data.size();                         \
    fwrite(&n, sizeof(int), 1, f);               \
    for (auto &el : data) {                      \
      fwrite(&el.value, sizeof(el.value), 1, f); \
    }                                            \
    fclose(f);                                   \
  }

/// \todo CHECK THIS!!!
#ifdef RHLLC_MPI
#error "todo module"
#include "../solve_module/MPI_utils.h"
static inline size_t WriteFileSolutionMPI(const std::string &main_dir, const grid_t &grid, const bound_size_t &size) {
  // int np, myid; MPI_GET_INF(np, myid);
  int myid = claster_cfg.id_hllc;
  int left = size.left;
  int right = size.right;

#ifdef ILLUM
  if (claster_cfg.comm_illum != MPI_COMM_NULL) {
    // todo: это на узле с излучением!
    std::vector<Type> illum;
    GetDirectionIllumFromFace(grid.size, 0, grid.Illum, illum);
    WriteSimpleFileBin(main_dir + F_ILLUM, illum);

    // int size_illum = illum.size();
    // WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_ILLUM).c_str(), (illum.data() + left), size_illum, left, right);
    return 0;
  }

#if !defined USE_CUDA

  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_ENERGY).c_str(), grid.cells, illum_val.energy, left, right);
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_STREAM).c_str(), grid.cells, illum_val.stream, left, right);
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_IMPULS).c_str(), grid.cells, illum_val.impuls, left, right);
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_DIVSTREAM).c_str(), grid.cells, illum_val.div_stream, left, right);
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_DIVIMPULS).c_str(), grid.cells, illum_val.div_impuls, left, right);
#else
#pragma warning "shift/sizeof(double) or just shift?"
  WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_ENERGY).c_str(), grid.energy, grid.size, left, right);
  WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_STREAM).c_str(), grid.stream, grid.size, left, right);
  WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_IMPULS).c_str(), grid.impuls, grid.size, left, right);
  WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_DIVSTREAM).c_str(), grid.divstream, grid.size, left, right);
  WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_DIVIMPULS).c_str(), grid.divimpuls, grid.size, left, right);
#endif
#endif

#if (defined HLLC || defined RHLLC)
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_DENSITY).c_str(), grid.cells, phys_val.d, left, right);
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_PRESSURE).c_str(), grid.cells, phys_val.p, left, right);
  WRITE_FILE_VECTOR_MPI(claster_cfg.comm_hllc, (main_dir + F_VELOCITY).c_str(), grid.cells, phys_val.v, left, right);
#endif

  return 0;
}

static inline size_t WriteFileSolutionSplit(const std::string &main_dir, const grid_t &grid, const bound_size_t &size) {
  // int np, myid; MPI_GET_INF(np, myid);
  int myid = claster_cfg.id_hllc;
  int left = size.left;
  int right = size.right;

#ifdef ILLUM
  if (MPI_IS_ILLUM) {
    std::vector<Type> illum;
    GetDirectionIllumFromFace(grid.size, 0, grid.Illum, illum);
    WriteSimpleFileBin(main_dir + F_ILLUM, illum);

#if !defined USE_CUDA

    WriteFileVectorMPI(MPI_COMM_ILLUM, main_dir + F_ENERGY, grid.cells, offsetof(illum_value_t, energy), sizeof(Type), left, right);
    WriteFileVectorMPI(MPI_COMM_ILLUM, main_dir + F_STREAM, grid.cells, offsetof(illum_value_t, stream), sizeof(Vector3), left, right);
    WriteFileVectorMPI(MPI_COMM_ILLUM, main_dir + F_IMPULS, grid.cells, offsetof(illum_value_t, impuls), sizeof(Matrix3), left, right);
    WriteFileVectorMPI(MPI_COMM_ILLUM, main_dir + F_DIVSTREAM, grid.cells, offsetof(illum_value_t, div_stream), sizeof(Type), left, right);
    WriteFileVectorMPI(MPI_COMM_ILLUM, main_dir + F_DIVIMPULS, grid.cells, offsetof(illum_value_t, div_impuls), sizeof(Vector3), left, right);
#else
#pragma error "new config like no use_cuda"
#pragma warning "shift/sizeof(double) or just shift?"
    WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_ENERGY).c_str(), grid.energy, grid.size, left, right);
    WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_STREAM).c_str(), grid.stream, grid.size, left, right);
    WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_IMPULS).c_str(), grid.impuls, grid.size, left, right);
    WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_DIVSTREAM).c_str(), grid.divstream, grid.size, left, right);
    WRITE_FILE_MPI(claster_cfg.comm_hllc, (main_dir + F_DIVIMPULS).c_str(), grid.divimpuls, grid.size, left, right);
#endif
#endif
  }

#if (defined HLLC || defined RHLLC)

  if (MPI_IS_HLLC) {
    WriteFileVectorMPI(MPI_COMM_HLLC, main_dir + F_DENSITY, grid.cells, offsetof(flux_t, d), sizeof(Type), left, right);
    WriteFileVectorMPI(MPI_COMM_HLLC, main_dir + F_PRESSURE, grid.cells, offsetof(flux_t, p), sizeof(Type), left, right);
    WriteFileVectorMPI(MPI_COMM_HLLC, main_dir + F_VELOCITY, grid.cells, offsetof(flux_t, v), sizeof(Vector3), left, right);
  }
#endif

  return 0;
}

#else

static inline size_t WriteFileSolutionOrder(const std::string &main_dir, const grid_t &grid) {

#ifdef ILLUM
#include "../solve_module/illum/illum_utils.h"
  std::vector<Type> illum;
  GetDirectionIllumFromFace(grid.size, 0, grid.Illum, illum);
  WriteSimpleFileBin(main_dir + "Illum.bin", illum);

#if !defined USE_CUDA
  WRITE_FILE_ELEM((main_dir + "energy.bin").c_str(), grid.cells, illum_val.energy);

  WRITE_FILE_ELEM((main_dir + "stream.bin").c_str(), grid.cells, illum_val.stream);

  WRITE_FILE_ELEM((main_dir + "impuls.bin").c_str(), grid.cells, illum_val.impuls);

  WRITE_FILE_ELEM((main_dir + "divstream.bin").c_str(), grid.cells, illum_val.div_stream);

  WRITE_FILE_ELEM((main_dir + "divimpuls.bin").c_str(), grid.cells, illum_val.div_impuls);
#else
  WriteSimple((main_dir + "energy.bin").c_str(), grid.size, grid.energy);

  WriteSimple((main_dir + "stream.bin").c_str(), grid.size, grid.stream);

  WriteSimple((main_dir + "impuls.bin").c_str(), grid.size, grid.impuls);

  WriteSimple((main_dir + "divstream.bin").c_str(), grid.size, grid.divstream);

  WriteSimple((main_dir + "divimpuls.bin").c_str(), grid.size, grid.divimpuls);
#endif
#endif

#if (defined HLLC || defined RHLLC)

  WRITE_FILE_ELEM((main_dir + F_DENSITY).c_str(), grid.cells, phys_val.d);

  WRITE_FILE_ELEM((main_dir + F_PRESSURE).c_str(), grid.cells, phys_val.p);

  WRITE_FILE_ELEM((main_dir + F_VELOCITY).c_str(), grid.cells, phys_val.v);

#endif

  return 0;
}

#endif

size_t WriteFileSolution(const std::string &main_dir, const grid_t &grid) {

#ifdef RHLLC_MPI
  return WriteFileSolutionMPI(main_dir, grid, hllc_loc_size[claster_cfg.id_hllc]);
#else
  if (get_mpi_id() == 0) {
    return WriteFileSolutionOrder(main_dir, grid);
  }
#endif
  return e_completion_success;
}

} // namespace bin
} // namespace files_sys
