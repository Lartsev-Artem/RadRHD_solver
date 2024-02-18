#include "global_def.h"
#include "global_value.h"

#include "convert_face_to_cell.h"
#include "mpi_ext.h"
#include "solvers_config.h"
#include "solvers_struct.h"
#include "writer_bin.h"

int files_sys::bin::WriteNormals(const std::string &name_file_normals, std::vector<Normals> &normals) {
  FILE *f;
  OPEN_FILE(f, name_file_normals.c_str(), "wb");

  int n = (int)normals.size();
  fwrite(&n, sizeof(int), 1, f);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < CELL_SIZE; j++) {
      fwrite(&normals[i].n[j], sizeof(Vector3), 1, f);
    }
  }
  fclose(f);

  return e_completion_success;
}

static inline size_t files_sys::bin::WriteFileSolutionMPI(const std::string &main_dir, const grid_t &grid, const bound_size_t &size) {
  // int np, myid; MPI_GET_INF(np, myid);
  int myid = get_mpi_id();
  int left = grid.mpi_cfg->disp_cells[myid];
  int right = left + grid.mpi_cfg->send_cells[myid];

#ifdef ILLUM
  if (claster_cfg.comm_illum != MPI_COMM_NULL) {
    // todo: это на узле с излучением!
    std::vector<Type> illum;
    GetDirectionDataFromFace(grid.size, 0, grid.Illum, 0.0, illum);
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

// не видет параметры под макросом
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "writer_txt.h"
static inline int WriteFileSolutionOrder(const std::string &main_dir, const grid_t &grid) {

#ifdef ILLUM
  std::vector<Type> illum;
#ifndef SEPARATE_GPU
  GetDirectionDataFromFace(grid.size, 0, grid.Illum, 0.0, illum);
#else
  GetCellDataBySelectedDirection(grid.size, grid.size_dir, 0, grid.Illum, illum);
#endif
  files_sys::bin::WriteSimple(main_dir + F_ILLUM, illum);
  // files_sys::txt::WriteSimple(main_dir + F_ILLUM + ".txt", illum);

#if !defined USE_CUDA
  WRITE_FILE_ELEM((main_dir + F_ENERGY).c_str(), grid.cells, illum_val.energy);

  WRITE_FILE_ELEM((main_dir + F_STREAM).c_str(), grid.cells, illum_val.stream);

  WRITE_FILE_ELEM((main_dir + F_IMPULS).c_str(), grid.cells, illum_val.impuls);

  WRITE_FILE_ELEM((main_dir + F_DIVSTREAM).c_str(), grid.cells, illum_val.div_stream);

  WRITE_FILE_ELEM((main_dir + F_DIVIMPULS).c_str(), grid.cells, illum_val.div_impuls);
#else
  files_sys::bin::WriteSimple((main_dir + F_ENERGY).c_str(), grid.size, grid.energy);

  files_sys::bin::WriteSimple((main_dir + F_STREAM).c_str(), grid.size, grid.stream);

  files_sys::bin::WriteSimple((main_dir + F_IMPULS).c_str(), grid.size, grid.impuls);

  files_sys::bin::WriteSimple((main_dir + F_DIVSTREAM).c_str(), grid.size, grid.divstream);

  files_sys::bin::WriteSimple((main_dir + F_DIVIMPULS).c_str(), grid.size, grid.divimpuls);
#endif
#endif

#if (defined HLLC || defined RHLLC)

  WRITE_FILE_ELEM((main_dir + F_DENSITY).c_str(), grid.cells, phys_val.d);

  WRITE_FILE_ELEM((main_dir + F_PRESSURE).c_str(), grid.cells, phys_val.p);

  WRITE_FILE_ELEM((main_dir + F_VELOCITY).c_str(), grid.cells, phys_val.v);

#endif

  return 0;
}

#pragma GCC diagnostic pop

int files_sys::bin::WriteSolution(const std::string &main_dir, const grid_t &grid) {

  if (get_mpi_id() == 0) {
    return WriteFileSolutionOrder(main_dir, grid);
  }
  return e_completion_success;
}
