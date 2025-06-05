#include "global_def.h"
#include "global_value.h"

#include "convert_face_to_cell.h"
#include "mpi_ext.h"
#include "solvers_config.h"
#include "solvers_struct.h"
#include "writer_bin.h"

#ifdef ILLUM
int files_sys::bin::WriteAllIllumDirections(const std::string &main_dir, const grid_t &grid)
{

  if (grid.Illum == nullptr)
    return e_completion_fail;

  std::vector<Type> illum;

  for (size_t i = 0; i < grid.size_dir; i++)
  {
#ifndef SEPARATE_GPU
    GetDirectionDataFromFace(grid.size, i, grid.Illum, 0.0, illum);
#else
    GetCellDataBySelectedDirection(grid.size, grid.size_dir, i, grid.Illum, illum);
#endif
    files_sys::bin::WriteSimple(main_dir + "Illum" + std::to_string(i) + ".bin", illum);
  }

  return e_completion_success;
}
#endif

int files_sys::bin::WriteNormals(const std::string &name_file_normals, std::vector<Normals> &normals)
{
  FILE *f;
  OPEN_FILE(f, name_file_normals.c_str(), "wb");

  int n = (int)normals.size();
  fwrite(&n, sizeof(int), 1, f);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < CELL_SIZE; j++)
    {
      fwrite(&normals[i].n[j], sizeof(Vector3), 1, f);
    }
  }
  fclose(f);

  return e_completion_success;
}

// не видет параметры под макросом
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

int files_sys::bin::WriteSolution(const std::string &main_dir, const grid_t &grid)
{

  if (get_mpi_id() != 0)
  {
    return e_completion_success;
  }

#ifdef ILLUM
  if (grid.Illum != nullptr)
  {
    std::vector<Type> illum;
#ifndef SEPARATE_GPU
    GetDirectionDataFromFace(grid.size, 0, grid.Illum, 0.0, illum);
#else
    GetCellDataBySelectedDirection(grid.size, grid.size_dir, 0, grid.Illum, illum);
#endif
    files_sys::bin::WriteSimple(main_dir + F_ILLUM, illum);
    // files_sys::txt::WriteSimple(main_dir + F_ILLUM + ".txt", illum);
  }

#ifdef ON_FULL_ILLUM_ARRAYS
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

#ifndef RAD_RHD
  files_sys::bin::WriteSimple((main_dir + F_DIVSTREAM).c_str(), grid.size, grid.divstream);

  files_sys::bin::WriteSimple((main_dir + F_DIVIMPULS).c_str(), grid.size, grid.divimpuls);
#endif // RAD_RHD
#endif // USE_CUDA
#endif // ON_FULL_ILLUM_ARRAYS
#endif // ILLUM

#if (defined HLLC || defined RHLLC)

  WRITE_FILE_ELEM((main_dir + F_DENSITY).c_str(), grid.cells, phys_val.d);

  WRITE_FILE_ELEM((main_dir + F_PRESSURE).c_str(), grid.cells, phys_val.p);

  WRITE_FILE_ELEM((main_dir + F_VELOCITY).c_str(), grid.cells, phys_val.v);

#endif
  return e_completion_success;
}

#pragma GCC diagnostic pop

#ifdef USE_MPI

int files_sys::bin::WriteSolutionMPI(const std::string &main_dir, const grid_t &grid)
{

  int myid = get_mpi_id();
  IdType left = grid.loc_shift;
  IdType right = left + grid.loc_size;

#ifdef ILLUM

  if (grid.Illum != nullptr)
  {
    std::vector<Type> illum;
    if (myid == 0)
    {
#ifndef SEPARATE_GPU
      GetDirectionDataFromFace(grid.size, 0, grid.Illum, 0.0, illum);
#else
      GetCellDataBySelectedDirection(grid.size, grid.size_dir, 0, grid.Illum, illum);
#endif
      files_sys::bin::WriteSimple(main_dir + F_ILLUM, illum);
    }
  }

#ifdef ON_FULL_ILLUM_ARRAYS
#if !defined USE_CUDA

  WRITE_FILE_ELEM_MPI(MPI_COMM_WORLD, (main_dir + F_ENERGY).c_str(), grid.cells, illum_val.energy, left, right);
  WRITE_FILE_ELEM_MPI(MPI_COMM_WORLD, (main_dir + F_STREAM).c_str(), grid.cells, illum_val.stream, left, right);
  WRITE_FILE_ELEM_MPI(MPI_COMM_WORLD, (main_dir + F_IMPULS).c_str(), grid.cells, illum_val.impuls, left, right);
  WRITE_FILE_ELEM_MPI(MPI_COMM_WORLD, (main_dir + F_DIVSTREAM).c_str(), grid.cells, illum_val.div_stream, left, right);
  WRITE_FILE_ELEM_MPI(MPI_COMM_WORLD, (main_dir + F_DIVIMPULS).c_str(), grid.cells, illum_val.div_impuls, left, right);
#else
#pragma warning "shift/sizeof(double) or just shift?"
  WriteSimpleMPI(main_dir + F_ENERGY, grid.size, grid.energy, left, right);
  WriteSimpleMPI(main_dir + F_STREAM, grid.size, grid.stream, left, right);
  WriteSimpleMPI(main_dir + F_IMPULS, grid.size, grid.impuls, left, right);
#ifndef RAD_RHD
  WriteSimpleMPI(main_dir + F_DIVSTREAM, grid.size, grid.divstream, left, right);
  WriteSimpleMPI(main_dir + F_DIVIMPULS, grid.size, grid.divimpuls, left, right);
#endif
#endif

#endif
#endif // ILLUM

#if (defined HLLC || defined RHLLC)
  WriteSimpleMPI(main_dir + F_FLUX, grid.size, grid.cells.data() + left, left, right, sizeof(flux_t), MPI_phys_val_t);
#endif

  return 0;
}
#endif
