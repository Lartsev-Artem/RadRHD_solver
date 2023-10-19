#ifdef ILLUM
#include "illum_add_prebuild.h"
#include "global_value.h"
#include "illum_add_dir.h"

#include "reader_txt.h"

#include "graph_main.h"
#include "trace_main.h"

int illum::additional_direction::PreBuildAddDirections(const int add_directions_cnt) {

  grid_directions_t grid;
  if (files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid))
    return e_completion_fail;

  if (get_mpi_id() == 0) {
    if (MakeDirectionReInterpolation(glb_files.base_address, glb_files.add_dir_address, add_directions_cnt, grid))
      return e_completion_fail;
  }

  MPI_BARRIER(MPI_COMM_WORLD);

  // перенаправляем вывод в новый каталог
  glb_files.name_file_sphere_direction = glb_files.add_dir_address + F_ADDITIONAL_DIRECTION_GRID;
  glb_files.graph_address = glb_files.add_dir_address;
  glb_files.illum_geo_address = glb_files.add_dir_address;
  glb_files.Build();

  if (graph::RunGraphModule())
    return e_completion_fail;

  if (trace::RunTracesModule())
    return e_completion_fail;

  return e_completion_success;
}
#endif //! ILLUM