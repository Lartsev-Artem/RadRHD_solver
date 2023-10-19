#if defined SOLVERS && defined ILLUM
#include "illum_add_main.h"
#include "global_value.h"

#include "global_types.h"
#include "illum_add_dir.h"
#include "illum_calc_cpu.h"
#include "illum_init_data.h"

#include "convert_face_to_cell.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

int illum::additional_direction::RunModule() {

  //переключаемся на новую сферу направлений
  glb_files.name_file_sphere_direction = glb_files.add_dir_address + F_ADDITIONAL_DIRECTION_GRID;
  glb_files.illum_geo_address = glb_files.add_dir_address;
  glb_files.graph_address = glb_files.add_dir_address;
  glb_files.Build();

  grid_t grid;
  grid_directions_t grid_direction;

  std::vector<IntId> neighbours;

  std::vector<BasePointTetra> vec_x;
  std::vector<std::vector<State>> face_states;
  std::vector<std::vector<cell_local>> vec_x0;
  std::vector<std::vector<IntId>> sorted_id_cell;
  std::vector<std::vector<IntId>> inner_bound_code;

  uint32_t err = 0;
  err |= files_sys::bin::ReadSimple(glb_files.name_file_neigh, neighbours);
  err |= files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);
  err |= files_sys::bin::ReadRadiationTrace(grid_direction.size, glb_files, vec_x, face_states, vec_x0, sorted_id_cell, inner_bound_code);

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);

  if (err) {
    RETURN_ERR("Error reading \n");
  }

  grid.InitMemory(grid.cells.size(), grid_direction.size);

  if (illum::InitRadiationState(glb_files.base_address, grid)) {
    DIE_IF(_solve_mode.class_vtk == e_grid_cfg_radiation); //в иных случаях допускает пропуск инициализации
  }

#ifdef USE_CUDA
  grid.Illum = new Type[grid_direction.loc_size * grid.size * CELL_SIZE];
  grid.scattering = new Type[grid_direction.loc_size * grid.size];
#endif

  int loc_dir = 0;
  for (size_t i = grid_direction.loc_shift; i < grid_direction.loc_shift + grid_direction.loc_size; i++, loc_dir++) {
    files_sys::bin::ReadSimple(glb_files.add_dir_address + F_SCATTERING + std::to_string(i) + ".bin", &grid.scattering[loc_dir * grid.size]);
  }

  cpu::CalculateAdditionalIllum(grid_direction, face_states, neighbours, inner_bound_code,
                                vec_x0, vec_x, sorted_id_cell, grid);

  std::vector<Type> illum(grid.size);
  for (int i = 0; i < grid_direction.loc_size; i++) {
    GetDirectionDataFromFace(grid.size, i, grid.Illum, 0.0, illum);
    files_sys::bin::WriteSimple(glb_files.add_dir_address + F_ILLUM + std::to_string(grid_direction.loc_shift + i) + ".bin", illum);
  }

#ifdef USE_CUDA
  delete[] grid.Illum;
  delete[] grid.scattering;
#endif

  return e_completion_success;
}

#endif //! SOLVERS