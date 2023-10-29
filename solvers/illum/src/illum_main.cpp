#if defined SOLVERS && defined ILLUM
#include "illum_main.h"

#include "global_types.h"
#include "illum_add_dir.h"
#include "illum_calc_cpu.h"
#include "illum_init_data.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

int illum::RunIllumModule() {

  WRITE_LOG("Start RunIllumModule()\n");

  if (get_mpi_id() != 0) {
    return e_completion_success;
  }
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
  cuda::interface::InitDevice(glb_files.base_address, grid_direction, grid);
#endif

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationTrace(grid_direction.size, glb_files, vec_x, face_states, vec_x0, sorted_id_cell, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

  cpu::CalculateIllum(grid_direction, face_states, neighbours, inner_bound_code,
                      vec_x0, vec_x, sorted_id_cell, grid);

  cpu::CalculateIllumParam(grid_direction, grid);

#ifdef USE_CUDA
  cuda::interface::CudaWait();
#endif

  files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);

  if (_solve_mode.max_number_of_iter > 1) //иначе интеграл рассеяния не расчитывался (1-означает считать без рассеяния, но переслать данные на карту)
  {
    additional_direction::SaveInterpolationScattering(glb_files.add_dir_address, grid_direction, grid);
  }

#ifdef USE_CUDA
  cuda::interface::ClearHost(grid);
  cuda::interface::ClearDevice();
#endif

  WRITE_LOG("End RunIllumModule()\n");
  return e_completion_success;
}

#endif //! SOLVERS