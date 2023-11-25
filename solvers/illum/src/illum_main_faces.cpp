#if defined SOLVERS && defined ILLUM
#include "illum_main.h"
#ifdef TRANSFER_CELL_TO_FACE

#include "global_types.h"
#include "illum_add_dir.h"
#include "illum_calc_cpu.h"
#include "illum_init_data.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

int illum::RunIllumFacesModule() {

  WRITE_LOG("Start RunIllumModule()\n");

  if (get_mpi_id() != 0) {
    return e_completion_success;
  }
  grid_t grid;
  grid_directions_t grid_direction;

  std::vector<align_cell_local> vec_x0;

  std::vector<std::vector<graph_pair_t>> sorted_graph;
  std::vector<std::vector<IntId>> sorted_id_bound_face;
  std::vector<std::vector<IntId>> inner_bound_code;

  uint32_t err = 0;
  err |= files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);

  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_direction);

  if (illum::InitRadiationState(glb_files.base_address, grid)) {
    DIE_IF(_solve_mode.class_vtk == e_grid_cfg_radiation); //в иных случаях допускает пропуск инициализации
  }

#ifdef USE_CUDA
  cuda::interface::InitDevice(glb_files.base_address, grid_direction, grid);
#endif

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationFaceTrace(grid_direction.size, glb_files, vec_x0, sorted_graph, sorted_id_bound_face, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

  WRITE_LOG("Start Illum solver()\n");

  cpu::CalculateIllumFace(grid_direction, inner_bound_code,
                          vec_x0, sorted_graph, sorted_id_bound_face, grid);

  cpu::CalculateIllumParam(grid_direction, grid);

#ifdef USE_CUDA
  cuda::interface::CudaWait();
#endif

  files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);

  // if (_solve_mode.max_number_of_iter > 1) //иначе интеграл рассеяния не расчитывался (1-означает считать без рассеяния, но переслать данные на карту)
  // {
  //   D_LD;
  //   additional_direction::SaveInterpolationScattering(glb_files.add_dir_address, grid_direction, grid);
  // }

#ifdef USE_CUDA
  cuda::interface::ClearHost(grid);
  cuda::interface::ClearDevice();
#endif

  WRITE_LOG("End RunIllumModule()\n");
  return e_completion_success;
}

#endif
#endif //! SOLVERS