#if defined SOLVERS && defined ILLUM && defined USE_MPI
#include "illum_main.h"

#include "global_types.h"
#include "illum_add_dir.h"
#include "illum_calc_gpu_async.h"
#include "illum_init_data.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

int illum::RunIllumExtraModule() {

#ifndef USE_CUDA
  D_LD;
#endif

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
  grid.local_Illum.resize(grid_direction.loc_size * grid.size);

  if (illum::InitRadiationState(glb_files.base_address, grid)) {
    DIE_IF(_solve_mode.class_vtk == e_grid_cfg_radiation); //в иных случаях допускает пропуск инициализации
  }

  WRITE_LOG("Init memory\n");

#ifdef USE_CUDA
  cuda::interface::InitDeviceExtraSize(glb_files.base_address, grid_direction, grid);
  cuda::interface::SetStreams();
  WRITE_LOG("Init mpi device\n");
#endif

  gpu_async::NewInitSender(grid_direction, grid); //после инициализации видеокарты, т.к. структура сетки инициализируется и там
  WRITE_LOG("Init mpi sender\n");

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationTrace(grid_direction.size, glb_files, vec_x, face_states, vec_x0, sorted_id_cell, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

  MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

  gpu_async::NewCalculateIllum(grid_direction, face_states, neighbours, inner_bound_code, vec_x0, vec_x, sorted_id_cell, grid);

  WRITE_LOG("end calculate illum\n");

/// \todo ...
#if 0 // def USE_CUDA
  if (get_mpi_id() == 0) {
    cuda::interface::CalculateAllParam(grid_direction, grid);
  }
  cuda::interface::CudaWait();
  cuda::interface::CudaSyncStream(cuda::e_cuda_params);
#endif

  if (get_mpi_id() == 0) {
    files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);
  }

  if (_solve_mode.max_number_of_iter > 1) //иначе интеграл рассеяния не расчитывался
  {
    additional_direction::SaveInterpolationScattering(glb_files.add_dir_address, grid_direction, grid);
  }

/// \todo сделать соответствующие функции
#ifdef USE_CUDA
  // cuda::interface::ClearDevice();
  // cuda::interface::ClearHost(grid);
#endif

  WRITE_LOG("end proc illum\n");
  MPI_BARRIER(MPI_COMM_WORLD);
  return e_completion_success;
}

#endif //! SOLVERS