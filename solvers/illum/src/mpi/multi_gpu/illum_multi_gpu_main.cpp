#if defined SOLVERS && defined ILLUM && defined USE_MPI && defined USE_CUDA
#include "illum_main.h"

#if defined SEPARATE_GPU && !defined SPECTRUM

#include "global_types.h"
#include "illum_calc_gpu_async.h"
#include "illum_init_data.h"
#include "illum_mpi_sender.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

namespace cuda_sep = cuda::interface::separate_device;

int illum::RunIllumMultiGpuModule() {

  WRITE_LOG("Start RunIllumMultiGpuModule()\n");

  grid_t grid;
  grid_directions_t grid_direction;

  std::vector<align_cell_local> vec_x0;

  std::vector<std::vector<graph_pair_t>> sorted_graph;
  boundary_faces_by_directions_t boundary_faces;
  std::vector<std::vector<IntId>> inner_bound_code;

  uint32_t err = 0;
  err |= files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);

  if (err) {
    RETURN_ERR("Error reading \n");
  }

  grid.InitMemory(grid.cells.size(), grid_direction);
  grid.InitFullPhysData();

  if (illum::InitRadiationState(glb_files.base_address, grid)) {
    DIE_IF(_solve_mode.class_vtk == e_grid_cfg_radiation); //в иных случаях допускает пропуск инициализации
  }

  WRITE_LOG("Init memory\n");

  cuda_sep::InitDevice(grid_direction, grid);
  cuda::interface::SetStreams();
  WRITE_LOG("Init mpi device\n");

  separate_gpu::InitSender(MPI_COMM_WORLD, grid_direction, grid); //после инициализации видеокарты, т.к. структура сетки инициализируется и там

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationFaceTrace(grid_direction.size, glb_files, vec_x0, sorted_graph, boundary_faces, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

#pragma omp parallel for
  for (int i = 0; i < grid.size; i++) {
    grid.cells[i].cell_data->alpha = 1;
    grid.cells[i].cell_data->betta = 0.5;
    grid.cells[i].cell_data->T = 1e7;
  }

  MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

  separate_gpu::CalculateIllum(grid_direction, inner_bound_code,
                               vec_x0, sorted_graph, boundary_faces, grid);

  WRITE_LOG("end calculate illum\n");

  cuda::interface::CudaWait();
  cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);

  // if (get_mpi_id() == 0) {
  //   cuda_sep::CalculateAllParamAsync(grid_direction, grid, cuda::e_cuda_params);
  // }

  cuda::interface::CudaSyncStream(cuda::e_cuda_params);

  MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

  if (get_mpi_id() == 0) {
    // files_sys::bin::WriteSimple(glb_files.solve_address + "scat.bin", grid.size, grid.scattering);
    files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);
  }
  // files_sys::bin::WriteSimple(glb_files.solve_address + "FullIllum" + std::to_string(get_mpi_id()) + "bin", grid.size * grid_direction.size, grid.Illum);

  cuda_sep::ClearDevice();
  cuda_sep::ClearHost(grid);

  WRITE_LOG("end proc illum\n");
  MPI_BARRIER(MPI_COMM_WORLD);
  return e_completion_success;
}

#endif //! SEPARATE_GPU
#endif //! SOLVERS