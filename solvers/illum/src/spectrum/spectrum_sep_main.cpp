#if defined SOLVERS && defined ILLUM
#include "illum_main.h"
#if defined TRANSFER_CELL_TO_FACE && defined SEPARATE_GPU && defined SPECTRUM

#include "global_types.h"
#include "illum_init_data.h"
#include "spec_all.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"
#include "illum_mpi_sender.h"

int illum::spec::RunIllumSpectrumModule() {

  WRITE_LOG("Start RunIllumSpectrumModuleMpiSep()\n");

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

  WRITE_LOG("Init memory\n");

  // cuda_sep::InitDevice(grid_direction, grid);
  cuda::interface::SetStreams();
  WRITE_LOG("Init mpi device\n");

  spectrum_gpu::InitSender(MPI_COMM_WORLD, grid_direction, grid); //после инициализации видеокарты, т.к. структура сетки инициализируется и там

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationFaceTrace(grid_direction.size, glb_files, vec_x0, sorted_graph, sorted_id_bound_face, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

  WRITE_LOG("Start Illum solver()\n");

  spec::CalculateIllumFaceMpi(grid_direction, inner_bound_code,
                              vec_x0, sorted_graph, sorted_id_bound_face, grid);

  WRITE_LOG("end calculate illum\n");

  cuda::interface::CudaWait();
  cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);
  cuda::interface::CudaSyncStream(cuda::e_cuda_params);

  MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

  if (get_mpi_id() == 0) {
    files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);
  }

  //  cuda_sep::ClearDevice();
  //  cuda_sep::ClearHost(grid);

  WRITE_LOG("End RunIllumSpectrumModule()\n");
  return e_completion_success;
}

#endif
#endif //! SOLVERS