#if defined SOLVERS && defined ILLUM && defined USE_MPI && defined USE_CUDA
#include "illum_main.h"

#if defined SEPARATE_GPU && defined SPECTRUM && !defined SAVE_FULL_SPECTRUM
#include "spectrum_calc.h"
#include "spectrum_utils.h"

#include "global_value.h"

#include "illum_init_data.h"
#include "illum_mpi_sender.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_txt.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

namespace cuda_sep = cuda::interface::separate_device;

int illum::RunSpectrumModule(int count_states) {

  WRITE_LOG("Start RunSpectrumModule()\n");

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
  err |= files_sys::txt::ReadTableFunc(glb_files.tab_func_address + F_COOLING_FUNC, t_cooling_function);

  if (err) {
    RETURN_ERR("Error reading \n");
  }

  grid.InitMemory(grid.cells.size(), grid_direction);
  grid.InitFullPhysData();
  grid.InitFrq();

  WRITE_LOG("Init memory\n");

  cuda_sep::InitDevice(grid_direction, grid);
  cuda::interface::SetStreams();
  WRITE_LOG("Init mpi device\n");

  separate_gpu::InitSender(MPI_COMM_WORLD, grid_direction, grid); //после инициализации видеокарты, т.к. структура сетки инициализируется и там

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationFaceTrace(grid_direction.size, glb_files, vec_x0, sorted_graph, boundary_faces, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

  for (int st = 0; st < count_states; st++) {

    if (spectrum::InitPhysState(st, grid) != e_completion_success) {
      break; //не удалось инициализировать газовое состояние
    }
    cuda_sep::SendVelocity(grid);

    MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

    spectrum::CalculateSpectrum(grid_direction, inner_bound_code,
                                vec_x0, sorted_graph, boundary_faces, grid);

    if (get_mpi_id() == 0) {
      std::vector<Type> frq(grid.spectrum.size());
      for (size_t i = 0; i < frq.size(); i++) {
        frq[i] = (grid.frq_grid[i] + grid.frq_grid[i + 1]) / 2;
      }

      files_sys::txt::WriteSimple(glb_files.solve_address + std::to_string(st) + F_SPECTRUM,
                                  frq, grid.spectrum);
    }
  }
  WRITE_LOG("end calculate spectrum\n");
  MPI_BARRIER(MPI_COMM_WORLD);

  cuda_sep::ClearDevice();
  cuda_sep::ClearHost(grid);

  WRITE_LOG("end proc spectrum\n");
  return e_completion_success;
}

#endif //! SEPARATE_GPU && SPECTRUM
#endif //! SOLVERS

//#include "plunk.h"
// if (get_mpi_id() == 0) {
//   std::vector<Type> frq(grid.spectrum.size());
//   for (size_t i = 0; i < frq.size(); i++) {
//     frq[i] = (grid.frq_grid[i] + grid.frq_grid[i + 1]) / 2;
//     grid.spectrum[i] = exp(B_Plank_log(40000000, grid.frq_grid[i + 1], grid.frq_grid[i]));
//   }
//   files_sys::txt::WriteSimple(glb_files.solve_address + "init_spectrum.txt",
//                               frq, grid.spectrum);
// }
