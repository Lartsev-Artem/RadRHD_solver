#if defined SOLVERS && defined ILLUM && SPECTRUM
#include "illum_main.h"
#include "spectrum_full.h"
#if defined TRANSFER_CELL_TO_FACE && defined SEPARATE_GPU && defined SAVE_FULL_SPECTRUM

#include "spectrum_utils.h"

#include "global_types.h"
#include "global_value.h"

#include "illum_init_data.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"
#include "writer_txt.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

#include "illum_mpi_sender.h"

namespace cuda_sep = cuda::interface::separate_device;

int illum::RunSpectrumModule(int count_states) {

  WRITE_LOG("Start RunSpectrumModule()\n");

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
  err |= files_sys::txt::ReadTableFunc(glb_files.tab_func_address + F_COOLING_FUNC, t_cooling_function);

  if (err) {
    RETURN_ERR("Error reading \n");
  }

  grid.InitFrq();
  grid.InitMemory(grid.cells.size(), grid_direction);
  grid.InitFullPhysData();

  WRITE_LOG("Init memory\n");

  cuda_sep::InitDevice(grid_direction, grid);
  cuda::interface::SetStreams();
  WRITE_LOG("Init mpi device\n");

  spectrum_gpu::InitSender(MPI_COMM_WORLD, grid_direction, grid); //после инициализации видеокарты, т.к. структура сетки инициализируется и там

  //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
  if (files_sys::bin::ReadRadiationFaceTrace(grid_direction.size, glb_files, vec_x0, sorted_graph, sorted_id_bound_face, inner_bound_code))
    RETURN_ERR("Error reading trace part\n");

  WRITE_LOG("Start Illum solver()\n");

  for (int st = 0; st < count_states; st++) {

    if (spectrum::InitPhysState(st, grid) != e_completion_success) {
      break; //не удалось инициализировать газовое состояние
    }
    cuda_sep::SendVelocity(grid);
    memset(grid.scattering, 0, grid.size * grid_direction.loc_size * grid.size_frq);
    memset(grid.Illum, 0, grid.size * grid_direction.size * grid.size_frq);
    for (size_t i = 0; i < grid.inter_coef_all.size(); i++) {
      for (size_t j = 0; j < grid.inter_coef_all[i].size(); j++) {
        for (size_t k = 0; k < grid.inter_coef_all[i][j].size(); k++) {
          grid.inter_coef_all[i][j][k] = 0;
        }
      }
    }

    MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

    full_spectrum::CalculateIllum(grid_direction, inner_bound_code,
                                  vec_x0, sorted_graph, sorted_id_bound_face, grid);

    WRITE_LOG("end calculate illum\n");

    cuda::interface::CudaWait();
    cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);
    cuda::interface::CudaSyncStream(cuda::e_cuda_params);

    MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

    if (get_mpi_id() == 0) {
      std::vector<Type> frq(grid.spectrum.size());
      // IdType cell_shift = grid.size_dir * grid.size_frq * 944;
      for (size_t i = 0; i < frq.size(); i++) {
        frq[i] = (grid.frq_grid[i] + grid.frq_grid[i + 1]) / 2;

        grid.spectrum[i] = 0;
        const int dir = 0;
        for (size_t id = 0; id < grid.size; id++) {
          grid.spectrum[i] += grid.Illum[id * (grid.size_dir * grid.size_frq) + grid.size_frq * dir + i];
        }
        grid.spectrum[i] /= grid.size;
      }

      files_sys::txt::WriteSimple(glb_files.solve_address + std::to_string(st) + F_SPECTRUM,
                                  frq, grid.spectrum);
    }
  }

  cuda_sep::ClearDevice();
  cuda_sep::ClearHost(grid);

  WRITE_LOG("End spectrum()\n");
  return e_completion_success;
}

#endif
#endif //! SOLVERS