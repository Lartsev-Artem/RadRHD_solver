#if defined SOLVERS && defined ILLUM && RAD_RHD && defined USE_MPI
#include "solvers_config.h"
#include "timer.h"

#if defined TRANSFER_CELL_TO_FACE && defined SEPARATE_GPU && !defined SPECTRUM
#include "radRHD_main.h"
#include "radRHD_utils.h"

#include "global_types.h"
#include "global_value.h"

#include "illum_calc_gpu_async.h"
#include "illum_init_data.h"
#include "illum_mpi_sender.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

#include "rhllc_calc.h"
#include "rhllc_flux_stab.h"
#include "rhllc_init.h"
#include "rhllc_mpi.h"
#include "rhllc_utils.h"

#ifndef USE_CUDA
#error "Need gpu support"
#endif

namespace cuda_sep = cuda::interface::separate_device;

int rad_rhd::RunRadRHDMpiModule() {

  WRITE_LOG("Start RunRadRHDMpiModule()\n");

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

  grid.InitMemory(grid.cells.size(), grid_direction);
  grid.InitFullPhysData();
  WRITE_LOG("Init memory\n");

  InitMPiStruct();

  // hllc init
  {
    int np = get_mpi_np();
    DIE_IF(rhllc::Init(glb_files.hllc_init_value, grid.cells));

    std::vector<int> metis;
    if (files_sys::txt::ReadSimple(glb_files.base_address + F_SEPARATE_METIS, metis)) {
      RETURN_ERR("Error reading metis \n");
    }

    grid.mpi_cfg = new mpi_hllc_t;
    rhllc_mpi::InitMpiConfig(metis, grid, grid.mpi_cfg);
    metis.clear();
  }

  // illum init
  {
    cuda_sep::InitDevice(grid_direction, grid);
    cuda::interface::SetStreams();
    WRITE_LOG("Init mpi device\n");

    illum::separate_gpu::InitSender(MPI_COMM_WORLD, grid_direction, grid); //после инициализации видеокарты, т.к. структура сетки инициализируется и там
    WRITE_LOG("Init mpi sender\n");

    //перенесено ниже,т.к. читается долго, а потенциальных ошибок быть не должно
    if (files_sys::bin::ReadRadiationFaceTrace(grid_direction.size, glb_files, vec_x0, sorted_graph, sorted_id_bound_face, inner_bound_code))
      RETURN_ERR("Error reading trace part\n");
  }

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;
  files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  Timer timer;

  MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

  const int myid = get_mpi_id();

  while (t < _hllc_cfg.T) {

    timer.start_timer();

    rhllc_mpi::Hllc3dStab(_hllc_cfg.tau, grid);
    rhllc_mpi::StartPhysCast(*grid.mpi_cfg, grid);
    rhllc_mpi::SyncAndCalcPhysCast(grid);

    illum::separate_gpu::CalculateIllum(grid_direction, inner_bound_code, vec_x0, sorted_graph, sorted_id_bound_face, grid);

    cuda::interface::CudaSyncStream(cuda::e_cuda_params);
    cuda::interface::CudaWait();
    D_L;
    rhllc_mpi::AddRadFlux(grid);
    D_L;
    rhllc_mpi::HllcConvToPhys(grid);
    D_L;
    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;

    if (cur_timer >= _hllc_cfg.save_timer) {
      DIE_IF(files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);

      WRITE_LOG("t= %lf, time_step= %d\n", t, res_count);
      cur_timer = 0;
    }

    _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg, grid.cells);

    WRITE_LOG("it= %lf, time_step=%lf\n", t, timer.get_delta_time_sec());

    MPI_BARRIER(MPI_COMM_WORLD);
  }

  WRITE_LOG("end calculate illum\n");

  cuda::interface::CudaWait();
  cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);

  // if (get_mpi_id() == 0) {
  //   cuda_sep::CalculateAllParamAsync(grid_direction, grid, cuda::e_cuda_params);
  // }

  cuda::interface::CudaSyncStream(cuda::e_cuda_params);

  MPI_BARRIER(MPI_COMM_WORLD); //ждём пока все процессы проинициализируют память

  DIE_IF(files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);

  cuda_sep::ClearDevice();
  cuda_sep::ClearHost(grid);

  WRITE_LOG("end proc illum\n");
  MPI_BARRIER(MPI_COMM_WORLD);
  return e_completion_success;
}

#endif //! TRANSFER_CELL_TO_FACE
#endif //! SOLVERS
