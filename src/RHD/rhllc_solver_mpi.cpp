#if defined RHLLC && defined USE_MPI
#include "mpi_ext.h"
#include "solvers_struct.h"

#include "hydro_mpi_cfg.h"
#include "rhllc_3d_mpi.h"
#include "rhllc_ini_states.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "reader_json.h"
#include "reader_txt.h"
#include "writer_bin.h"

using namespace rrhd;

int main(int argc, char *argv[])
{
  MPI_START(argc, argv);
  INIT_ENVIRONMENT(argc, argv);

  grid_t grid;
  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err)
  {
    RETURN_ERR("Error geo reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_directions_t(0));
  DIE_IF(rhllc::Init(glb_files.hllc_init_value,
                     rhllc::ini::Soda, grid));

  {
    InitMPiStruct();
    std::vector<int> metis;
    if (files_sys::txt::ReadSimple(glb_files.base_address + F_SEPARATE_METIS, metis))
    {
      RETURN_ERR("Error reading metis \n");
    }

    hydro_mpi::InitMpiConfig(metis, grid);
    metis.clear();
  }

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;
  _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg);

  WRITE_LOG("Init hllc success. tau = %lf, res_count=%d, T=%lf\n",
            _hllc_cfg.tau, res_count, _hllc_cfg.T);

  files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  WRITE_LOG("Start mpi_hllc solver\n");
  MPI_BARRIER(MPI_COMM_WORLD);

  Timer timer;

  while (t < _hllc_cfg.T)
  {
    rhllc_mpi::Hllc3d(_hllc_cfg.tau, grid);

    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;

    if (cur_timer >= _hllc_cfg.save_timer)
    {
      WRITE_LOG("t= %lf, step= %d, time_step=%ld ms\n", t, res_count, timer.get_delta_time_ms());
      DIE_IF(files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);
      timer.start_timer();
      cur_timer = 0;
    }

    _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg);
    MPI_BARRIER(MPI_COMM_WORLD);
  }

  files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid);

  DEINIT_ENVIRONMENT(argc, argv);
  MPI_END;
  return e_completion_success;
}

#endif