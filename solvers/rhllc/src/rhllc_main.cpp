
#if defined SOLVERS && defined RHLLC
#include "rhllc_main.h"

#include "rhllc_calc.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "writer_bin.h"

#include <chrono>
namespace tick = std::chrono;

int rhllc::RunRhllcModule() {
  WRITE_LOG("Start RunRhllcModule()\n");
  grid_t grid;

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_directions_t(0));

  DIE_IF(rhllc::Init(glb_files.hllc_init_value, grid.cells));

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  Timer timer;

  while (t < _hllc_cfg.T) {
    Hllc3dStab(_hllc_cfg.tau, grid);

    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;

    if (cur_timer >= _hllc_cfg.save_timer) {
      WRITE_LOG("t= %lf, step= %d, time_step=%lfs\n", t, res_count, timer.get_delta_time_sec());
      DIE_IF(files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);
      timer.start_timer();
      cur_timer = 0;
    }

    _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);
  }

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid);

  WRITE_LOG("End RunRhllcModule()\n");
  return e_completion_success;
}

#endif //! SOLVERS