#if defined SOLVERS && defined HLLC
#include "hllc_main.h"
#include "global_types.h"

#include "hllc_calc.h"
#include "hllc_init.h"
#include "hllc_utils.h"

#include "reader_bin.h"
#include "writer_bin.h"

#include <chrono>
namespace tick = std::chrono;

int hllc::RunHllcModule() {

  grid_t grid;

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_directions_t(0));

  DIE_IF(Init(glb_files.hllc_init_value, grid.cells));

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  auto start_clock = tick::steady_clock::now();

  while (t < _hllc_cfg.T) {
    Hllc3d(_hllc_cfg.tau, grid);

    if (cur_timer >= _hllc_cfg.save_timer) {
      DIE_IF(files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);

      WRITE_LOG("t= %lf, step= %d, time_step=%lf\n", t, res_count, (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
      cur_timer = 0;
    }

    _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);
    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;
  }

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid);
  return e_completion_success;
}

#endif //! &&  HLLC