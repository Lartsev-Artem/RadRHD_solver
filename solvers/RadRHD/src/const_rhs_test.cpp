#if defined SOLVERS && defined RAD_RHD && defined USE_CUDA
#include "global_value.h"
#include "radRHD_main.h"
#include "radRHD_utils.h"

#include "rhllc_calc.h"
#include "rhllc_flux_stab.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

#include <chrono>
namespace tick = std::chrono;

int rad_rhd::RadRHD_ConstRadStateTest() {
  WRITE_LOG("Start RadRHDTest()\n");
  grid_t grid;

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), 0);

  DIE_IF(rhllc::Init(glb_files.hllc_init_value, grid.cells));

  grid_directions_t grid_direction;
  files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);
  cuda::interface::InitDevice(glb_files.base_address, grid_direction, grid);

#pragma omp parallel for
  for (int i = 0; i < grid.size; i++) {
    if (grid.cells[i].geo.center[0] < 0.5) {
      grid.energy[i] = 10;
      grid.stream[i] = Vector3(5, 0, 0);
      grid.impuls[i] = Matrix3::Zero();
      grid.impuls[i](0, 0) = grid.impuls[i](1, 1) = grid.impuls[i](2, 2) = grid.energy[i];
    } else {
      grid.energy[i] = 0;
      grid.stream[i] = Vector3(0, 0, 0);
      grid.impuls[i] = Matrix3::Zero();
      grid.impuls[i](0, 0) = grid.impuls[i](1, 1) = grid.impuls[i](2, 2) = grid.energy[i] / 3;
    }
    grid.cells[i].phys_val.d = 1;
    grid.cells[i].phys_val.p = 10;
    grid.cells[i].phys_val.v = Vector3::Zero();
    grid.cells[i].illum_val.absorp_coef = 0.5;
    grid.cells[i].illum_val.scat_coef = 0.5;
  }

  rhllc::HllcPhysToConv(grid.cells);

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg, grid.cells);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  auto start_clock = tick::steady_clock::now();

  while (t < _hllc_cfg.T) {
    rhllc::Hllc3dStab(_hllc_cfg.tau, grid);

    if (cur_timer >= _hllc_cfg.save_timer) {
      DIE_IF(files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);

      WRITE_LOG("t= %lf, step= %d, time_step=%lf\n", t, res_count, (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
      cur_timer = 0;
    }

#pragma omp parallel for
    for (int cell = 0; cell < grid.size; cell++) {

      Vector4 G;
      GetRadSource(cell, grid, G);

      constexpr Type ds = 1;
      grid.cells[cell].conv_val.p += ds * _hllc_cfg.tau * G[0];
      grid.cells[cell].conv_val.v[0] += ds * _hllc_cfg.tau * G[1];
      grid.cells[cell].conv_val.v[1] += ds * _hllc_cfg.tau * G[2];
      grid.cells[cell].conv_val.v[2] += ds * _hllc_cfg.tau * G[3];

      rhllc::GetPhysValueStab(grid.cells[cell].conv_val, grid.cells[cell].phys_val);
    }

    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;
    _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg, grid.cells);
  }

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid);

  WRITE_LOG("End RadRHDTest()\n");
  return e_completion_success;
}

#endif //! SOLVERS
