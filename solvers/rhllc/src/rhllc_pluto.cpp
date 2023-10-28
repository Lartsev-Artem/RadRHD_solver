#if defined RHLLC && defined SOLVERS
#include "rhllc_calc.h"
#include "rhllc_flux.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "writer_bin.h"
using namespace rhllc;

int hllc_solver(const Type tau, grid_t &grid) {

#pragma omp parallel default(none) firstprivate(tau) shared(grid, glb_files)
  {
    const int size_grid = grid.size;

    flux_all_t bound_val;
    const int size_face = grid.faces.size();

// потоки
#pragma omp for
    for (int i = 0; i < size_face; i++) {
      face_t &f = grid.faces[i];
      BoundConditions(f, grid.cells, bound_val);
      GetFluxStab(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f);
    }

#pragma omp for
    for (int i = 0; i < size_grid; i++) {
      elem_t &el = grid.cells[i];
      flux_t sumF;
      for (int j = 0; j < CELL_SIZE; j++) {
        if (el.geo.sign_n[j]) {
          sumF += grid.faces[el.geo.id_faces[j]].f;
        } else {
          sumF -= grid.faces[el.geo.id_faces[j]].f;
        }
      }
      sumF *= (tau / el.geo.V);
      el.conv_val -= sumF;

      if (GetPhysValueStab(el.conv_val, el.phys_val)) {
        D_LD;
        DIE_IF(PhysPressureFix(el.conv_val, el.phys_val));
      }
    }

  } // omp

  return 0;
}

int RunRhllcPluto() {
  grid_t grid;

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), 0);

  // DIE_IF(rhllc::Init(glb_files.hllc_init_value, grid.cells));

#pragma omp parallel for
  for (int i = 0; i < grid.cells.size(); i++) {
    elem_t &el = grid.cells[i];
    Type x = el.geo.center[0];
    if (x < 0.5) {
      el.phys_val.d = 1;
      el.phys_val.p = 1;
      el.phys_val.v = Vector3(0.9, 0, 0);
    } else {
      el.phys_val.d = 1;
      el.phys_val.p = 10;
      el.phys_val.v = Vector3(0, 0, 0);
    }
  }

#pragma omp parallel for
  for (int i = 0; i < grid.cells.size(); i++) {
    GetConvValueStab(grid.cells[i].phys_val, grid.cells[i].conv_val);
  }

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  while (t < _hllc_cfg.T) {
    hllc_solver(_hllc_cfg.tau, grid);

    if (cur_timer >= _hllc_cfg.save_timer) {
      DIE_IF(files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);

      WRITE_LOG("t= %lf, step= %d\n", t, res_count);
      cur_timer = 0;
    }

    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;
  }

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid);
  return e_completion_success;
}

#endif