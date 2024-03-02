#if defined RHLLC && defined SOLVERS
#include "rhllc_calc.h"
#include "rhllc_flux_stab.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "writer_bin.h"

void rhllc::Hllc3dStab(const Type tau, grid_t &grid) {

  rhllc::max_signal_speed = 0;

#pragma omp parallel default(none) firstprivate(tau) shared(grid, glb_files, rhllc::max_signal_speed)
  {
    const int size_grid = grid.size;
    Type max_speed = 0;
    flux_all_t bound_val;
    const int size_face = grid.faces.size();
// потоки
#pragma omp for
    for (int i = 0; i < size_face; i++) {
      face_t &f = grid.faces[i];
      BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, GetFluxStab(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
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
        DIE_IF(PhysPressureFix(el.conv_val, el.phys_val));
      }
    }

    if (max_speed > rhllc::max_signal_speed) {
#pragma omp critical
      {
        rhllc::max_signal_speed = std::max(max_signal_speed, max_speed);
      }
    }

  } // omp
  WRITE_LOG("speed=%lf\n", rhllc::max_signal_speed)
}
#endif