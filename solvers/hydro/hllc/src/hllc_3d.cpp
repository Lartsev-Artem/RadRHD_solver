#if defined SOLVERS && defined HLLC
#include "hllc_calc.h"
#include "hllc_flux.h"
#include "hllc_utils.h"

#include <omp.h>

void hllc::Hllc3d(const Type tau, grid_t &grid) {

#pragma omp parallel default(none) firstprivate(tau) shared(grid, glb_files)
  {
    flux_all_t bound_val;
    flux_t phys_val;
    Matrix3 T;
    elem_t *cell;

//потоки
#pragma omp for
    for (int i = 0; i < grid.faces.size(); i++) {
      face_t &f = grid.faces[i];
      BoundConditions(f, grid.cells, bound_val);
      GetFlux(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f);
    }

// ячейки
#pragma omp for
    for (int i = 0; i < grid.cells.size(); i++) {
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

      GetPhysValue(el.conv_val, el.phys_val);
    }
  } // omp

  return;
}

#endif