#if defined RHLLC && defined SOLVERS
//&& !defined RHLLC_MPI && NUMBER_OF_MEASUREMENTS == 3
#include "rhllc_calc.h"
#include "rhllc_flux.h"
#include "rhllc_utils.h"

#include <omp.h>

void rhllc::Hllc3d(const Type tau, grid_t &grid) {

#pragma omp parallel default(none) firstprivate(tau) shared(grid, glb_files)
  {
    const int size_grid = grid.size;

#ifdef ILLUM
    // востановление физических переменных
#pragma omp for
    for (int i = 0; i < size_grid; i++) {
      int back = GetPhysValueSave(grid.cells[i].conv_val, grid.cells[i].phys_val);
      if (back) {
        D_LD;
        if (back == 2) {

        } else {
          // printf("try id= %d\n", i);
        }
        // если был персчёт

        GetConvValue(grid.cells[i].phys_val, grid.cells[i].conv_val);
      }
    }
#endif

    flux_all_t bound_val;
    const int size_face = grid.faces.size();

    // потоки
#pragma omp for
    for (int i = 0; i < size_face; i++) {
      face_t &f = grid.faces[i];
      BoundConditions(f, grid.cells, bound_val);
      GetFlux(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f);
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
      el.conv_val -= sumF * (tau / el.geo.V);

      GetPhysValue(el.conv_val, el.phys_val); // востановление физических переменных
    }

  } // omp
}

#if 0
int MPI_RHLLC_3d(const int myid, const Type tau, grid_t &grid) {
#ifdef ILLUM

  if (myid == 0) {
#pragma omp parallel default(none) shared(tau, grid, glb_files)
    {
      const int size_grid = grid.size;
#ifdef ILLUM
      // востановление физических переменных
#pragma omp for
      for (int i = 0; i < size_grid; i++) {
        int back = rhllc_get_phys_value_ost1098(grid.cells[i].conv_val, grid.cells[i].phys_val);
        if (back) {
          if (back == 2) {
            D_LD;
          } else {
            // printf("try id= %d\n", i);
          }
          // если был пернсчёт
          rhllc_get_conv_value_ost1098(grid.cells[i].phys_val, grid.cells[i].conv_val);
        }
      }
#endif

      flux_all_t bound_val;
      const int size_face = grid.faces.size();

      // потоки
#pragma omp for
      for (int i = 0; i < size_face; i++) {
        face_t &f = grid.faces[i];
        rhllc_bound_conditions(f, grid.cells, bound_val);
        flux_t_calc(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f);
      }
    }

  } // myid ==0

  auto calc{[&grid, myid, tau](const int left, const int right) {
    if (myid == 0) {

#pragma omp for
      for (int i = left; i < right; i++) {
        elem_t &el = grid.cells[i];
        flux_t sumF;
        for (int j = 0; j < CELL_SIZE; j++) {
          if (el.geo.sign_n[j]) {
            sumF += grid.faces[el.geo.id_faces[j]].f;
          } else {
            sumF -= grid.faces[el.geo.id_faces[j]].f;
          }
        }
        el.conv_val -= sumF * (tau / el.geo.V);

        rhllc_get_phys_value_ost1098(el.conv_val, el.phys_val); // востановление физических переменных
        phys_local[i] = el.phys_val;
      }
    }
  }};

  const int size_grid = grid.size;

  for (int i = 0; i < disp_hllc[0].size(); i++) {
#pragma omp parallel default(none) firstprivate(i, tau) shared(grid, send_hllc, disp_hllc, phys_local, calc)
    {
      calc(disp_hllc[0][i], disp_hllc[0][i] + send_hllc[0][i]);
    }

    SendPhysValue(phys_local.data() + disp_hllc[0][i], send_hllc[0][i], i);
  }

#else
#pragma error "bad prj config"
#endif // USE_MPI

  return 0;
}

#endif
#endif // RHLLC_3d
