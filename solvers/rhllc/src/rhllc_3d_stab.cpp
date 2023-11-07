#if defined RHLLC && defined SOLVERS
#include "rhllc_calc.h"
#include "rhllc_flux_stab.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "writer_bin.h"

#pragma GCC target("avx2")
#pragma GCC optimize("O3")

#include <bits/stdc++.h>
#include <x86intrin.h>

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
      max_speed = std::max(max_speed, GetFluxStabVec(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
    }

#if 0
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

      // WRITE_LOG("%d %lf %lf %lf %lf %lf\n", i, el.conv_val.d,
      // el.conv_val.v[0], el.conv_val.v[1], el.conv_val.v[2], el.conv_val.p);
      if (GetPhysValueStab(el.conv_val, el.phys_val)) {
        DIE_IF(PhysPressureFix(el.conv_val, el.phys_val));
      }
    }
#else
#pragma omp for
    for (int i = 0; i < size_grid; i++) {
      elem_t &el = grid.cells[i];
      Type sumF = 0;
      __m256d F = _mm256_set1_pd(0);
      for (int j = 0; j < CELL_SIZE; j++) {
        if (el.geo.sign_n[j]) {
          F = _mm256_add_pd(F, _mm256_loadu_pd(&grid.faces[el.geo.id_faces[j]].f.d));
          sumF += grid.faces[el.geo.id_faces[j]].f.p;
        } else {
          F = _mm256_sub_pd(F, _mm256_loadu_pd(&grid.faces[el.geo.id_faces[j]].f.d));
          sumF -= grid.faces[el.geo.id_faces[j]].f.p;
        }
      }
      Type tv = tau / el.geo.V;
      F = _mm256_mul_pd(F, _mm256_set1_pd(tv));
      sumF *= (tv);

      _mm256_storeu_pd(&el.conv_val.d, _mm256_sub_pd(_mm256_loadu_pd(&el.conv_val.d), F));
      el.conv_val.p -= sumF;

      // WRITE_LOG("%d %lf %lf %lf %lf %lf\n", i, el.conv_val.d,
      // el.conv_val.v[0], el.conv_val.v[1], el.conv_val.v[2], el.conv_val.p);

      if (GetPhysValueStab(el.conv_val, el.phys_val)) {
        DIE_IF(PhysPressureFix(el.conv_val, el.phys_val));
      }
    }
#endif

    if (max_speed > rhllc::max_signal_speed) {
#pragma omp critical
      {
        rhllc::max_signal_speed = std::max(max_signal_speed, max_speed);
      }
    }

  } // omp
}
#endif