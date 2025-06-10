#if defined RHLLC && defined USE_MPI
#include "rhllc_3d_mpi.h"

#include "global_value.h"
#include "reader_txt.h"

#include "rhllc_bound_cond.h"
#include "rhllc_flux.h"
#include "rhllc_utils.h"

using namespace rrhd;

void rhllc_mpi::Hllc3d(const Type tau, grid_t &grid) {
  rhllc::max_signal_speed = 0;
  MPI_Request rq_max_speed = MPI_REQUEST_NULL;
  const int myid = get_mpi_id();

  hydro_mpi::StartExchangeBoundaryCells(grid.mpi_cfg);

#pragma omp parallel default(none) firstprivate(myid) shared(grid, glb_files, rhllc::max_signal_speed)
  {
    const int size_grid = grid.size;
    Type max_speed = 0;
    flux_all_t bound_val;
    const int size_face = grid.faces.size();

    int reg_f_begin = grid.mpi_cfg->maps[myid].face.reg_l;
    int reg_f_end = grid.mpi_cfg->maps[myid].face.reg_r;

// потоки
#pragma omp for // nowait
    for (int i = reg_f_begin; i < reg_f_end; i++) {
      face_t &f = grid.faces[i];
      RRHD_ASSERT(f.geo.is_regular);
      rhllc::BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, rhllc::GetFlux(
                                          grid.cells[f.geo.id_l].conv_val, bound_val.conv_val,
                                          grid.cells[f.geo.id_l].phys_val, bound_val.phys_val,
                                          f));
    }
// @note если здесь будем зависать, можно вычислять регулярные ячейки до синхронизации
#pragma omp master
    {
      hydro_mpi::SyncExchangeBoundaryCells(grid.mpi_cfg);
    }

#pragma omp barrier

    int f_begin = grid.mpi_cfg->maps[myid].face.left_id;
    int f_end = grid.mpi_cfg->maps[myid].face.right_id;
// потоки на граничных ячейках узла
#pragma omp for
    for (int i = f_begin; i < reg_f_begin; i++) {
      face_t &f = grid.faces[i];
      // DIE_IF(f.geo.is_regular);

      rhllc::BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, rhllc::GetFlux(
                                          grid.cells[f.geo.id_l].conv_val, bound_val.conv_val,
                                          grid.cells[f.geo.id_l].phys_val, bound_val.phys_val,
                                          f));
    }
#pragma omp for
    for (int i = reg_f_end; i < f_end; i++) {
      face_t &f = grid.faces[i];
      // DIE_IF(f.geo.is_regular);

      rhllc::BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, rhllc::GetFlux(
                                          grid.cells[f.geo.id_l].conv_val, bound_val.conv_val,
                                          grid.cells[f.geo.id_l].phys_val, bound_val.phys_val,
                                          f));
    }
    // ищем максимум по потокам
    if (max_speed > rhllc::max_signal_speed) {
#pragma omp critical
      {
        rhllc::max_signal_speed = std::max(rhllc::max_signal_speed, max_speed);
      }
    }
  }

  // ищем максимум по всем узлам
  Type max_speed = rhllc::max_signal_speed;
  MPI_Iallreduce(&max_speed, &rhllc::max_signal_speed, 1, MPI_DOUBLE, MPI_MAX, grid.mpi_cfg->comm, &rq_max_speed);

#pragma omp parallel default(none) firstprivate(tau, myid) shared(grid, glb_files)
  {
    int start = grid.mpi_cfg->maps[myid].cell.left_id;
    int end = grid.mpi_cfg->maps[myid].cell.right_id;

#pragma omp for
    for (int i = start; i < end; i++) {
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

      if (rhllc::GetPhysValue(el.conv_val, el.phys_val)) {
        DIE_IF(rhllc::PhysPressureFix(el.conv_val, el.phys_val));
      }
    }

  } // omp

  MPI_Wait(&rq_max_speed, MPI_STATUS_IGNORE);
  rhllc_mpi_log("rhllc_speed=%lf\n", rhllc::max_signal_speed)
}

#endif //!  RHLLC &&  SOLVERS && USE_MPI