/**
 * @file rhllc_mpi_test.cpp
 * @author your name (you@domain.com)
 * @brief Тест газодинамического расчёта с технологией MPI + OpenMP
 * @note Разбиение по группам должно быть организовано так, чтобы каждая подобласть имела
 * не более двух соседей (расширение возможно, но не реализовано)
 * @version 0.1
 * @date 2025-06-05
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "global_value.h"

#include "reader_json.h"
#include "reader_txt.h"
#include "reader_bin.h"

#include "rhllc_flux.h"
#include "rhllc_ini_states.h"
#include "rhllc_bound_cond.h"
#include "rhllc_utils.h"
#include "rhllc_3d_mpi.h"
#include "rhllc_3d.h"

#include "writer_bin.h"

#include <omp.h>
using namespace rrhd;

static int init_state(grid_t &grid)
{
  {
    std::vector<Vector3> centers;
    if (files_sys::bin::ReadSimple(glb_files.base_address + F_CENTERS, centers))
      RETURN_ERR("Default rhllc value not set\n");

#pragma omp parallel for
    for (int i = 0; i < grid.size; i++)
    {
      elem_t &el = grid.cells[i];

      Type x = centers[i][0];
      if (x < 0.5)
      {
        el.phys_val.d = 1;
        el.phys_val.p = 1;
        el.phys_val.v = Vector3(0.9, 0, 0);
      }
      else
      {
        el.phys_val.d = 1;
        el.phys_val.p = 10;
        el.phys_val.v = Vector3(0, 0, 0);
      }
    }

    rhllc::HllcPhysToConv(grid);
  }

  return e_completion_success;
}

static int SingleRhllc(int argc, char **argv, grid_t &grid)
{
  glb_files.solve_address += "_def";
  WRITE_LOG("Start SingleRhllc()\n");

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err)
  {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_directions_t(0));

  DIE_IF(init_state(grid));

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  Timer timer;
  Timer timer_all;

  while (t < _hllc_cfg.T)
  {
    rhllc::Hllc3d(_hllc_cfg.tau, grid);

    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;

    if (cur_timer >= _hllc_cfg.save_timer)
    {
      WRITE_LOG("t= %lf, step= %d, time_step=%lld ms\n", t, res_count, timer.get_delta_time_ms());
      DIE_IF(files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);
      timer.start_timer();
      cur_timer = 0;
    }

    _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg);
  }

  WRITE_LOG("Full_time= %lld ms\n", timer_all.get_delta_time_ms());

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid);

  WRITE_LOG("End SingleRhllc()\n");
  return e_completion_success;
}

int main(int argc, char **argv)
{
  MPI_START(argc, argv);
  INIT_ENVIRONMENT(argc, argv);

  grid_t default_grid;
  if (get_mpi_id() == 0)
  {
    SingleRhllc(argc, argv, default_grid);
  }
  MPI_BARRIER(MPI_COMM_WORLD);

  INIT_ENVIRONMENT(argc, argv);

  grid_t grid;

  uint32_t err = 0;
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err)
  {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_directions_t(0));

  DIE_IF(init_state(grid));

  // mpi_init:
  {
    InitMPiStruct();

    std::vector<int> metis;
    if (files_sys::txt::ReadSimple(glb_files.base_address + F_SEPARATE_METIS, metis))
      RETURN_ERR("Error reading metis \n");

    rhllc_mpi::InitMpiConfig(metis, grid);
    metis.clear();
  }

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  rhllc::max_signal_speed = 1;
  _hllc_cfg.tau = rhllc::GetTimeStep(_hllc_cfg);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  MPI_BARRIER(MPI_COMM_WORLD);

  Timer timer;
  Timer timer_all;

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

  WRITE_LOG("Full_time= %lld ms\n", timer_all.get_delta_time_ms());
  files_sys::bin::WriteSolutionMPI(glb_files.solve_address + std::to_string(res_count++), grid);

  rhllc_mpi::StartPhysCast(grid.mpi_cfg, grid);
  rhllc_mpi::SyncPhysCast(grid.mpi_cfg);
  if (get_mpi_id() == 0)
  {
    Vector3 sum = Vector3::Zero();
    for (size_t i = 0; i < grid.size; i++)
    {
      sum[0] += fabs(grid.cells[i].phys_val.d - default_grid.cells[i].phys_val.d);
      sum[1] += (grid.cells[i].phys_val.v - default_grid.cells[i].phys_val.v).norm();
      sum[2] += fabs(grid.cells[i].phys_val.p - default_grid.cells[i].phys_val.p);
    }

    WRITE_LOG("Result Sum_error (mpi vs single): %lf\n", sum.norm());
  }

  DEINIT_ENVIRONMENT(argc, argv);
  MPI_END;
  return e_completion_success;
}
