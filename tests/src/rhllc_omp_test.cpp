/**
 * @file rhllc_test.cpp
 * @brief Тест газодинамического расчёта на релятивистском аналоге задачи сода
 * @version 0.1
 * @date 2023-10-08
 *
 * \details тест проводит газодинамический расчёт. Ожидаемое время выполнения 5минут
 * на выходе бинарные файлы. Требуется дополнительно вызвать rebuild_solve
 *
 * \note требуется поддержка VTK
 *
 */

#include "global_value.h"
#include "reader_json.h"
#include "rhllc_main.h"
#include "solvers_struct.h"

#include "reader_bin.h"

#include "rhllc_flux_stab.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include <omp.h>

static int Rhllc_test()
{
  grid_t grid;

  uint32_t err = 0;
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err)
  {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), grid_directions_t(0));

  DIE_IF(rhllc::Init(glb_files.hllc_init_value, grid.cells));

  Type t = 0.0;
  _hllc_cfg.tau = 0.001;
  _hllc_cfg.T = 1;

  Timer timer_all;
  Timer timer;
  struct
  {
    Type flux;
    Type godunov;
    Type all;
  } times;
  times.all = times.flux = times.godunov = 0;

  while (t < _hllc_cfg.T)
  {
    // rhllc::Hllc3dStab(_hllc_cfg.tau, grid);
    {
#pragma omp parallel default(none) shared(grid, _hllc_cfg, glb_files, timer, times)
      {
        Type tau = _hllc_cfg.tau;
        const int size_grid = grid.size;
        Type max_speed = 0;
        flux_all_t bound_val;
        const int size_face = grid.faces.size();
#pragma omp master
        {
          timer.start_timer();
        }
// потоки
#pragma omp for
        for (int i = 0; i < size_face; i++)
        {
          face_t &f = grid.faces[i];
          rhllc::BoundConditions(f, grid.cells, bound_val);
          max_speed = std::max(max_speed, rhllc::GetFlux(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
        }

#pragma omp master
        {
          times.flux += timer.get_delta_time_ns();
          timer.start_timer();
        }

#pragma omp for
        for (int i = 0; i < size_grid; i++)
        {
          elem_t &el = grid.cells[i];
          flux_t sumF;
          for (int j = 0; j < CELL_SIZE; j++)
          {
            if (el.geo.sign_n[j])
            {
              sumF += grid.faces[el.geo.id_faces[j]].f;
            }
            else
            {
              sumF -= grid.faces[el.geo.id_faces[j]].f;
            }
          }
          sumF *= (tau / el.geo.V);
          el.conv_val -= sumF;

          if (rhllc::GetPhysValue(el.conv_val, el.phys_val))
          {
            DIE_IF(rhllc::PhysPressureFix(el.conv_val, el.phys_val));
          }
        }

#pragma omp master
        {
          times.godunov += timer.get_delta_time_ns();
        }

      } // omp
    } //! Hllc3dStab

    t += _hllc_cfg.tau;
  } // while

  times.all = timer_all.get_delta_time_ns();
  WRITE_LOG("Flux= %lf ms, Godunov= %lf ms, All=%lf ms\n",
            times.flux / 1e6, times.godunov / 1e6, times.all / 1e6);

  return e_completion_success;
} //! RunRhllcModule

int main(int argc, char **argv)
{

  MPI_START(argc, argv);

  std::string file_config = "config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  WRITE_LOG("Start: %s\n", argv[0]);

#pragma omp parallel default(none) shared(glb_files)
  {
#pragma omp master
    {
      WRITE_LOG("Config: threads=%d, maxth=%d\n", omp_get_num_threads(), omp_get_max_threads());
    }
  }
  _solve_mode.start_point = 0;
  _hllc_cfg.CFL = 0.7;
  _hllc_cfg.h_min = 1e-5;
  _hllc_cfg.save_timer = 0.1;

  Rhllc_test();

  MPI_END;
  return e_completion_success;
}
