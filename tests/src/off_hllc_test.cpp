/**
 * @file hllc_test.cpp
 * @brief Тест газодинамического расчёта на задаче сода
 * @version 0.1
 * @date 2023-10-08
 *
 * \details тест проводит газодинамический расчёт. Ожидаемое время выполнения 5минут
 * на выходе бинарные файлы. Требуется дополнительно вызвать rebuild_solve
 *
 * \note требуется поддержка VTK
 *
 */
#include "file_module.h"
#include "global_value.h"
#include "make_internal_format.h"

#include "hllc_main.h"
#include "solvers_struct.h"

#ifdef USE_VTK

int main(int argc, char **argv) {

  std::string prj_dir = "/home/artem/projects/solver/";

  glb_files.base_address = prj_dir + "tests/build/";
  glb_files.name_file_vtk = prj_dir + "tests/data/brick.vtk";
  glb_files.solve_address = glb_files.base_address;
  glb_files.Build();

  fs::create_directory(glb_files.base_address);
  fs::copy_file(glb_files.name_file_vtk, glb_files.base_address + "test.vtk", fs::copy_options::overwrite_existing);

  if (BuildDataFromVTK(glb_files) != e_completion_success) {
    return e_completion_fail;
  }

  _solve_mode.start_point = 0;
  _hllc_cfg.CFL = 0.7;
  _hllc_cfg.h_min = 1e-5;
  _hllc_cfg.save_timer = 0.1;
  _hllc_cfg.T = 0.5;

  hllc::RunHllcModule();

  return e_completion_success;
}

#endif