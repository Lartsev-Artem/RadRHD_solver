#if defined SOLVERS && defined HLLC
#include "hllc_init.h"
#include "hllc_utils.h"

#include "global_value.h"
#include "reader_bin.h"

#include <omp.h>

void hllc::SetCfgDefault(hllc_value_t &hllc_set) {

  hllc_set.h_min = 0.005; // default
#ifdef ILLUM
  hllc_set.tau = 1e-8;
  hllc_set.CFL = 0.001;
  hllc_set.save_timer = 1e-5;
  hllc_set.T = 0.5;
#else // ILUM

  hllc_set.tau = 1e-5;
  hllc_set.CFL = 0.5;
  hllc_set.save_timer = 0.01;
  hllc_set.T = 0.4;

#endif // ILUM
}

static int SetHllcValueDefault(std::vector<elem_t> &cells) {

  std::vector<Vector3> centers;
  if (files_sys::bin::ReadSimple(glb_files.base_address + "centers.bin", centers))
    RETURN_ERR("Default hllc value not set\n");

#pragma omp parallel for
  for (int i = 0; i < cells.size(); i++) {
    elem_t &el = cells[i];
#if GEOMETRY_TYPE == Cone_JET
    const Type betta = 0.1;
    const Type a = 1;
    const Type b = 0.001;
    Type x = centers[i][0];
    el.phys_val.d = a * exp(-x * x / betta) + b;
    el.phys_val.p = a * exp(-x * x / betta) + (1e-5);
    el.phys_val.v = Vector3(1e-4, 0, 0);

#elif GEOMETRY_TYPE == Cube
    Type x = centers[i][0];
    if (x < 0.5) {
      el.phys_val.d = 1;
      el.phys_val.p = 1;
      el.phys_val.v = Vector3(0, 0, 0);
    } else {
      el.phys_val.d = 0.125;
      el.phys_val.p = 0.1;
      el.phys_val.v = Vector3(0, 0, 0);
    }
#elif GEOMETRY_TYPE == Sphere || GEOMETRY_TYPE == Step
    el.phys_val.d = 0.1;
    el.phys_val.p = 0.01;
    el.phys_val.v = Vector3(0, 0, 0);
#else
    el.phys_val.d = 0.1;
    el.phys_val.p = 0.01;
    el.phys_val.v = Vector3(0, 0, 0);
#endif
  }

  return e_completion_success;
}
int hllc::Init(std::string &file_init_value, std::vector<elem_t> &cells) {

  if (files_sys::bin::ReadHllcInit(file_init_value, cells) != e_completion_success) {
    if (SetHllcValueDefault(cells) != e_completion_success) {
      WRITE_LOG("Default rhllc value not set\n");
      return e_completion_fail;
    }
    WRITE_LOG("SetDefault hllc value\n");
  }

  HllcPhysToConv(cells);

  _hllc_cfg.tau = GetTimeStep(_hllc_cfg, cells);
  return e_completion_success;
}

#endif