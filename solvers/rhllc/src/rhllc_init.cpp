#if defined RHLLC && defined SOLVERS
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "global_value.h"

#include "reader_bin.h"

Type Density(const Vector3 &p) {
  return 1e-8;
  const Type x = p[0];
  const Type x1 = 0.04;
  const Type x2 = 0.06;
  const Type x3 = 0.1;

  const Type a = 1e-9;
  const Type b = 1e-14;
  const Type betta0 = 0.012;
  const Type betta1 = 0.1;
  if (x < x1) {
    const Type val = (x1 - x2) / betta0;
    return a * exp(-val * val);
  }

  if (x < x2) {
    const Type val = (x - x2) / betta0;
    return a * exp(-val * val);
  }

  if (x < x3) {
    return a;
  }

  const Type val = (x - x3) / betta1;
  return a * exp(-val * val) + b;
}
Type Pressure(const Vector3 &p) {
  return 100;
  const Type x = p[0];
  const Type x1 = 0.04;
  const Type x2 = 0.06;
  const Type x3 = 0.1;

  const Type A = 200;
  const Type betta2 = 0.03;
  const Type betta3 = 0.02;
  const Type b = 1e-2;

  if (x < x1) {
    return A;
  }

  if (x < x2) {
    const Type val = (x - x1) / betta3;
    return A * exp(-val * val);
  }

  if (x < x3) {
    const Type val = (x2 - x1) / betta3;
    return A * exp(-val * val);
  }

  const Type val0 = (x2 - x1) / betta3;
  const Type coef = A * exp(-val0 * val0);

  const Type val = (x - x3) / betta2;
  return coef * exp(-val * val) + b;
}
Vector3 Velocity(const Vector3 &p) {
  return Vector3(kC_Light * 1e-5, 0, 0);
}

static int SetRHllcValueDefault(std::vector<elem_t> &cells) {
  std::vector<Vector3> centers;
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_CENTERS, centers))
    RETURN_ERR("Default rhllc value not set\n");

  int i = 0;
  for (auto &el : cells) {
#if GEOMETRY_TYPE == Cylinder
    Vector3 x = centers[i];
    if (Vector2(x[1], x[2]).norm() < 0.03 && x[0] < 0.1) {
      el.phys_val.d = 0.1;
      el.phys_val.p = 0.01;
      el.phys_val.v = Vector3(0.99, 0, 0);
    } else {
      el.phys_val.d = 10;
      el.phys_val.p = 0.01;
      el.phys_val.v = Vector3(0, 0, 0);
    }
#elif GEOMETRY_TYPE == Cube
    Type x = centers[i][0];
    if (x < 0.5) {
      el.phys_val.d = 1;
      el.phys_val.p = 1;
      el.phys_val.v = Vector3(0.9, 0, 0);
    } else {
      el.phys_val.d = 1;
      el.phys_val.p = 10;
      el.phys_val.v = Vector3(0, 0, 0);
    }
#elif GEOMETRY_TYPE == Cone_JET
    Vector3 x = centers[i];
    if (Vector2(x[1], x[2]).norm() < 0.01 && x[0] < 0.05) {
      el.phys_val.d = 0.1;
      el.phys_val.p = 0.01;
      el.phys_val.v = Vector3(0.99, 0, 0);
    } else {
      el.phys_val.d = 10;
      el.phys_val.p = 0.01;
      el.phys_val.v = Vector3(0, 0, 0);
    }

    // const Vector3 x = centers[i];
    // el.phys_val.d = Density(x) / DENSITY;
    // el.phys_val.p = Pressure(x) / PRESSURE;
    // el.phys_val.v = Velocity(x) / VELOCITY;
#elif GEOMETRY_TYPE == Cone

    {
      el.phys_val.d = 0.5;
      el.phys_val.p = 0.1;
      el.phys_val.v = Vector3(0, 0, 0);
    }

#elif GEOMETRY_TYPE == Sphere
    el.phys_val.d = 0.1;
    el.phys_val.p = 0.1;
    el.phys_val.v = Vector3(0, 0, 0);
#else
    el.phys_val.d = 0.1;
    el.phys_val.p = 1;
    el.phys_val.v = Vector3(1e-4, 0, 0);
#endif // Cube

    i++;
  } // for

  return e_completion_success;
}

// ------------------
void rhllc::SetCfgDefault(hllc_value_t &hllc_set) {

  hllc_set.h_min = 0.004; // default

#ifdef ILLUM
  hllc_set.h_min = 0.0007166575761593; // jet
  hllc_set.tau = 1e-7;
  hllc_set.CFL = 0.001;
  hllc_set.save_timer = 0.02;
  hllc_set.T = 1;
#else // ILLUM

  hllc_set.h_min = 0.0005369329546790;
  hllc_set.tau = 1e-5;
  hllc_set.CFL = 0.7;
  hllc_set.save_timer = 0.05; // hllc_set.h* hllc_set.CFL; //1e-5;
  hllc_set.T = 1;             // hllc_set.print_timer * 2 - 1e-6;

#endif // ILLUM
}

int rhllc::Init(std::string &file_init_value, std::vector<elem_t> &cells) {

  if (files_sys::bin::ReadHllcInit(file_init_value, cells) != e_completion_success) {
    if (SetRHllcValueDefault(cells) != e_completion_success) {
      WRITE_LOG("Default rhllc value not set\n");
      return e_completion_fail;
    }
    WRITE_LOG("SetDefault hllc value\n");
  }

  HllcPhysToConv(cells);

  _hllc_cfg.tau = GetTimeStep(_hllc_cfg, cells);
  return e_completion_success;
}

#endif //! RHLLC &&  SOLVERS
