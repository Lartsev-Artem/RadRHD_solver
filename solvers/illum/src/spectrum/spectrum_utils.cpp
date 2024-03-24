#if defined SPECTRUM
#include "spectrum_utils.h"

#include "gas_state.h"
#include "global_value.h"

#include "compton.h"
#include "illum_utils.h"
#include "plunk.h"

#include "reader_bin.h"

Type illum::spectrum::get_full_illum(const IdType num_dir, const grid_t &grid) {
  const IdType N = grid.size;
  const IdType M = grid.size_dir;
  Type sumI = 0;
  for (size_t cell = 0; cell < N; cell++) {
    sumI += grid.Illum[M * cell + num_dir];
  }
  return sumI * kRadiation / kDist;
}

#include <array>
static int vel_idx = 0;
static std::array<Vector3, 3> VelArr = {Vector3(0, 0, 0), Vector3(0.99, 0, 0), Vector3(-0.99, 0, 0)};

int illum::spectrum::InitPhysState(const int num, grid_t &grid) {

  std::vector<Type> density;
  std::vector<Type> pressure;
  std::vector<Vector3> velocity;

  uint32_t err = 0;
  err |= files_sys::bin::ReadSimple(glb_files.solve_address + std::to_string(num) + F_DENSITY, density);
  err |= files_sys::bin::ReadSimple(glb_files.solve_address + std::to_string(num) + F_PRESSURE, pressure);
  err |= files_sys::bin::ReadSimple(glb_files.solve_address + std::to_string(num) + F_VELOCITY, velocity);

  if (err) {
    WRITE_LOG_ERR("Error reading phys state\n");

    if (vel_idx >= VelArr.size()) {
      RETURN_ERR("No static state\n");
    }

#pragma omp parallel for
    for (int i = 0; i < grid.size; i++) {
      grid.cells[i].phys_val = flux_t(1e-11 / kDensity, VelArr[vel_idx], GetPressure(1e-11 / kDensity, 5000.0));
      grid.cells[i].cell_data->Init(&grid.cells[i].phys_val);
    }
    vel_idx++;
    WRITE_LOG("Init static phys state\n");
    return e_completion_success;
  }

  if (!((density.size() == pressure.size()) && (density.size() == velocity.size()) && (density.size() == grid.size))) {
    RETURN_ERR("Error data size in InitPhysState\n");
  }

#pragma omp parallel for
  for (int i = 0; i < grid.size; i++) {
    grid.cells[i].phys_val = flux_t(density[i], velocity[i], pressure[i]);
    grid.cells[i].cell_data->Init(&grid.cells[i].phys_val);
  }

  return e_completion_success;
}

Type illum::GetRhsOpt(const Vector3 x, const Type S, elem_t &cell, Type &k,
                      Type frq0, Type frq1) {

  full_phys_data_t *phys = cell.cell_data;

  Type betta;
  if (LIKELY(phys->vel > kC_LightInv)) {
    betta = (get_scat_coef(0.5 * (frq1 + frq0), phys->vel, phys->cosf, phys->lorenz) / (kM_hydrogen)) * (phys->val->d * kDensity) * kDist;
  } else {
    betta = (get_scat_coef(0.5 * (frq1 + frq0))) * (phys->val->d * kDensity / kM_hydrogen) * kDist;
  }

  Type alpha = 0; // phys->alpha;
  Type Q = B_Plank(phys->T, phys->logT, frq1, frq0) / kRadiation;
  Type SS = (phys->val->d * kDensity / kM_hydrogen) * S * kDist;

#ifdef LOG_SPECTRUM
  if (log_enable) {
    log_spectrum("S=%e, SS=%e, sig: %e %e %e\n",
                 S, SS, 0.5 * (frq1 + frq0), phys->vel, phys->lorenz);
  }

#endif

  k = alpha + betta;

#ifndef SAVE_FULL_SPECTRUM
  if (k < numeric_limit_abs_coef) {
    return (alpha * Q + SS);
  }
#endif

#ifdef DEBUG //
  Type res = (alpha * Q + SS) / k;
  if (res < 0 || std::isnan(res) || std::isinf(res)) {
    EXIT_ERR("res=%e %e %e %e \n", res, alpha, Q, SS);
  }
#endif

  return (alpha * Q + SS) / k;
}

#endif //! SPECTRUM