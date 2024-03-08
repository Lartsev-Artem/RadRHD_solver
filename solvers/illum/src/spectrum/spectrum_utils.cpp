#if defined SPECTRUM
#include "spectrum_utils.h"

#include "gas_state.h"
#include "global_value.h"

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
static std::array<Vector3, 3> VelArr = {Vector3(0, 0, 0), Vector3(0, 0.99, 0), Vector3(-0.99, 0, 0)};

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
#endif //! SPECTRUM