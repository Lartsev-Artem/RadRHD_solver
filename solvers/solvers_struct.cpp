#include "solvers_struct.h"

#include "dbgdef.h"
#include "mpi_ext.h"
#include <omp.h>

solve_mode_t _solve_mode;
hllc_value_t _hllc_cfg;

Type TableFunc::operator()(Type x, Type y) {
  if (x < min_x) {
    x = min_x;
  }

  if (y < min_y) {
    y = min_y;
  }

  if (x > max_x) {
    x = max_x;
  }

  if (y > max_y) {
    y = max_y;
  }

  int i = std::min((int)((round(x) - min_x) / step_x), Nx - 1);
  int j = std::min((int)((round(y) - min_y) / step_y), Ny - 1);
  return data[Ny * i + j];
}

flux_t flux_t::operator+(const flux_t &x) const {
  return flux_t(d + x.d, v + x.v, p + x.p);
}
void flux_t::operator=(const flux_t &x) {
  d = x.d;
  v = x.v;
  p = x.p;
}

void flux_t::operator+=(const flux_t &x) {
  d += x.d;
  v += x.v;
  p += x.p;
}
void flux_t::operator-=(const flux_t &x) {
  d -= x.d;
  v -= x.v;
  p -= x.p;
}
void flux_t::operator*=(const Type x) {
  d *= x;
  v *= x;
  p *= x;
}

void flux_t::operator/=(const Type x) {
  d /= x;
  v /= x;
  p /= x;
}

Type flux_t::operator[](const int i) const {
  return *((Type *)((uint8_t *)&(*this) + sizeof(Type) * i));
}
Type &flux_t::operator[](const int i) {
  return *((Type *)((uint8_t *)&(*this) + sizeof(Type) * i));
}

#ifndef USE_CUDA
void grid_t::InitMemory(const IdType num_cells, const IdType num_directions) {

  DIE_IF(cells.size() != num_cells);

  size = cells.size();

#if defined ILLUM
  Illum = new Type[num_directions * size * CELL_SIZE];
  scattering = new Type[num_directions * size];
  memset(Illum, 0.0, sizeof(Type) * num_directions * size * CELL_SIZE);
  memset(scattering, 0.0, sizeof(Type) * num_directions * size);

  inter_coef_all.resize(omp_get_max_threads());
  for (size_t i = 0; i < inter_coef_all.size(); i++) {
    inter_coef_all[i].resize(size * CELL_SIZE);
  }

  if (get_mpi_id() == 0) {
    for (IdType i = 0; i < size; i++) {
      cells[i].illum_val.illum.resize(num_directions * CELL_SIZE, 0);
    }
  }
#endif
}
#ifdef ILLUM
grid_t::~grid_t() {
  delete[] Illum;
  delete[] scattering;
  inter_coef_all.clear();
}
#endif
#else // CUDA
#include "cuda_interface.h"
void grid_t::InitMemory(const IdType num_cells, const grid_directions_t &dir_grid) {

  DIE_IF(cells.size() != num_cells);

  // loc_size = hllc_loc_size[myid].right - hllc_loc_size[id].left;

  size_dir = dir_grid.size;
  size = num_cells;
  loc_size = size;
  loc_shift = 0;
  size_face = faces.size();

  inter_coef_all.resize(omp_get_max_threads());
  for (size_t i = 0; i < inter_coef_all.size(); i++) {
#ifndef ILLUM_ON_CELL
    inter_coef_all[i].resize(size * CELL_SIZE);
#else
    inter_coef_all[i].resize(size_face);
#endif

#ifdef SEPARATE_GPU
    local_Illum.resize(dir_grid.loc_size * size, 0);
#endif
  }
}
grid_t::~grid_t() {
  inter_coef_all.clear();
}
#endif // NOT USE_CUDA
