#include "solvers_struct.h"

#include "dbgdef.h"
#include "mpi_ext.h"
#include <omp.h>

solve_mode_t _solve_mode;
hllc_value_t _hllc_cfg;

Type TableFunc::operator()(Type x, Type y) {

  DIE_IF(data.size() == 0);

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
  try {
    return data[Ny * i + j];
  } catch (...) {
    WRITE_LOG_ERR("x=%lf, y=%lf, i=%d, j=%d, Nx=%d, Ny=%d\n", x, y, i, j, Nx, Ny);
    D_LD;
  }
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

#ifdef SPECTRUM
#include "plunk.h"
void grid_t::InitMemory(const IdType num_cells, const grid_directions_t &dir_grid) {

  DIE_IF(cells.size() != num_cells);

  // loc_size = hllc_loc_size[myid].right - hllc_loc_size[id].left;

  size_dir = dir_grid.size;
  size = num_cells;
  loc_size = size;
  loc_shift = 0;
  size_face = faces.size();

  get_splitting_spectrum(spectrum);
  size_frq = spectrum.size();

  inter_coef_all.resize(omp_get_max_threads());
  for (size_t i = 0; i < inter_coef_all.size(); i++) {
#ifndef ILLUM_ON_CELL
    inter_coef_all[i].resize(size * CELL_SIZE);
#else
    inter_coef_all[i].resize(size_face);
#endif

#ifdef SEPARATE_GPU
    loc_illum_wptr = 0;
    loc_illum_rptr = 0;
    local_Illum.resize(loc_illum_size * size, 0);
#endif
  }
}
#else
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
  }
}
#endif
grid_t::~grid_t() {
  inter_coef_all.clear();
#ifdef USE_MPI
  if (mpi_cfg) {
    delete mpi_cfg;
  }
#endif

#ifdef ILLUM
  for (auto &el : cells) {
    if (el.cell_data) {
      delete el.cell_data;
    }
  }
#endif
}
#endif // NOT USE_CUDA

#ifdef ILLUM
void grid_t::InitFullPhysData() {
  // for (auto &el : cells) {
  //   el.cell_data = new full_phys_data_t;
  // }
  for (size_t i = 0; i < size; i++) {
    cells[i].cell_data = new full_phys_data_t;
  }
}
#endif

#include "gas_state.h"
#include "global_value.h"
#ifdef ILLUM
#ifdef SPECTRUM
void full_phys_data_t::InitDirection(const Vector3 &dir) {
  cosf = 0;
  if (LIKELY(vel > kC_LightInv)) {
    cosf = val->v.dot(dir) / vel;
  }
}
#endif

void full_phys_data_t::Init(const flux_t *src) {

  val = src;
  T = GetTemperature(src->d, src->p);
  logT = log(T);

  Type vel2 = val->v.dot(val->v);
  vel = sqrt(vel2);

  lorenz = 1. / sqrt(1. - vel2);

  Type L = t_cooling_function(log(val->d) + LOG(kDensity), logT);
  Type log_alpha = L - (LOG(kStefanBoltzmann4) + 4 * logT) + LOG(kDist);
  alpha = exp(log_alpha);
  betta = (kSigma_thomson / kM_hydrogen * kDist) * val->d;
}
#endif