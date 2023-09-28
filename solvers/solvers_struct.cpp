#ifdef SOLVERS
#include "solvers_struct.h"

#include "dbgdef.h"
#include "mpi_ext.h"

solve_mode_t _solve_mode;

flux_t flux_t::operator+(const flux_t &x) {
  d += x.d;
  v += x.v;
  p += x.p;
  return *this;
}
flux_t flux_t::operator+=(const flux_t &x) {
  d += x.d;
  v += x.v;
  p += x.p;
  return *this;
}
flux_t flux_t::operator-=(const flux_t &x) {
  d -= x.d;
  v -= x.v;
  p -= x.p;
  return *this;
}
flux_t flux_t::operator*(const Type x) {
  d *= x;
  v *= x;
  p *= x;
  return *this;
}
flux_t flux_t::operator-(const flux_t &x) {
  d -= x.d;
  v -= x.v;
  p -= x.p;
  return *this;
}
flux_t flux_t::operator/(const Type x) {
  d /= x;
  v /= x;
  p /= x;
  return *this;
}

// это временно для свзяи со старым кодом
Type flux_t::operator[](const int i) const {
  switch (i) {
  case 0:
    return d;
  case 1:
    return v[0];
  case 2:
    return v[1];
  case 3:
    return v[2];
  case 4:
    return p;

  default:
    D_LD;
  }
}
Type &flux_t::operator[](const int i) {
  switch (i) {
  case 0:
    return d;
  case 1:
    return v[0];
  case 2:
    return v[1];
  case 3:
    return v[2];
  case 4:
    return p;

  default:
    D_LD;
  }
}
Type flux_t::operator()(const int i) {
  switch (i) {
  case 0:
    return d;
  case 1:
    return v[0];
  case 2:
    return v[1];
  case 3:
    return v[2];
  case 4:
    return p;

  default:
    D_LD;
  }
}

illum_value_t::illum_value_t(const int num_dir) {

  absorp_coef = 0;
  rad_en_loose_rate = 0;

#ifndef USE_CUDA
  illum.resize(num_dir * CELL_SIZE, 0);
  energy = 0;
  stream = Vector3::Zero();
  impuls = Matrix3::Zero();
  div_stream = 0;
  div_impuls = Vector3::Zero();
#endif
}

#ifndef USE_CUDA
void grid_t::InitMemory(const uint32_t num_cells, const uint32_t num_directions) {

  DIE_IF(cells.size() != num_cells);

  size = cells.size();

#if defined ILLUM
  Illum = new Type[num_directions * size * CELL_SIZE];
  scattering = new Type[num_directions * size];

  if (get_mpi_id() == 0) {
    for (int i = 0; i < size; i++) {
      cells[i].illum_val.illum.resize(num_directions * CELL_SIZE, 0);
    }
  }
#endif
}
#ifdef ILLUM
grid_t::~grid_t() {
  delete[] Illum;
  delete[] scattering;
}
#endif
#endif // NOT USE_CUDA

#endif //! SOLVERS