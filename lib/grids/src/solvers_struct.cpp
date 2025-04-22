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
    int i = std::min((int)((round(x) - min_x) / step_x), Nx - 1);

    Type y1 = exp(data[Ny * i + (Ny - 1)]);
    Type y0 = exp(data[Ny * i + (Ny - 2)]);

    Type x1 = exp(max_y);
    Type x0 = exp(max_y - step_y);

    Type k = (y0 - y1) / (x0 - x1);
    Type m = (-x1 * y0 + x0 * y1) / (x0 - x1);
    // WRITE_LOG("(%lf %lf) (%lf %lf), %lf %lf: %e (def=%lf)\n", x0, y0, x1, y1, x, y, log(k * exp(y) + m), data[Ny * i + Ny - 1]);
    return log(k * exp(y) + m);
    // y > max_y
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

void grid_t::InitMemory(const IdType num_cells, const grid_directions_t &dir_grid) {

  DIE_IF(cells.size() != num_cells);

  // loc_size = hllc_loc_size[myid].right - hllc_loc_size[id].left;

  size_dir = dir_grid.size;
  size = num_cells;
  loc_size = size;
  loc_shift = 0;
  size_face = faces.size();

  inter_coef_all.resize(omp_get_max_threads());

#ifndef SAVE_FULL_SPECTRUM
  for (size_t i = 0; i < inter_coef_all.size(); i++) {
#ifndef ILLUM_ON_CELL
    inter_coef_all[i].resize(size * CELL_SIZE);
#else
    inter_coef_all[i].resize(size_face);
#endif
  }
#else
  for (size_t i = 0; i < inter_coef_all.size(); i++) {
    inter_coef_all[i].resize(size_face);
    for (size_t j = 0; j < size_face; j++) {
      inter_coef_all[i][j].resize(size_frq);
    }
  }
#endif
}

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
  for (size_t i = 0; i < size; i++) {
    cells[i].cell_data = new full_phys_data_t;
  }
}
#endif

#ifdef SPECTRUM
#include "plunk.h"
void grid_t::InitFrq() {
  get_splitting_spectrum(frq_grid);
  size_frq = frq_grid.size() - 1;
  spectrum.resize(size_frq, 0);
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

void full_phys_data_t::Init(const flux_t *src, const elem_t *cell) {

  val = src;
  if (cell) {
    if (cell->geo.center[0] < 0.5) {
      T = 1e17;
    } else {
      T = 1000;
    }
  } else {
    T = GetTemperature(src->d, src->p);
  }

  logT = log(T);

  Type vel2 = val->v.dot(val->v);
  vel = sqrt(vel2);

  lorenz = 1. / sqrt(1. - vel2);

  Type L = t_cooling_function(log(val->d) + LOG(kDensity), logT);
  Type log_alpha = L - (LOG(kStefanBoltzmann4) + 4 * logT) + LOG(kDist);
  alpha = exp(log_alpha);
  betta = (kSigma_thomson / kM_hydrogen * kDist) * val->d * kDensity;
}
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
  betta = (kSigma_thomson / kM_hydrogen * kDist) * val->d * kDensity;
}

#include "reader_bin.h"
#include "reader_txt.h"
#include "trace_nodes.h"
int TracerData::Init(const global_files_t& files)
  {
    uint32_t err = 0;
    err |= files_sys::txt::ReadSimple(files.base_address + F_INTERNAL_BOUND, inter_boundary_face_id);
    err |= files_sys::txt::ReadInitBoundarySetInFaces(files.base_address + F_FACE_ID, inter_faces);


    err |= files_sys::bin::ReadSimple(files.base_address + F_NEIGHBOR, neighbours);
    err |= files_sys::bin::ReadSimple(files.base_address + F_TRACE_GRID, grid);
    err |= files_sys::bin::ReadSimple(files.base_address + F_TRACE_VERTEX, vertexs);
    err |= files_sys::bin::ReadNormals(files.base_address + F_NORMALS, normals);    

   { 
    err |= files_sys::bin::ReadGridGeo(files.name_file_geometry_cells, geo_grid.cells);
    err |= files_sys::bin::ReadGridGeo(files.name_file_geometry_faces, geo_grid.faces);

    DIE_IF(err);

    geo_grid.size_dir = 1;
    geo_grid.size = geo_grid.cells.size();
    geo_grid.loc_size = geo_grid.size;
    geo_grid.loc_shift = 0;
    geo_grid.size_face = geo_grid.faces.size();
    geo_grid.inter_coef_all.resize(1);
    geo_grid.inter_coef_all[0].resize(geo_grid.size_face);
    geo_grid.Illum = new Type[geo_grid.size];
  }
  
    vec_x.resize(vertexs.size());
    DIE_IF(trace::GetInterpolationNodes(vertexs, vec_x));
  
    graph.resize(normals.size(), 0);
    return e_completion_success;
  }
#endif