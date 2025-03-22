#if defined ILLUM && defined SOLVERS

#include "illum_utils.h"
#include "dbgdef.h"
#include "gas_state.h"
#include "global_value.h"
#include "plunk.h"

#include "compton.h"
#include "spectrum_utils.h"

TableFunc t_cooling_function;

#if GEOMETRY_TYPE == Cone
// static Type boundary_value = 1e52 / (PI * 0.1 * kDist * 0.1 * kDist) / kRadiation;
static Type boundary_value = 1e52 / (0.0001 * kDist * kDist) / kRadiation; // 0.0001 --- характерная площадь грани (зависит от сетки)
#elif GEOMETRY_TYPE == Sphere
static Type boundary_value = 1;
#else
static Type boundary_value = 0;
#endif

void illum::set_boundary_value(Type val) {
  boundary_value = val;
}

/**
 * @brief Точное решение уравнения переноса (вдоль луча) в ячейке
 *
 * @param s путь
 * @param Q источник
 * @param S интеграл рассеяния
 * @param I_0 начальное значение з.Коши
 * @param k коэффициент поглощения (alpha + betta)
 * @param betta коэффициент рассеяния
 * @return значение луча в конце отрезка интегрирования
 * @note осуществляется предельный переход для среды без поглощения
 */
static inline Type GetI(Type s, Type Q, Type S, Type I_0, Type alpha, Type betta) {
  Type k = alpha + betta;
  if (s * k > 1e-10) {
    Type src = (alpha * Q + betta * S) / k;
    return (exp(-k * s) * (I_0 - src) + src);
  } else
    return (1 - s * k) * (I_0 + s * (alpha * Q + S * betta));
}

Type illum::BoundaryConditions(const IdType type_bound, const IntId type_obj, const Vector3 &inter_coef) {
  switch (type_bound) {

  case e_bound_outer_surface:
  case e_bound_free:
  case e_bound_lock:
    return 0;
  case e_bound_out_source:
#if GEOMETRY_TYPE == Cone
    return boundary_value;
#endif
    return 0;

#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
  case e_bound_inner_source: {
    switch (type_obj) {
    case e_ray_intersect_none: // внутренняя граница не может быть не определённой
    case e_ray_intersect_rosh:
      D_LD;

    case e_ray_intersect_disk:
      return 0;

    case e_ray_intersect_sphere:
      return 20;

    default:
      return 0;
      Vector3 copy(inter_coef);
      return GetIllumeFromInFace(0, copy); //значение на ячейке (код не не границы)
    }
  }
#else
  case e_bound_inner_source:
#if GEOMETRY_TYPE == Sphere
    return boundary_value;
#endif
    return 0;
#endif
  default:
    D_LD;
  }
}

Type illum::GetIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, elem_t &cell) {
  switch (_solve_mode.class_vtk) {

  case e_grid_cfg_default:
  case e_grid_cfg_static_illum: {

#if GEOMETRY_TYPE == Sphere
    Type Q = 10;
    Type alpha = 10;
    Type betta = 5;
    Type S = int_scattering;

    if ((x - Vector3(0, 0, 0)).norm() > 0.3) // излучающий шар
    {
      Q = 0;
      alpha = 1;
      betta = 5;
    }
    return std::max(0.0, GetI(s, Q, S, I_0, alpha, betta));
#endif
#if GEOMETRY_TYPE == Cone
    Type Q = 0;
    Type alpha = 0.5;
    Type betta = 0.5;
    Type S = int_scattering;

    if (x[0] < 0.06) // излучающий слой
    {
      Q = 10;
      alpha = 1;
      betta = 2;
    }
    return std::max(0.0, GetI(s, Q, S, I_0, alpha, betta));
#endif

    D_LD;
  }

  case e_grid_cfg_radiation: // test task
  {
    Type S = int_scattering;
#if GEOMETRY_TYPE == TEST_ELLIPSE
    Type Q = cell.illum_val.rad_en_loose_rate; // Q=alpha*Ie
    Type alpha = cell.illum_val.absorp_coef;
    Type betta = alpha / 2; // просто из головы

    return std::max(0.0, GetI(s, Q, S, I_0, alpha, betta));
#elif GEOMETRY_TYPE == MAIN_ELLIPSE

    const Type ss = s * kDistAccretor;         // числа --- переход к размерным параметрам
    Type Q = cell.illum_val.rad_en_loose_rate; // Q=alpha*Ie
    Type alpha = cell.phys_val.d * cell.illum_val.absorp_coef;
    Type betta = cell.phys_val.d * (kSigma_thomson / kM_hydrogen);
    return std::max(0.0, GetI(s, Q, S, I_0, alpha, betta));
    // I = I_0 * exp(-alpha * ss) + Q * (1 - exp(-alpha * ss)) / alpha;
    // I = I_0 * (1 - alpha * ss) + Q * ss;
#else
    D_LD;
#endif
  }

#ifdef RAD_RHD
    ///\todo Здесь пересчёт в размерные единицы, потом вычисление излучения и возврат к безразмерным
  case e_grid_cfg_full_init: // HLLC + Illum для конуса
  {
    // переход к размерным параметрам
    Type S = int_scattering;
    Type d = cell.phys_val.d;
    Type p = cell.phys_val.p;
    Type T = GetTemperature(d, p); // размерная

    Type T2 = T * T;
    Type T4 = T2 * T2;

    Type L = t_cooling_function(log(d * kMass / (kDist * kDist * kDist)), log(T));
    Type alpha = exp(L) / (4 * kStefanBoltzmann * T4) * kDist;

    Type Q = B_Plank(T);

    Type betta = (kSigma_thomson / kM_hydrogen * kDist) * d;

    cell.illum_val.absorp_coef = alpha;
    cell.illum_val.scat_coef = betta;

    return std::max(0.0, GetI(s, Q, S, I_0, alpha, betta));
  }
#endif

  default:
    D_LD;
  }
}

#if 0
Type illum::GetIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, elem_t &cell) {

  Type f0, f1;
  Type S = int_scattering;
  Type d = cell.phys_val.d;
  Type p = cell.phys_val.p;
  Type T = GetTemperature(d, p); // размерная

#if GEOMETRY_TYPE == Sphere
  T = 27;
#endif

  Type T2 = T * T;
  Type T4 = T2 * T2;

  Type L = t_cooling_function(log(d * kMass / (kDist * kDist * kDist)), log(T));
  Type kSB = SigmaSB(T, f0, f1);
  Type alpha = exp(L) / (4 * kSB * T4) * kDist;

  Type Q = kSB * T4 / PI;

  Type betta = 0; // (kSigma_thomson / kM_hydrogen * kDist) * d;

  Type Q = 10;
  Type alpha = 10;
  Type betta = 5;
  Type S = int_scattering;

  if ((x - Vector3(0, 0, 0)).norm() > 0.3) // излучающий шар
  {
    Q = 0;
    alpha = 1;
    betta = 5;
  }
  return std::max(0.0, GetI(s, Q, S, I_0, alpha, betta));
}
#endif
static const BaseTetra_t tetra;
Type illum::ReCalcIllum(const IdType num_dir, const std::vector<Vector3> &inter_coef, grid_t &grid, IdType mpi_dir_shift) {

  Type norm = -1;
  const IdType shift_dir = num_dir * grid.size;

  for (IdType num_cell = 0; num_cell < grid.size; num_cell++) {

    for (IdType i = 0; i < CELL_SIZE; i++) {

      Vector3 Il = inter_coef[num_cell * CELL_SIZE + i];

      const Type curI = (Il[0] + Il[1] + Il[2]) / 3; //  среднее на грани (в идеале переход к ax+by+c)
      const IdType id = mpi_dir_shift + CELL_SIZE * (shift_dir + num_cell) + i;

      // if (curI < 1e-15) // защита от деления на ноль
      //        norm = 1;
      norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));

      grid.Illum[id] = curI;

#ifndef USE_CUDA
      grid.cells[num_cell].illum_val.illum[num_dir * CELL_SIZE + i] = curI; //на каждой грани по направлениям
#endif
    }
  }

  return norm;
}

Type illum::GetIllumeFromInFace(const IdType neigh_id, Vector3 &inter_coef
#ifdef INTERPOLATION_ON_FACES
                                ,
                                const Vector2 &x0
#endif
) {

#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
  if (neigh_id != e_bound_inner_source) //при использовании трассировки сквозь границу, внутренняя грань определена до этого
#endif
  {
    if (neigh_id < 0) {
      Type I_x0 = illum::BoundaryConditions(neigh_id);
      inter_coef = Vector3(I_x0, I_x0, I_x0);
      return I_x0;
    }
  }

  Vector3 &coef = inter_coef;

#ifdef INTERPOLATION_ON_FACES

  //#pragma error("Unsupported config")

  Type I_x0 = x0[0] * coef[0] + x0[1] * coef[1] + coef[2];

#else
  /// \note сейчас храним значения а не коэффициента интерполяции
  Type I_x0 = (coef[0] + coef[1] + coef[2]) / 3.;
#endif

  if (I_x0 < 0) {
    D_L;
    return 0;
  }
  return I_x0;
}

#ifdef TRANSFER_CELL_TO_FACE
Type illum::ReCalcIllum(const IdType num_dir, const std::vector<Type> &inter_coef, grid_t &grid, IdType mpi_dir_shift) {

  Type norm = -1;
  const IdType shift_dir = num_dir * grid.size;

#ifdef ILLUM_ON_CELL //на ячейках
  for (size_t cell = 0; cell < grid.size; cell++) {
    Type curI = 0;
    for (size_t j = 0; j < CELL_SIZE; j++) {
      curI += inter_coef[grid.cells[cell].geo.id_faces[j]];
    }
    curI /= CELL_SIZE;

    IdType id = mpi_dir_shift + (shift_dir + cell);
    norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));
    grid.Illum[id] = curI;
  }
#else //на гранях
  for (size_t cell = 0; cell < grid.size; cell++) {
    for (size_t j = 0; j < CELL_SIZE; j++) {
      Type curI = inter_coef[grid.cells[cell].geo.id_faces[j]];

      IdType id = mpi_dir_shift + CELL_SIZE * (shift_dir + cell) + j;
      norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));
      grid.Illum[id] = curI;
    }
  }
#endif
  return norm;
}
#endif

#if defined SEPARATE_GPU
Type illum::separate_gpu::ReCalcIllum(const IdType num_dir, const std::vector<Type> &inter_coef, grid_t &grid, const IdType dir_disp) {
  Type norm = -1;
  const IdType shift_dir = num_dir + dir_disp;

  for (size_t cell = 0; cell < grid.size; cell++) {
    Type curI = 0;
    for (size_t j = 0; j < CELL_SIZE; j++) {
      curI += inter_coef[grid.cells[cell].geo.id_faces[j]];
    }
    curI /= CELL_SIZE;

    IdType id = cell * grid.size_dir + shift_dir;
    norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));
    grid.Illum[id] = curI;
  }

  return norm;
}

// 0-чтение, 0-сразу выкинуть из кэша
#define PREFETCH_FACE_COEF(_cell_id) //                                         \
  __builtin_prefetch(&inter_coef[grid.cells[_cell_id].geo.id_faces[0]], 0, 0); \
  __builtin_prefetch(&inter_coef[grid.cells[_cell_id].geo.id_faces[1]], 0, 0); \
  __builtin_prefetch(&inter_coef[grid.cells[_cell_id].geo.id_faces[2]], 0, 0); \
  __builtin_prefetch(&inter_coef[grid.cells[_cell_id].geo.id_faces[3]], 0, 0);

Type illum::separate_gpu::ReCalcIllumOpt(const IdType num_dir, const std::vector<Type> &inter_coef, grid_t &grid, const IdType dir_disp) {
  Type norm = -1;
  const IdType shift_dir = num_dir + dir_disp;

  PREFETCH_FACE_COEF(0);

  for (size_t cell = 0; cell < grid.size - 1; cell++) {

    PREFETCH_FACE_COEF(cell + 1);
    const IdType id = cell * grid.size_dir + shift_dir;
    //__builtin_prefetch(&grid.Illum[id], 1, 0);

#if 0 
    Type curI = 0;
#pragma loop unroll(CELL_SIZE)
    for (size_t j = 0; j < CELL_SIZE; j++) {
      curI += inter_coef[grid.cells[cell].geo.id_faces[j]];
    }
    curI /= CELL_SIZE;
#else
    Type curI = 0;
    Type S = 0;
#pragma loop unroll(CELL_SIZE)
    for (size_t j = 0; j < CELL_SIZE; j++) {
      curI += inter_coef[grid.cells[cell].geo.id_faces[j]] *
              grid.faces[grid.cells[cell].geo.id_faces[j]].geo.S;
      S += grid.faces[grid.cells[cell].geo.id_faces[j]].geo.S;
    }
    curI /= S;
#endif

    norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));
    grid.Illum[id] = curI;
  }

  //последняя итерация
  {
    const size_t cell = grid.size - 1;
    const IdType id = cell * grid.size_dir + shift_dir;
    //__builtin_prefetch(&grid.Illum[id], 1, 0);
    Type curI = 0;
#pragma loop unroll(CELL_SIZE)
    for (size_t j = 0; j < CELL_SIZE; j++) {
      curI += inter_coef[grid.cells[cell].geo.id_faces[j]];
    }
    curI /= CELL_SIZE;

    norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));
    grid.Illum[id] = curI;
  }

  return norm;
}
#endif

#pragma GCC target("avx2")
#pragma GCC optimize("O3")

#include <bits/stdc++.h>
#include <x86intrin.h>
static inline __m256d _mm256_exp_pd(__m256d x) {
  alignas(32) Type X[4];
  _mm256_store_pd(X, x);
  return _mm256_setr_pd(
      exp(X[0]),
      exp(X[1]),
      exp(X[2]),
      0);
}

#ifdef TRANSFER_CELL_TO_FACE
Type illum::GetIllum(const Type *I0, const Type *s, const Type k, const Type rhs) {

  // alignas(32) Type Icur[4];
#if 0 // ndef LOG_SPECTRUM
  __m256d RHS = _mm256_set1_pd(rhs);
  __m256d S = _mm256_load_pd(s);
  __m256d I_0 = _mm256_load_pd(I0);

  __m256d EXP = _mm256_exp_pd(_mm256_mul_pd(_mm256_set1_pd(-k), S));
  __m256d I = _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(I_0, RHS), EXP), RHS);
  _mm256_store_pd(Icur, _mm256_max_pd(_mm256_set1_pd(0.0), I));

#elif 0
  for (size_t i = 0; i < 3; i++) {
    Icur[i] = exp(-k * s[i]) * (I0[i] - rhs) + rhs;
    if (log_enable) {
      log_spectrum("k=%e,s=%e, I0=%e, rhs %e, I=%e\n", k, s[i], I0[i], rhs, Icur[i]);
    }
  }

  if (log_enable) {

    log_spectrum("Ires=%e\n", (Icur[0] + Icur[1] + Icur[2]) / 3.);
    log_enable = 0;
  }
#else
  Type Icur = exp(-k * s[0]) * (I0[0] - rhs) + rhs;
  Icur += (exp(-k * s[1]) * (I0[1] - rhs) + rhs);
  Icur += (exp(-k * s[2]) * (I0[2] - rhs) + rhs);
  return (Icur * 0.3333333333333333);
#endif
  // return (Icur[0] + Icur[1] + Icur[2]) / 3.;
}

Type illum::GetIllumLimit(const Type *I0, const Type *s, const Type k, const Type rhs) {
  // (1 - s * k) * (I_0 + s * (alpha * Q + S * betta));

  alignas(32) Type Icur[4];
  for (size_t i = 0; i < 3; i++) {
    Icur[i] = (1 - s[i] * k) * (I0[i] - rhs) + rhs;

#ifdef LOG_SPECTRUM
    if (log_enable) {
      log_spectrum("Lim k=%e s=%e, I0=%e, rhs %e, I=%e\n", k, s[i], I0[i], rhs, Icur[i]);
    }
#endif
  }

  // __m256d RHS = _mm256_set1_pd(rhs);
  // __m256d S = _mm256_load_pd(s);
  // __m256d I_0 = _mm256_load_pd(I0);
  // __m256d K = _mm256_set1_pd(k);

  // __m256d sk = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(S, _mm256_set1_pd(k)));
  // __m256d I = _mm256_mul_pd(sk, _mm256_add_pd(I_0, _mm256_mul_pd(S, RHS)));
  // _mm256_store_pd(Icur, _mm256_max_pd(_mm256_set1_pd(0.0), I));

#ifdef LOG_SPECTRUM
  if (log_enable) {
    log_spectrum("Lim Ires=%e\n", (Icur[0] + Icur[1] + Icur[2]) / 3.);
    log_enable = 0;
  }
#endif
  return (Icur[0] + Icur[1] + Icur[2]) / 3.;
}
#pragma GCC pop("O3");

Type illum::GetRhs(const Vector3 x, const Type int_scattering, elem_t &cell, Type &k) {
  switch (_solve_mode.class_vtk) {

  case e_grid_cfg_default:
  case e_grid_cfg_static_illum: {

#if GEOMETRY_TYPE == Sphere
    Type Q = 10;
    Type alpha = 10;
    Type betta = 5;
    const Type S = int_scattering;
    if ((x - Vector3::Zero()).norm() > 0.3) // излучающий шар
    {
      Q = 0;
      alpha = 1;
      betta = 5;
    }
    k = alpha + betta;
    cell.illum_val.absorp_coef = alpha;
    cell.illum_val.scat_coef = betta;

    return (alpha * Q + betta * S) / k;
#endif
#if GEOMETRY_TYPE == Cone
    Type Q = 0;
    Type alpha = 0.5;
    Type betta = 0.5;
    Type S = int_scattering;

    if (x[0] < 0.06) // излучающий слой
    {
      Q = 10;
      alpha = 1;
      betta = 2;
    }
    k = alpha + betta;
    return (alpha * Q + betta * S) / k;
#endif

    D_LD;
  }

  case e_grid_cfg_radiation: // test task
  {
    Type S = int_scattering;
#if GEOMETRY_TYPE == TEST_ELLIPSE
    Type Q = cell.illum_val.rad_en_loose_rate; // Q=alpha*Ie
    Type alpha = cell.illum_val.absorp_coef;
    Type betta = alpha / 2; // просто из головы

    k = alpha + betta;
    return (alpha * Q + betta * S) / k;
#elif GEOMETRY_TYPE == MAIN_ELLIPSE

    // const Type ss = s * kDistAccretor;         // числа --- переход к размерным параметрам
    Type Q = cell.illum_val.rad_en_loose_rate; // Q=alpha*Ie
    Type alpha = cell.phys_val.d * cell.illum_val.absorp_coef;
    Type betta = cell.phys_val.d * (kSigma_thomson / kM_hydrogen);
    k = alpha + betta;
    return (alpha * Q + betta * S) / k;
#else
    D_LD;
#endif
  }

#ifdef RAD_RHD
    ///\todo Здесь пересчёт в размерные единицы, потом вычисление излучения и возврат к безразмерным
  case e_grid_cfg_full_init: // HLLC + Illum для конуса
  {
    // переход к размерным параметрам
    Type S = int_scattering;
    Type d = cell.phys_val.d;
    Type p = cell.phys_val.p;
    Type T = GetTemperature(p, d); // размерная

    Type T2 = T * T;
    Type T4 = T2 * T2;

    Type L = t_cooling_function(log(d * kMass / (kDist * kDist * kDist)), log(T));
    Type alpha = exp(L) / (4 * kStefanBoltzmann * T4) * kDist;

    Type Ie = B_Plank(T);
    Type Q = alpha * Ie;

    Type betta = (kSigma_thomson / kM_hydrogen * kDist) * d;

    cell.illum_val.absorp_coef = alpha;
    cell.illum_val.scat_coef = betta;

    k = alpha + betta;
    return (alpha * Q + betta * S) / k;
  }
#endif

  default:
    D_LD;
  }
}

Type illum::GetRhsOpt(const Vector3 x, const Type int_scattering, elem_t &cell, Type &k) {

  // переход к размерным параметрам
  const Type S = int_scattering;
  const Type alpha = cell.cell_data->alpha;
  const Type betta = cell.cell_data->betta;

  // cell.illum_val.absorp_coef = alpha;
  // cell.illum_val.scat_coef = betta;

  Type Q = B_Plank(cell.cell_data->T) / kRadiation;

  k = alpha + betta;
  return (alpha * Q + betta * S) / k;
}

#endif // TRANSFER_CELL_TO_FACE

#endif //! defined ILLUM && defined SOLVERS