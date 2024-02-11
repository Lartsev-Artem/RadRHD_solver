#include "solvers_config.h"
#ifdef SPECTRUM
#include "compton.h"
#include "global_value.h"
#include "illum_utils.h"
#include "plunk.h"

#include "spec_all.h"
#define SOURCE_TEMPERATURE 5000.0 ///< температура источника, заданного как АЧТ

static inline Type GetI(Type s, Type Q, Type S, Type I_0, Type k) {
  if (s * k > 1e-10) {
    Type src = (Q + S) / k;
    return (exp(-k * s) * (I_0 - src) + src);
  } else
    return (1.0 - s * k) * (I_0 + s * (Q + S));
}

#include <omp.h>
Type illum::spec::GetIllum(const Vector3 &dir, const Vector3 &x,
                           const Type s,
                           const Type I_0,
                           const Type int_scattering,
                           const Type frq0, const Type frq1,
                           elem_t &cell) {
  switch (_solve_mode.class_vtk) {

  case e_grid_cfg_default:
  case e_grid_cfg_static_illum:
  case e_grid_cfg_radiation: // test task
  case e_grid_cfg_full_init: // HLLC + Illum для конуса
  {
    Type S = int_scattering;
    Type d = 1;   // cell.phys_val.d;
    Type p = 0.1; // cell.phys_val.p;
    cell.phys_val.v = Vector3(0, 0, 0);
    Type T = SOURCE_TEMPERATURE; // GetTemperature(d, p); // размерная

    // Type T2 = T * T;
    // Type T4 = T2 * T2;

    Type alpha = 0;
    //надо предрасчитывать логаримф руками - log(kRadiation);
    Type Q = B_Plank_log(T, frq1, frq0);
    /// \warning тут что то не так с логарифмами
    // if (Q > -400.0) {
    //   Type L = 0; // t_cooling_function(log(d) + LOG(kDensity), log(T));
    //   Type alpha = exp(L) / (4 * PI * Q) * kDist;
    //   Type alpha_log = L - Q - LOG(PI4); // log(4. * PI);

    //   alpha = exp(alpha_log + LOG(kDist));
    //   Q = exp(alpha_log + Q);
    //   Q *= (kDist / kRadiation);
    // }
    alpha = 0;
    Q = alpha * exp(Q);

    Type v = cell.phys_val.v.norm();
    Type betta;
    if (v > 1e-10) {
      Type cosf = cell.phys_val.v.dot(dir) / v;
      betta = (get_scat_coef(0.5 * (frq1 + frq0), v, cosf) / (kM_hydrogen * kDist)) * d;
    } else {
      betta = (get_scat_coef(0.5 * (frq1 + frq0)) / (kM_hydrogen * kDist)) * d;
    }

    //    WRITE_LOG("Q=%lf, b=%lf\n", Q, betta);

    cell.illum_val.absorp_coef = alpha; //интегральный (сделать обнуление на старте)

    //здесь нужна интеграл и по частоте и по направлениям. Значит храним по локальному направлению(массив из числа потоков) и интегральный
    // cell.illum_val.scat_coef_loc[omp_get_num_threads()] += betta;

    return std::max(0.0, GetI(s, Q, betta * S, I_0, alpha + betta));
  }
  default:
    D_LD;
  }
}

Type illum::spec::BoundaryConditions(const IdType type_bound, const Type frq0, const Type frq1) {
  switch (type_bound) {

  case e_bound_free:
  case e_bound_lock:
  case e_bound_out_source:
    return 0;
  case e_bound_inner_source:
    return (B_Plank(SOURCE_TEMPERATURE, frq1, frq0)); // B_Plank(SOURCE_TEMPERATURE, frq1, frq0);

  default:
    D_LD;
  }
}

#ifndef SEPARATE_GPU
/// \todo spectrum: новая размерность в grid под частота
Type illum::spec::ReCalcIllum(const IdType num_dir, const IdType num_frq, const std::vector<Type> &inter_coef, grid_t &grid, IdType mpi_dir_shift) {

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
#else
Type illum::spec::ReCalcIllum(const IdType num_dir, const std::vector<std::vector<Type>> &inter_coef, grid_t &grid, const IdType dir_disp) {
  Type norm = -1;
  const IdType shift_dir = num_dir * grid.size * grid.size_frq;

  for (IdType frq = 0; frq < grid.size_frq; frq++) {

    for (IdType cell = 0; cell < grid.size; cell++) {

      Type curI = 0;
      for (size_t j = 0; j < CELL_SIZE; j++) {
        curI += inter_coef[grid.cells[cell].geo.id_faces[j]][frq]; //тут печаль с кэшами
      }
      curI /= CELL_SIZE;

      IdType id = (shift_dir + cell * grid.size_frq + frq); ///\todo это вопрос. см отправку
      norm = std::max(norm, fabs((grid.local_Illum[id] - curI) / curI));
      grid.local_Illum[id] = curI; //здесь по направлениям

      grid.Illum[cell * grid.size_frq * grid.size_dir + (grid.size_frq * (num_dir + dir_disp)) + frq] = curI; //здесь по ячейкам
    }
  }

  return norm;
}
#endif

#endif //! SPECTRUM