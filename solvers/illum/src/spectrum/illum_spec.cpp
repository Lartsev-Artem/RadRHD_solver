#include "solvers_config.h"
#ifdef SPECTRUM
#include "illum_utils.h"

///@todo compile_time log_const

#include "illum_rad_func.h"

#include "compton.h"
#include "plunk.h"

static inline Type GetI(Type s, Type Q, Type S, Type I_0, Type k) {
  if (s * k > 1e-10) {
    Type src = (Q + S) / k;
    return (exp(-k * s) * (I_0 - src) + src);
  } else
    return (1 - s * k) * (I_0 + s * (Q + S));
}

Type GetIllum(const Vector3 &dir, const Vector3 &x,
              const Type s,
              const Type I_0,
              const Type int_scattering,
              const Type frq, const Type frq0,
              elem_t &cell, Type &k) {
  switch (_solve_mode.class_vtk) {

  case e_grid_cfg_default:
  case e_grid_cfg_static_illum:
  case e_grid_cfg_radiation: // test task
  case e_grid_cfg_full_init: // HLLC + Illum для конуса
  {
    Type S = int_scattering;
    Type d = cell.phys_val.d;
    Type p = cell.phys_val.p;
    Type T = illum::GetTemperature(d, p); // размерная

    Type T2 = T * T;
    Type T4 = T2 * T2;

    Type alpha = 0;
    //надо предрасчитывать логаримф руками - log(kRadiation);
    Type Q = B_Plank_log(T, frq, frq0);
    if (Q > -400.0) {
      Type L = t_cooling_function(log(d * kDensity), log(T));
      // Type alpha = exp(L) / (4 * PI * Q) * kDist;
      Type alpha_log = L - Q - log(4. * PI);

      alpha = exp(alpha_log + log(kDist));
      Q = exp(alpha_log + Q);
      Q *= (kDist / kRadiation);
    }

    Type v = cell.phys_val.v.norm();
    Type cosf = cell.phys_val.v.dot(dir) / v;
    Type betta = get_scat_coef(frq, v, cosf) / (kM_hydrogen * kDist) * d;

    cell.illum_val.absorp_coef = alpha;
    cell.illum_val.scat_coef = betta;

    return std::max(0.0, GetI(s, Q, S, I_0, alpha + betta));

  default:
    D_LD;
  }
  }
}
#endif //! SPECTRUM