#if defined SOLVERS && defined RHLLC
#include "rhllc_flux_stab.h"
#include "rhllc_utils.h"

#include "global_value.h"
#include "linear_alg.h"

#define MAX_ITER 40

void rhllc::GetConvValueStab(flux_t &W, flux_t &U) {

  const Type h = 1 + kGamma_g * W.p / W.d;

  Type g = W.v.dot(W.v);
  if (g >= 1.0) {
    WRITE_LOG("Warning fix mod_vel=%lf\n", g);
    constexpr Type beta_fix = 0.9999;
    g = beta_fix / sqrt(g);
    W.v *= g;

    g = beta_fix * beta_fix;
  }

  g = 1.0 / (1.0 - g);
  Type scrh = W.d * h * g;
  g = sqrt(g);

  U.d = W.d * g;
  U.v = scrh * W.v;
  U.p = scrh - W.p;
}

int rhllc::PhysPressureFix(flux_t &U, flux_t &W) {

  /* ----------------------------------------------
     1. Solve f(u) = 0 with secant method
     ---------------------------------------------- */

  Type p = kMinPressure;
  Type D = U.d;
  Type m2 = U.v.dot(U.v);
  Type m = sqrt(m2);
  Type umax = m / D;
  Type u0 = umax;

  Type lor = sqrt(1.0 + u0 * u0);
  Type plor = p * lor;
  Type Dh = D + plor * kGamma_g;

  Type f0 = m / Dh - u0;

  Type u1 = (-D + sqrt(D * D + 4.0 * m * kGamma_g * p)) / (2.0 * kGamma_g * p);
  int done = 0;

  Type f1, du;
  for (int k = 1; k < MAX_ITER; k++) {
    lor = sqrt(1.0 + u1 * u1);
    plor = p * lor;
    Dh = D + plor * kGamma_g;

    f1 = m / Dh - u1;

    if (done == 1)
      break;
    du = (u1 - u0) / (f1 - f0) * f1;
    u0 = u1;
    f0 = f1;

    u1 -= du;
    u1 = std::min(u1, umax);
    u1 = std::max(u1, 0.0);
    if (fabs(f1) < 1.e-9)
      done = 1;
  }

  if (done != 1) // k>MAX_ITER
    return 2;

  /* ----------------------------------------------
     2. Solution u has been found.
        Update to converged value of u.
     ---------------------------------------------- */

  lor = sqrt(1.0 + u1 * u1);
  plor = p * lor;
  Dh = D + plor * kGamma_g;

  W.d = U.d / lor;
  W.p = p;

  U.p = Dh * lor - p; /* Redefine energy */

  f0 = 1.0 / (Dh * lor); /* = 1 / W */
  W.v = U.v * f0;

  return e_completion_success;
}

int rhllc::GetPhysValueStab(const flux_t &U, flux_t &W) {

  constexpr double tol = 1.e-11;
  constexpr double eps2 = 1.e-12; /* Maximum 1/gamma^2 */

  Type D = U.d;
  Type E = U.p;
  Type m2 = U.v.dot(U.v);

  Type Q = E - sqrt(m2 + D * D);

  if (Q < 0.0)
    return 1; /* Equation does not admit a solution */

  Type m = sqrt(m2);
  Type p = std::max(m - E, 0.0);
  Type D_1 = 1.0 / D;

  double pmin = sqrt(m2 / (1.0 - eps2)) - E;

  double alpha, alpha2, lor2, lor, tau, theta;
  double h, dh_dp, dh_dtau, yp, dp, dyp;

  /* ----------------------------------------------
     1. Solve f(p) = 0 by Newton's method
     ---------------------------------------------- */

  p = std::max(p, pmin);
  int iter;
  for (iter = 0; iter < MAX_ITER; iter++) {

    alpha = E + p;
    alpha2 = alpha * alpha;
    lor2 = 1.0 - m2 / alpha2;
    lor2 = 1.0 / lor2;

    DIE_IF(lor2 < 1.0);
    lor = sqrt(lor2);

    tau = lor * D_1;
    theta = p * tau;

    h = 1.0 + kGamma_g * theta;
    dh_dp = kGamma_g * tau;
    dh_dtau = kGamma_g * p;

    yp = D * h * lor - E - p;
    dyp = D * lor * dh_dp - m2 * lor2 * lor / (alpha2 * alpha) * (lor * dh_dtau + D * h) - 1.0;
    dp = yp / dyp;
    p -= dp;

    if (p < pmin)
      p = pmin;
    if (fabs(dp) < tol * E)
      break;
  }

  /* ----------------------------------------------
     1b. Check if solution is consistent
     ---------------------------------------------- */

  if (p < 0.0)
    return 2;
  if (std::isnan(p))
    return 4;
  if (iter >= MAX_ITER || fabs(yp / (E + p)) > 1.e-4 || p < (m - E)) {
    WRITE_LOG("GetPhysValueStab: solution may be inaccurate %d\n", iter);
    return 3;
  }

  /* ----------------------------------------------
     2. Solution has been found, update to
        converged value of p.
     ---------------------------------------------- */

  alpha = E + p;
  alpha2 = alpha * alpha;
  lor2 = alpha2 / (alpha2 - m2);
  lor = sqrt(lor2);
  tau = lor * D_1;

  // h = 1.0 + kGamma_g * p * tau;

  W.d = 1.0 / tau;
  W.p = p;
  W.v = U.v / (U.p + p); /* = 1 / W */

  return e_completion_success;
}

/**
 * @brief Get the Signal Speed
 *
 * @param cs2 квадрат скорости звука
 * @param Vn нормальная компонента скорости
 * @param V_norm квадрат нормы скорости
 * @param cmin
 * @param cmax
 */
static void inline GetSignalSpeed(Type cs2, Type Vn, Type V2_norm, Type &cmin, Type &cmax) {
  const Type Vn2 = Vn * Vn;
  const Type Vt2 = V2_norm - Vn2;
  const Type sroot_L = sqrt(cs2 * (1 - Vn2 - Vt2 * cs2) * (1 - V2_norm));
  cmax = (Vn * (1.0 - cs2) + sroot_L) / (1.0 - V2_norm * cs2);
  cmin = (Vn * (1.0 - cs2) - sroot_L) / (1.0 - V2_norm * cs2);
}

Type rhllc::GetFluxStab(const flux_t &conv_val_l, const flux_t &conv_val_r,
                        const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f) {
  Matrix3 T;
  GetRotationMatrix(f.geo.n, T);

  flux_t U_L = conv_val_l;
  U_L.v = T * conv_val_l.v;

  flux_t U_R = conv_val_r;
  U_R.v = T * conv_val_r.v;

  flux_t W_L = phys_val_l;
  W_L.v = T * phys_val_l.v;

  flux_t W_R = phys_val_r;
  W_R.v = T * phys_val_r.v;

  //==================== Кэшируем физические переменные слева и справа============================//
  // нормальная скорость
  const Vector3 &Vel_L = W_L.v; // T * velocity[num_cell];
  const Vector3 &Vel_R = W_R.v; // T * velocity[neig / 4];

  const Type d_L = W_L.d;
  const Type d_R = W_R.d;

  const Type p_L = W_L.p;
  const Type p_R = W_R.p;

  const Type VV_L = Vel_L.dot(Vel_L);
  const Type VV_R = Vel_R.dot(Vel_R);

  //=========================Вычисляем релятивистикие параметры============================//
  const Type g_L = 1. / sqrt(1 - VV_L); // фактор Лоренца
  const Type g_R = 1. / sqrt(1 - VV_R);

  const Type h_L = 1 + kGamma_g * p_L / d_L; // энтальпия
  const Type h_R = 1 + kGamma_g * p_R / d_R;

  const Type cs_L2 = ((kGamma1 * p_L) / (d_L * h_L)); //квадрат скорости звука
  const Type cs_R2 = ((kGamma1 * p_R) / (d_R * h_R));

  Type cmax_L, cmax_R;
  Type cmin_L, cmin_R;
  GetSignalSpeed(cs_L2, Vel_L[0], VV_L, cmin_L, cmax_L);
  GetSignalSpeed(cs_R2, Vel_R[0], VV_R, cmin_R, cmax_R);

  const Type lambda_L = std::min(cmin_L, cmin_R);
  const Type lambda_R = std::max(cmax_L, cmax_R);

  /* --------------------------------------------------
          compute HLLC  flux
     -------------------------------------------------- */
  flux_t F;
  if (lambda_L >= 0.0) {
    F.d = U_L.d * Vel_L[0]; // D*v_x
    F.v = U_L.v * Vel_L[0]; // mx*vx+p
    F.p = U_L.v[0];

    F.v[0] += p_L;

    // WRITE_LOG("L: %lf %lf %lf %lf %lf\n", F.d, F.v[0], F.v[1], F.v[2], F.p);
    //  sweep->press[i] = pL[i];

  } else if (lambda_R <= 0.0) {

    F.d = U_R.d * Vel_R[0]; // D*v_x
    F.v = U_R.v * Vel_R[0]; // mx*vx+p
    F.p = U_R.v[0];
    F.v[0] += p_R;

    // WRITE_LOG("R: %lf %lf %lf %lf %lf\n", F.d, F.v[0], F.v[1], F.v[2], F.p);

  } else {

#if 0 // HLL
    flux_t F_L;               //(U_L.d * Vel_L[0], U_L.v * Vel_L[0], U_L.v[0]);
    flux_t F_R;               //(U_R.d * Vel_R[0], U_R.v * Vel_R[0], U_R.v[0]);
    F_L.d = U_L.d * Vel_L[0]; // D*v_x
    F_L.v = U_L.v * Vel_L[0]; // mx*vx+p
    F_L.p = U_L.v[0];

    F_R.d = U_R.d * Vel_R[0]; // D*v_x
    F_R.v = U_R.v * Vel_R[0]; // mx*vx+p
    F_R.p = U_R.v[0];

    Type scrh = 1.0 / (lambda_R - lambda_L);

    for (int i = 0; i < 5; i++) {
      F[i] = lambda_L * lambda_R * (U_R[i] - U_L[i]) + lambda_R * F_L[i] - lambda_L * F_R[i];
      F[i] *= scrh;
    }

    F.v[0] += (lambda_R * p_L - lambda_L * p_R) * scrh;
    // sweep->press[i] = (SR[i] * pL[i] - SL[i] * pR[i]) * scrh;
#else

#if 0 // SHOCK_FLATTENING == MULTID
// scrh = MAX(fabs(lambda_L), fabs(lambda_R));  
    if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i + 1] & FLAG_HLL)) {
      scrh = 1.0 / (SR[i] - SL[i]);
      NFLX_LOOP(nv) {
        sweep->flux[i][nv] = SL[i] * SR[i] * (uR[nv] - uL[nv]) + SR[i] * fL[i][nv] - SL[i] * fR[i][nv];
        sweep->flux[i][nv] *= scrh;
      }
      sweep->press[i] = (SR[i] * pL[i] - SL[i] * pR[i]) * scrh;
      continue;
    }
#endif

    /* ---------------------------------------
                       get u*
       --------------------------------------- */
    flux_t F_L(U_L.d * Vel_L[0], U_L.v * Vel_L[0], U_L.v[0]);
    flux_t F_R(U_R.d * Vel_R[0], U_R.v * Vel_R[0], U_R.v[0]);
    Type AL = lambda_L * U_L.p - F_L.p;
    Type AR = lambda_R * U_R.p - F_R.p;

    Type BL = lambda_L * U_L.v[0] - F_L.v[0] - p_L;
    Type BR = lambda_R * U_R.v[0] - F_R.v[0] - p_R;

    Type a = AR * lambda_L - AL * lambda_R;
    Type b = AL + BL * lambda_R - AR - BR * lambda_L;
    Type c = BR - BL;
    /*
          if (fabs(a) > 1.e-9){
            usp = 0.5*(- b + sqrt(b*b - 4.0*a*c))/a;
            usm = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a;
          }else{
            usp = usm = -c/b;
          }
    */
    Type scrh = -0.5 * (b + SIGN(b) * sqrt(b * b - 4.0 * a * c));
    Type us = c / scrh;

    Type ps = (AL * us - BL) / (1.0 - us * lambda_L);

    if (us >= 0.0) {
      flux_t usl;
      Type dif = 1.0 / (lambda_L - us);

      usl.d = U_L.d * (lambda_L - Vel_L[0]) * dif;
      usl.v[0] = (lambda_L * (U_L.p + ps) - U_L.v[0]) * us * dif;
      usl.v[1] = U_L.v[1] * (lambda_L - Vel_L[0]) * dif;
      usl.v[2] = U_L.v[2] * (lambda_L - Vel_L[0]) * dif;
      usl.p = U_L.p + (usl.v[0] - U_L.v[0]) / lambda_L;

      // F = F_L + lambda_L *(usl - U_L);
      usl -= U_L;
      usl *= lambda_L;
      F = F_L + usl;

      F.v[0] += p_L;

      // WRITE_LOG("UL: %lf %lf %lf %lf %lf\n", F.d, F.v[0], F.v[1], F.v[2], F.p);
    } else {
      flux_t usr;
      Type dif = 1.0 / (lambda_R - us);
      usr.d = U_R.d * (lambda_R - Vel_R[0]) * dif;
      usr.v[0] = (lambda_R * (U_R.p + ps) - U_R.v[0]) * us * dif;
      usr.v[1] = U_R.v[1] * (lambda_R - Vel_R[0]) * dif;
      usr.v[2] = U_R.v[2] * (lambda_R - Vel_R[0]) * dif;
      usr.p = U_R.p + (usr.v[0] - U_R.v[0]) / lambda_R;

      usr -= U_R;
      usr *= lambda_R;
      F = F_R + usr;
      F.v[0] += p_R;

      // WRITE_LOG("UR: %lf %lf %lf %lf %lf\n", F.d, F.v[0], F.v[1], F.v[2], F.p);
    }
#endif // HLL
  }
  //  WRITE_LOG("F: %lf, %lf, %lf, %lf, %lf\n", F.d, F.v[0], F.v[1], F.v[2], F.p);

  f.f = F;
  f.f.v = (T.transpose()) * F.v;
  f.f *= f.geo.S;

  return std::max(lambda_L, lambda_R);
}

#pragma GCC target("avx2")
#pragma GCC optimize("O3")

#include <bits/stdc++.h>
#include <x86intrin.h>

using namespace std;

Type rhllc::GetFluxStabVec(const flux_t &conv_val_l, const flux_t &conv_val_r,
                           const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f) {
  Matrix3 T;
  GetRotationMatrix(f.geo.n, T);

  alignas(16) Type U_L[5];
  alignas(16) Type U_R[5];
  alignas(16) Type F[5];

  const Vector3 &Uvl = T * conv_val_l.v;
  U_L[0] = conv_val_l[0];
  U_L[1] = Uvl[0];
  U_L[2] = Uvl[1];
  U_L[3] = Uvl[2];
  U_L[4] = conv_val_l[4];

  const Vector3 &Uvr = T * conv_val_r.v;
  U_R[0] = conv_val_r[0];
  U_R[1] = Uvr[0];
  U_R[2] = Uvr[1];
  U_R[3] = Uvr[2];
  U_R[4] = conv_val_r[4];

  const Vector3 &Wvl = T * phys_val_l.v;
  const Vector3 &Wvr = T * phys_val_r.v;

  alignas(16) Type W_K[10] =
      {
          phys_val_l[0],
          phys_val_r[0],

          Wvl[0],
          Wvr[0],

          Wvl[1],
          Wvr[1],

          Wvl[2],
          Wvr[2],

          phys_val_l[4],
          phys_val_r[4],
      };

  //==================== Кэшируем физические переменные слева и справа============================//

  const __m128d ones2 = _mm_set1_pd(1.0);
  __m128d dk = _mm_load_pd(W_K);
  __m128d vnk = _mm_load_pd(W_K + 2);
  __m128d vk = _mm_setr_pd(Wvl.dot(Wvl), Wvr.dot(Wvr));
  __m128d pk = _mm_load_pd(W_K + 8);

  //=========================Вычисляем релятивистикие параметры============================//

  // 1. / sqrt(1 - VV_L); // фактор Лоренца
  __m128d gk = _mm_div_pd(ones2, _mm_sqrt_pd(_mm_sub_pd(ones2, vk)));

  // 1 + kGamma_g * p_L / d_L; // энтальпия
  __m128d hk = _mm_add_pd(ones2, _mm_mul_pd(_mm_set1_pd(kGamma_g), _mm_div_pd(pk, dk)));

  //((kGamma1 * p_L) / (d_L * h_L)); //квадрат скорости звука
  __m128d cs2 = _mm_div_pd(_mm_mul_pd(_mm_set1_pd(kGamma1), pk), _mm_mul_pd(dk, hk));

  //=========================Вычисляем скорости сигналов============================//
  __m128d Vn2 = _mm_mul_pd(vnk, vnk); //квадрат нормальной компоненты
  __m128d Vt2 = _mm_sub_pd(vk, Vn2);  //квадрат нормы - квадрат компоненты

  __m128d s1 = _mm_sub_pd(_mm_sub_pd(ones2, Vn2), _mm_mul_pd(Vt2, cs2)); //(1 - Vn2 - Vt2 * cs2)

  // sqrt(cs2 * (1 - Vn2 - Vt2 * cs2) * (1 - V2_norm));
  __m128d srootk = _mm_sqrt_pd(_mm_mul_pd(_mm_mul_pd(cs2, s1), _mm_sub_pd(ones2, vk)));

  __m128d s2 = _mm_mul_pd(vnk, _mm_sub_pd(ones2, cs2)); // Vn * (1.0 - cs2)
  __m128d div = _mm_div_pd(ones2, _mm_sub_pd(ones2, _mm_mul_pd(vk, cs2)));
  __m128d cmax = _mm_mul_pd(_mm_add_pd(s2, srootk), div); //(Vn * (1.0 - cs2) + sroot_L) / (1.0 - V2_norm * cs2);
  __m128d cmin = _mm_mul_pd(_mm_sub_pd(s2, srootk), div); //(Vn * (1.0 - cs2) - sroot_L) / (1.0 - V2_norm * cs2);

  //_mm256_store_si256((__m256i*) &c[i], z);

  alignas(16) Type cmax_k[2];
  alignas(16) Type cmin_k[2];

  _mm_store_pd(cmax_k, cmax);
  _mm_store_pd(cmin_k, cmin);

  alignas(16) Type lambda_K[2] = {std::min(cmin_k[0], cmin_k[1]), std::max(cmax_k[0], cmax_k[1])};

  /* --------------------------------------------------
          compute HLLC  flux
     -------------------------------------------------- */
  if (lambda_K[0] >= 0.0) {
    __m256d ul = _mm256_load_pd(U_L);
    __m256d fl = _mm256_mul_pd(ul, _mm256_set1_pd(Wvl[0]));

    _mm256_store_pd(F, fl);
    F[4] = U_L[1];
    F[1] += W_K[8];

    // WRITE_LOG("L: %lf %lf %lf %lf %lf\n", F[0], F[1], F[2], F[3], F[4]);

  } else if (lambda_K[1] <= 0.0) {

    __m256d ur = _mm256_load_pd(U_R);
    __m256d fr = _mm256_mul_pd(ur, _mm256_set1_pd(Wvr[0]));
    _mm256_store_pd(F, fr);

    F[4] = U_R[1];
    F[1] += W_K[9];

    // WRITE_LOG("R: %lf %lf %lf %lf %lf\n", F[0], F[1], F[2], F[3], F[4]);
  } else {

    /* ---------------------------------------
                       get u*
       --------------------------------------- */

    __m128d L2 = _mm_load_pd(lambda_K);

    __m128d Uv = _mm_setr_pd(U_L[1], U_R[1]);
    __m128d Up = _mm_setr_pd(U_L[4], U_R[4]);
    __m128d Fv2 = _mm_mul_pd(Uv, vnk);

    __m128d Ak = _mm_sub_pd(_mm_mul_pd(L2, Up), Uv);
    __m128d Bk = _mm_sub_pd(_mm_sub_pd(_mm_mul_pd(L2, Uv), pk), Fv2);

    alignas(16) Type A[2];
    _mm_store_pd(A, Ak);
    alignas(16) Type B[2];
    _mm_store_pd(B, Bk);

    Type a = A[1] * lambda_K[0] - A[0] * lambda_K[1];
    Type b = A[0] + B[0] * lambda_K[1] - A[1] - B[1] * lambda_K[0];
    Type c = B[1] - B[0];

    /*
          if (fabs(a) > 1.e-9){
            usp = 0.5*(- b + sqrt(b*b - 4.0*a*c))/a;
            usm = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a;
          }else{
            usp = usm = -c/b;
          }
    */
    Type scrh = -0.5 * (b + SIGN(b) * sqrt(b * b - 4.0 * a * c));
    Type us = c / scrh;

    Type ps = (A[0] * us - B[0]) / (1.0 - us * lambda_K[0]);

    __m256d ul = _mm256_load_pd(U_L);
    __m256d ur = _mm256_load_pd(U_R);

    __m256d fl = _mm256_mul_pd(ul, _mm256_set1_pd(Wvl[0]));
    __m256d fr = _mm256_mul_pd(ur, _mm256_set1_pd(Wvr[0]));

    __m128d dif2 = _mm_div_pd(ones2, _mm_sub_pd(L2, _mm_set1_pd(us)));

    alignas(16) Type sUb[2];
    _mm_store_pd(sUb, _mm_sub_pd(L2, vnk));
    alignas(16) Type f_usl[4];
    if (us >= 0.0) {
      Type dd = 1.0 / (lambda_K[0] - us);

      f_usl[0] = U_L[0] * (sUb[0]);
      f_usl[1] = (lambda_K[0] * (U_L[4] + ps) - U_L[1]) * us;
      f_usl[2] = U_L[2] * sUb[0];
      f_usl[3] = U_L[3] * sUb[0];

      __m256d f4 = _mm256_load_pd(f_usl);
      f4 = _mm256_mul_pd(f4, _mm256_set1_pd(dd));
      f4 = _mm256_sub_pd(f4, ul);
      f4 = _mm256_mul_pd(f4, _mm256_set1_pd(lambda_K[0]));
      f4 = _mm256_add_pd(f4, fl);

      _mm256_store_pd(F, f4);
      F[4] = (f_usl[1] * dd);
      F[1] += W_K[8];

      // WRITE_LOG("UL: %lf %lf %lf %lf %lf\n", F[0], F[1], F[2], F[3], F[4]);

    } else {
      Type dd = 1.0 / (lambda_K[1] - us);

      f_usl[0] = U_R[0] * (sUb[1]);
      f_usl[1] = (lambda_K[1] * (U_R[4] + ps) - U_R[1]) * us;
      f_usl[2] = U_R[2] * sUb[1];
      f_usl[3] = U_R[3] * sUb[1];
      __m256d f4 = _mm256_load_pd(f_usl);
      f4 = _mm256_mul_pd(f4, _mm256_set1_pd(dd));
      f4 = _mm256_sub_pd(f4, ur);
      f4 = _mm256_mul_pd(f4, _mm256_set1_pd(lambda_K[1]));
      f4 = _mm256_add_pd(f4, fr);

      _mm256_store_pd(F, f4);
      F[4] = (f_usl[1] * dd);
      F[1] += W_K[9];

      // WRITE_LOG("UR: %lf %lf %lf %lf %lf\n", F[0], F[1], F[2], F[3], F[4]);
    }
  }

  _mm256_storeu_pd(&f.f.d, _mm256_mul_pd(_mm256_load_pd(F), _mm256_set1_pd(f.geo.S)));
  f.f.p = F[4] * f.geo.S;

  f.f.v = (T.transpose()) * f.f.v;
  return std::max(lambda_K[0], lambda_K[1]);
}

#undef MAX_ITER
#endif //! SOLVERS