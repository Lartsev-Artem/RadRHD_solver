#if defined SOLVERS && defined RHLLC
#include "rhllc_flux.h"

#include "global_value.h"
#include "linear_alg.h"

#define MAX_ITER 20

void rhllc::GetConvValue(const flux_t &W, flux_t &U) {
  const Type d = W.d;
  const Type Gamma = 1. / sqrt(1 - W.v.dot(W.v));
  const Type h = 1 + kGamma_g * W.p / d;
  const Type dhGG = d * h * Gamma * Gamma;

  U.d = Gamma * W.d;
  U.v = dhGG * W.v;
  U.p = dhGG - W.p;
}
void rhllc::GetConvValueStab(flux_t &W, flux_t &U) {

  const Type h = 1 + kGamma_g * W.p / W.d;

  Type g = W.v.dot(W.v);
  if (g >= 1.0) {
    D_L;
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

void rhllc::GetPhysValue(const flux_t &U, flux_t &W) {

  Type Gamma0 = 1. / sqrt(1 - (W.v.dot(W.v)));
  Type p = W.p;

  const Type h = 1 + kGamma_g * p / W.d;

  Type W0 = W.d * h * Gamma0 * Gamma0; // U[0] * Gamma0 * h;

  Vector3 m = U.v;
  Type mm = m.dot(m);

  Vector3 v = W.v;

  Type D = U.d;
  Type E = U.p;

  int cc = 0;

  Type err = 1;
  do {
    err = W0;

    // W.d * h * Gamma0 * Gamma0 - p - E;
    Type fW = W0 - p - E;

    Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
    Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * kGamma_g));
    W0 -= (fW / dFdW);

    Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

    p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * kGamma_g);

    v = m / W0;

    err -= W0;
    cc++;

  } while (fabs(err / W0) > 1e-14);

  W.d = D / Gamma0;
  W.v = v;
  W.p = p;
}

bool rhllc::GetPhysValueSave(const flux_t &U, flux_t &W) {
  bool call_back_flag = false;
  Type Gamma0 = 1. / sqrt(1 - (W.v.dot(W.v)));
  Type p = W.p;

  const Type h = 1 + kGamma_g * p / W.d;

  Type W0 = W.d * h * Gamma0 * Gamma0; // U[0] * Gamma0 * h;

  Vector3 m = U.v;
  Type mm = m.dot(m);

  Vector3 v = W.v;

  Type D = U.d;
  Type E = U.p;

  int cc = 0;

  Type err = 1;
  do {
    err = W0;

    // W.d * h * Gamma0 * Gamma0 - p - E;
    Type fW = W0 - p - E;

    Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
    Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * kGamma_g));
    W0 -= (fW / dFdW);

    Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

    p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * kGamma_g);

    v = m / W0;

    err -= W0;
    cc++;

    if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d)) {
      //	EXIT(1);
      //	p = max(sqrt(mm) - E, 1e-20);
      //	if (std::isnan(p)) p = 1e-20;
      //	v = m / (E + p);
      //	if (v.norm() > 0.9999999995)
      //	{
      //		v /= 1.0005;
      //	}
      //	Gamma0 = 1. / sqrt(1 - v.dot(v));
      ////	printf("try resolve %.16lf\n", v.norm());
      //	//printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
      //	//D_LD;
      call_back_flag = true;
      break;
    }

  } while (fabs(err / W0) > 1e-14);

  /*if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d) || v.norm()>1)
  {
          printf("W= %lf,(%lf), %lf\n", W.d, W.v.norm(), W.p);
          printf("Error cell (p = %lf, d= %lf, %lf)\n", p, D , Gamma0);
          EXIT(1);
          call_back_flag = 2;
  }*/

  W.d = D / Gamma0;
  if (std::isnan(Gamma0) || std::isinf(Gamma0) || (D / Gamma0) < 1e-10) {
    W.d = 1e-10;
    call_back_flag = true;
  }
  W.v = v;

  W.p = p;
  if (std::isnan(p) || std::isinf(p) || p < 1e-20) {
    W.p = 1e-20;
    call_back_flag = true;
  }

  return call_back_flag;
}

int rhllc::PhysPressureFix(flux_t &U, flux_t &W)
/*!
 *
 * \return Error codes are:
 * - 0 = success
 * - 1 = v^2 > 1
 * - 2 = too many iterations
 *
 *********************************************************************** */
{
  int k, done = 0;
  double D, m, m2, p, Dh, lor, plor;
  double u0, u1, f0, f1, du, umax;

  /* ----------------------------------------------
     1. Solve f(u) = 0 with secant method
     ---------------------------------------------- */

  p = kMinPressure;
  D = U.d;
  m2 = U.v.dot(U.v);
  m = sqrt(m2);
  umax = m / D;
  u0 = umax;

  lor = sqrt(1.0 + u0 * u0);
  plor = p * lor;
  Dh = D + plor * kGamma_g;

  f0 = m / Dh - u0;

  u1 = (-D + sqrt(D * D + 4.0 * m * kGamma_g * p)) / (2.0 * kGamma_g * p);
  done = 0;
  for (k = 1; k < MAX_ITER; k++) {
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

  if (k >= MAX_ITER)
    return 2;

  /* ----------------------------------------------
     2. Solution u has been found.
        Update to converged value of u.
        Also, redefine conserved energy and entropy.
     ---------------------------------------------- */

  lor = sqrt(1.0 + u1 * u1);
  plor = p * lor;
  Dh = D + plor * kGamma_g;

  W.d = U.d / lor;
  W.p = p;

  U.p = Dh * lor - p; /* Redefine energy */

  f0 = 1.0 / (Dh * lor); /* = 1 / W */
  W.v = U.v * f0;

  return 0; /* -- success -- */
}

int rhllc::GetPhysValueStab(const flux_t &U, flux_t &W)
/*!
 *
 * \param [in,out]  u      array of conservative variables (entropy will
 *                         be redefined)
 * \param [out]     v      array of primitive variables
 *
 * \return Error codes:
 *  - 0 = success
 *  - 1 = solution does not exist
 *  - 2 = negative pressure
 *  - 3 = inaccurate solution (debug only)
 *  - 4 = NaN
 *********************************************************************** */
{
  int iter;
  double p, D_1, alpha, alpha2, lor2, lor, m, tol = 1.e-11;
  double tau, theta, h, dh_dp, dh_dtau;
  double yp, dyp, dp, scrh;
  double D, E, m2, Q;

  D = U.d;
  E = U.p;
  m2 = U.v.dot(U.v);

  Q = E - sqrt(m2 + D * D);

  if (Q < 0.0)
    return 1; /* Equation does not admit a solution */

  m = sqrt(m2);
  p = m - E;
  p = std::max(p, 0.0);
  D_1 = 1.0 / D;

  double eps2 = 1.e-12; /* Maximum 1/gamma^2 */
  double pmin = sqrt(m2 / (1.0 - eps2)) - E;

  /* ----------------------------------------------
     1. Solve f(p) = 0 by Newton's method
     ---------------------------------------------- */

  p = std::max(p, pmin);
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
    D_L;
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
  theta = p * tau;

  h = 1.0 + kGamma_g * theta;

  W.d = 1.0 / tau;
  W.p = p;
  scrh = 1.0 / (U.p + p); /* = 1 / W */

  W.v = U.v * scrh;

  return 0; /* -- success -- */
}

void rhllc::GetFlux(const flux_t &conv_val_l, const flux_t &conv_val_r,
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
  const Vector3 Vel_L(W_L.v); // T * velocity[num_cell];
  const Vector3 Vel_R(W_R.v); // T * velocity[neig / 4];

  const Type d_L = W_L.d;
  const Type d_R = W_R.d;

  const Type p_L = W_L.p;
  const Type p_R = W_R.p;

  const Type VV_L = Vel_L.dot(Vel_L);
  const Type VV_R = Vel_R.dot(Vel_R);

  //========================================================================================//

  //=========================Вычисляем релятивистикие параметры============================//
  const Type g_L = 1. / sqrt(1 - VV_L); // фактор Лоренца
  const Type g_R = 1. / sqrt(1 - VV_R);

  const Type h_L = 1 + kGamma_g * p_L / d_L; // энтальпия
  const Type h_R = 1 + kGamma_g * p_R / d_R;

  const Type cs_L2 = ((kGamma1 * p_L) / (d_L * h_L)); //квадрат скорости звука
  const Type cs_R2 = ((kGamma1 * p_R) / (d_R * h_R));

  const Type sigmaS_L = (cs_L2) / (g_L * g_L * (1 - cs_L2)); // что-то для расчета собственных чисел HLL
  const Type sigmaS_R = (cs_R2) / (g_R * g_R * (1 - cs_R2));

  //========================================================================================//

  const Type sqr_L = sqrt(sigmaS_L * (1 - Vel_L[0] * Vel_L[0] + sigmaS_L));
  const Type sqr_R = sqrt(sigmaS_R * (1 - Vel_R[0] * Vel_R[0] + sigmaS_R));

  // здесь встречалась альтернатива сравнения с нулем min(0,L), max(0,R)
  const Type lambda_L = std::min((Vel_L[0] - sqr_L) / (1 + sigmaS_L), (Vel_R[0] - sqr_R) / (1 + sigmaS_R));
  const Type lambda_R = std::max((Vel_L[0] + sqr_L) / (1 + sigmaS_L), (Vel_R[0] + sqr_R) / (1 + sigmaS_R));

  flux_t F;
  if (lambda_R <= 0) // если верно выполнить всегда
  {
    F.d = U_R.d * Vel_R[0];             // D*v_x
    F.v(0) = U_R.v[0] * Vel_R[0] + p_R; // mx*vx+p
    F.v(1) = U_R.v[1] * Vel_R[0];
    F.v(2) = U_R.v[2] * Vel_R[0];
    F.p = U_R.v[0];
    // continue;
  } else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
  {
    F.d = U_L.d * Vel_L[0];             // D*v_x
    F.v(0) = U_L.v[0] * Vel_L[0] + p_L; // mx*vx+p
    F.v(1) = U_L.v[1] * Vel_L[0];
    F.v(2) = U_L.v[2] * Vel_L[0];
    F.p = U_L.v[0];
    // continue;
  } else {
    //====================Расчёт потоков и приближений hll=========================================//
    flux_t F_L;
    flux_t F_R;
    flux_t U_hll;
    flux_t F_hll;

    F_R.d = U_R.d * Vel_R[0];             // D*v_x
    F_R.v(0) = U_R.v[0] * Vel_R[0] + p_R; // mx*vx+p
    F_R.v(1) = U_R.v[1] * Vel_R[0];
    F_R.v(2) = U_R.v[2] * Vel_R[0];
    F_R.p = U_R.v[0];

    F_L.d = U_L.d * Vel_L[0];             // D*v_x
    F_L.v(0) = U_L.v[0] * Vel_L[0] + p_L; // mx*vx+p
    F_L.v(1) = U_L.v[1] * Vel_L[0];
    F_L.v(2) = U_L.v[2] * Vel_L[0];
    F_L.p = U_L.v[0];

    for (int i = 0; i < CELL_SIZE + 1; i++) {
      F_hll[i] = (F_L[i] * lambda_R - F_R[i] * lambda_L + ((U_R[i] - U_L[i]) * lambda_R * lambda_L)) / (lambda_R - lambda_L);
      U_hll[i] = ((U_R[i] * lambda_R) - (U_L[i] * lambda_L) + (F_L[i] - F_R[i])) / (lambda_R - lambda_L);
    }
    /*F_hll = (F_L * lambda_R - F_R * lambda_L + ((U_R - U_L) * lambda_R * lambda_L)) / (lambda_R - lambda_L);
    U_hll = ((U_R * lambda_R) - (U_L * lambda_L) + (F_L - F_R)) / (lambda_R - lambda_L);*/

#ifdef ONLY_RHLL
    F = F_hll;
#else
    //=========================Поиск скорости промежуточной волны===============================//
    const Type a = F_hll.p;               // F_E^hll
    const Type b = -U_hll.p - F_hll.v[0]; // (E_hll + F_mx^hll)
    const Type c = U_hll.v[0];            // mx_hll

#if 1 // как описано в Mignone...
    Type quad = -0.5 * (b + SIGN(b) * sqrt(b * b - 4 * a * c));
    Type _lambda = c / quad;

#endif
    {
      if (_lambda >= 0.0) {
        //============================Поиск промежуточного давления ===================================//
        const Type _p = -F_hll.p * _lambda + F_hll.v[0];
        //============================================================================================//

        //==========================Финальный поток HLLC=============================================//
        flux_t _U_L;
        const Type dif_L = 1.0 / (lambda_L - _lambda);

        _U_L.d = (U_L.d * (lambda_L - Vel_L[0])) * dif_L;
        _U_L.v[0] = (U_L.v[0] * (lambda_L - Vel_L[0]) + _p - p_L) * dif_L;
        _U_L.v[1] = (U_L.v[1] * (lambda_L - Vel_L[0])) * dif_L;
        _U_L.v[2] = (U_L.v[2] * (lambda_L - Vel_L[0])) * dif_L;
        _U_L.p = (U_L.p * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) * dif_L;

        _U_L -= U_L;
        _U_L *= lambda_L;
        F = F_L + _U_L;

        //============================================================================================//
      } else //(_S <= 0)
      {
        //============================Поиск промежуточного давления ===================================//
        const Type _p = -F_hll.p * _lambda + F_hll.v[0];
        //============================================================================================//
        flux_t _U_R;
        const Type dif_R = 1.0 / (lambda_R - _lambda);

        _U_R.d = (U_R.d * (lambda_R - Vel_R[0])) * dif_R;
        _U_R.v[0] = (U_R.v[0] * (lambda_R - Vel_R[0]) + _p - p_R) * dif_R;
        _U_R.v[1] = (U_R.v[1] * (lambda_R - Vel_R[0])) * dif_R;
        _U_R.v[2] = (U_R.v[2] * (lambda_R - Vel_R[0])) * dif_R;
        _U_R.p = (U_R.p * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) * dif_R;

        _U_R -= U_R;
        _U_R *= lambda_R;
        F = F_R + _U_R;
      }
    }
#endif
  }
  // WRITE_LOG(" F0= " << F[0] << ' ' << F[1] << ' ' << F[2] << ' ' << F[3] << ' ' << F[4] << '\n');
  f.f = F;
  f.f.v = (T.transpose()) * F.v;
  // WRITE_LOG(" F1= " << f.f[0] << ' ' << f.f[1] << ' ' << f.f[2] << ' ' << f.f[3] << ' ' << f.f[4] << '\n');
  f.f *= f.geo.S;
}

/*FLUX
 fx[RHO] = u[RHO]*vn;
    fx[MX1] = u[MX1]*vn;
    fx[MX2] = u[MX2]*vn;
    fx[MX3] = u[MX3]*vn;
    fx[ENG] = u[MXn];
    state->prs[i] = v[PRS];
 */
void rhllc::GetFluxStab(const flux_t &conv_val_l, const flux_t &conv_val_r,
                        const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f)
// const Sweep *sweep, int beg, int end, double *cmax, Grid *grid
/*!
 * Solve the RHD Riemann problem using the HLLC Riemann solver.
 *
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
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
  const Vector3 Vel_L(W_L.v); // T * velocity[num_cell];
  const Vector3 Vel_R(W_R.v); // T * velocity[neig / 4];

  const Type d_L = W_L.d;
  const Type d_R = W_R.d;

  const Type p_L = W_L.p;
  const Type p_R = W_R.p;

  const Type VV_L = Vel_L.dot(Vel_L);
  const Type VV_R = Vel_R.dot(Vel_R);

  //========================================================================================//

  //=========================Вычисляем релятивистикие параметры============================//
  const Type g_L = 1. / sqrt(1 - VV_L); // фактор Лоренца
  const Type g_R = 1. / sqrt(1 - VV_R);

  const Type h_L = 1 + kGamma_g * p_L / d_L; // энтальпия
  const Type h_R = 1 + kGamma_g * p_R / d_R;

  const Type cs_L2 = ((kGamma1 * p_L) / (d_L * h_L)); //квадрат скорости звука
  const Type cs_R2 = ((kGamma1 * p_R) / (d_R * h_R));

  /* ----------------------------------------------------
       compute sound speed & fluxes at zone interfaces
     ---------------------------------------------------- */
  const Type Vn2_L = Vel_L[0] * Vel_L[0];
  const Type Vt2_L = VV_L - Vn2_L;
  const Type sroot_L = sqrt(cs_L2 * (1 - Vn2_L - Vt2_L * cs_L2) * (1 - VV_L));
  const Type cmax_L = (Vel_L[0] * (1.0 - cs_L2) + sroot_L) / (1.0 - VV_L * cs_L2);
  const Type cmin_L = (Vel_L[0] * (1.0 - cs_L2) - sroot_L) / (1.0 - VV_L * cs_L2);

  const Type Vn2_R = Vel_R[0] * Vel_R[0];
  const Type Vt2_R = VV_R - Vn2_R;
  const Type sroot_R = sqrt(cs_R2 * (1 - Vn2_R - Vt2_R * cs_R2) * (1 - VV_R));
  const Type cmax_R = (Vel_R[0] * (1.0 - cs_R2) + sroot_R) / (1.0 - VV_R * cs_R2);
  const Type cmin_R = (Vel_R[0] * (1.0 - cs_R2) - sroot_R) / (1.0 - VV_R * cs_R2);

  const Type lambda_L = std::min(cmin_L, cmin_R);
  const Type lambda_R = std::max(cmax_L, cmax_R);

  // fx[RHO] = u[RHO]*vn;
  // fx[MX1] = u[MX1]*vn;
  // fx[MX2] = u[MX2]*vn;
  // fx[MX3] = u[MX3]*vn;
  // fx[ENG] = u[MXn];
  // state->prs[i] = v[PRS];
  // Flux (stateL, beg, end);

  /* --------------------------------------------------
          compute HLLC  flux
     -------------------------------------------------- */
  // flux_t F_L(U_L.d * Vel_L[0], U_L.v * Vel_L[0], U_L.v[0]);
  // flux_t F_R(U_R.d * Vel_R[0], U_R.v * Vel_R[0], U_R.v[0]);
  // flux_t F_L(U_L.d * Vel_L[0], conv_val_l.v * Vel_L[0], U_L.v[0]);
  // flux_t F_R(U_R.d * Vel_R[0], conv_val_r.v * Vel_R[0], U_R.v[0]);

  flux_t F;
  if (lambda_L >= 0.0) {
    F.d = U_L.d * Vel_L[0]; // D*v_x
    F.v = U_L.v * Vel_L[0]; // mx*vx+p
    F.p = U_L.v[0];

    F.v[0] += p_L;
    //  sweep->press[i] = pL[i];

  } else if (lambda_R <= 0.0) {

    F.d = U_R.d * Vel_R[0]; // D*v_x
    F.v = U_R.v * Vel_R[0]; // mx*vx+p
    F.p = U_R.v[0];
    F.v[0] += p_R;

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

    // F_L *= lambda_R;
    // F_R *= lambda_L;
    // F_L -= F_R;
    // U_R -= U_L; //повернутые или нет
    // U_R *= (lambda_L * lambda_R);
    // F = U_R + F_L;
    // F *= scrh;

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

    // flux_t usl, usr;
    if (us >= 0.0) {
      flux_t Fusl;
      Fusl.d = U_L.d * (lambda_L - Vel_L[0]) / (lambda_L - us);
      Fusl.v[0] = (lambda_L * (U_L.p + ps) - U_L.v[0]) * us / (lambda_L - us);
      Fusl.v[1] = U_L.v[1] * (lambda_L - Vel_L[0]) / (lambda_L - us);
      Fusl.v[2] = U_L.v[2] * (lambda_L - Vel_L[0]) / (lambda_L - us);
      Fusl.p = U_L.p + (Fusl.v[0] - U_L.v[0]) / lambda_L;

      Fusl -= U_L;
      Fusl *= lambda_L;
      F = F_L + Fusl;
      // F = F_L + lambda_L *(Fusl - U_L);
      //  sweep->press[i] = pL[i];
      F.v[0] += p_L;
    } else {
      flux_t Fusr;
      Fusr.d = U_R.d * (lambda_R - Vel_R[0]) / (lambda_R - us);
      Fusr.v[0] = (lambda_R * (U_R.p + ps) - U_R.v[0]) * us / (lambda_R - us);
      Fusr.v[1] = U_R.v[1] * (lambda_R - Vel_R[0]) / (lambda_R - us);
      Fusr.v[2] = U_R.v[2] * (lambda_R - Vel_R[0]) / (lambda_R - us);
      Fusr.p = U_R.p + (Fusr.v[0] - U_R.v[0]) / lambda_R;

      Fusr -= U_R;
      Fusr *= lambda_R;
      F = F_R + Fusr;
      F.v[0] += p_R;
      // sweep->press[i] = pR[i];
    }
#endif // HLL

  } /* -- end loop on points -- */
  //  WRITE_LOG("F: %lf, %lf, %lf, %lf, %lf\n", F.d, F.v[0], F.v[1], F.v[2], F.p);

  f.f = F;
  f.f.v = (T.transpose()) * F.v;
  f.f *= f.geo.S;
}

#undef MAX_ITER
#endif //! SOLVERS