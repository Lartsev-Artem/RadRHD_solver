#if defined SOLVERS && defined RHLLC
#include "rhllc_flux.h"

#include "global_value.h"
#include "linear_alg.h"

void rhllc::GetConvValue(const flux_t &W, flux_t &U) {
  const Type d = W.d;
  const Type Gamma = 1. / sqrt(1 - W.v.dot(W.v));
  const Type h = 1 + kGamma_g * W.p / d;
  const Type dhGG = d * h * Gamma * Gamma;

  U.d = Gamma * W.d;
  U.v = dhGG * W.v;
  U.p = dhGG - W.p;
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

  const Type cs_L = sqrt((kGamma1 * p_L) / (d_L * h_L)); // скорость звука
  const Type cs_R = sqrt((kGamma1 * p_R) / (d_R * h_R));

  const Type sigmaS_L = (cs_L * cs_L) / (g_L * g_L * (1 - cs_L * cs_L)); // что-то для расчета собственных чисел HHL
  const Type sigmaS_R = (cs_R * cs_R) / (g_R * g_R * (1 - cs_R * cs_R));

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
#endif

    //============================================================================================//
#ifndef ONLY_RHLL
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

#endif //! SOLVERS