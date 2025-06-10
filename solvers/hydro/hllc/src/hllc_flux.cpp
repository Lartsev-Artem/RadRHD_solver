
#if defined SOLVERS && defined HLLC
#include "hllc_flux.h"

#include "global_value.h"
#include "hllc_utils.h"

#include "linear_alg.h"

void hllc::GetConvValue(const flux_t &W, flux_t &U) {
  const Type v = W.v.dot(W.v);
  const Type d = W.d;

  U.d = d;
  U.v = d * W.v;
  U.p = W.p / (kGamma1 - 1) + d * v / 2;
}

void hllc::GetPhysValue(const flux_t &U, flux_t &W) {
  const Type d = U.d;
  W.d = d;
  W.v = U.v / d;
  const Type vv = W.v.dot(W.v);
  W.p = (U.p - vv * d / 2.) * (kGamma1 - 1);
}

void hllc::GetFlux(const flux_t &conv_val_l, const flux_t &conv_val_r,
                   const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f) {
  Matrix3 T;
  GetRotationMatrix(f.geo.n, T);

  flux_t U_L = conv_val_l;
  U_L.v = T * conv_val_l.v;

  flux_t U_R = conv_val_r;
  U_R.v = T * conv_val_r.v;

  Type d_L = U_L.d; // density[num_cell];
  Type d_R = U_R.d; // density[neig];

  Vector3 vel = U_L.v / d_L;
  Type v = vel.dot(vel);
  Type p_L = (U_L.p - v * d_L / 2.) * (kGamma1 - 1);

  vel = U_R.v / d_R;
  v = vel.dot(vel);
  Type p_R = (U_R.p - v * d_R / 2.) * (kGamma1 - 1);

  const Type v_L = U_L[1] / d_L; // sqrt(U_L[1] * U_L[1] + U_L[2] * U_L[2] + U_L[3] * U_L[3]);  //velocity[num_cell].norm();
  const Type v_R = U_R[1] / d_R; // sqrt(U_R[1] * U_R[1] + U_R[2] * U_R[2] + U_R[3] * U_R[3]); //velocity[neig].norm();

  const Type a_L = sqrt(kGamma1 * p_L / d_L);
  const Type a_R = sqrt(kGamma1 * p_R / d_R);

  const Type P = (p_L + p_R) / 2;
  const Type Den = (d_L + d_R) / 2;
  const Type A = (a_L + a_R) / 2;
  const Type _p = std::max(0.0, P - (v_R - v_L) * Den * A / 2);

  // pressure-based wave speed estimates
  Type q_L = 1;
  const Type G = (kGamma1 + 1) / (2 * kGamma1);
  if (_p > p_L)
    q_L = sqrt(1 + G * (_p / p_L - 1));

  Type q_R = 1;
  if (_p > p_R)
    q_R = sqrt(1 + G * (_p / p_R - 1));

  const Type S_L = v_L - a_L * q_L;
  const Type S_R = v_R + a_R * q_R;

  flux_t F;
  if (S_R <= 0) // если верно выполнить всегда
  {
    F.d = U_R[1];
    F.v(0) = U_R[1] * U_R[1] / d_R + p_R;
    F.v(1) = U_R[1] * U_R[2] / d_R;
    F.v(2) = U_R[1] * U_R[3] / d_R;
    F.p = (U_R[4] + p_R) * U_R[1] / d_R;
    // continue;
  } else if (S_L >= 0) // выполнить либо по условию либо для всех границ
  {
    F.d = U_L[1];
    F.v(0) = U_L[1] * U_L[1] / d_L + p_L;
    F.v(1) = U_L[1] * U_L[2] / d_L;
    F.v(2) = U_L[1] * U_L[3] / d_L;
    F.p = (U_L[4] + p_L) * U_L[1] / d_L; // TU[4]*d_L
                                         // continue;
  } else {
    const Type roS_L = d_L * (S_L - v_L);
    const Type roS_R = d_R * (S_R - v_R);
    const Type _S = (p_R - p_L + roS_L * v_L - roS_R * v_R) / (roS_L - roS_R);

    const Type P_LR = (p_L + p_R + roS_L * (_S - v_L) + roS_R * (_S - v_R)) / 2.0;

    flux_t D(0, Vector3(1, 0, 0), _S);

    if (_S >= 0) {
      F.d = U_L[1];
      F.v(0) = U_L[1] * U_L[1] / d_L + p_L;
      F.v(1) = U_L[1] * U_L[2] / d_L;
      F.v(2) = U_L[1] * U_L[3] / d_L;
      F.p = (U_L[4] + p_L) * U_L[1] / d_L;

      U_L *= S_L;
      U_L -= F;
      U_L *= _S;

      D *= P_LR;
      D *= S_L;

      F = (U_L + D);
      F /= (S_L - _S);
    } else //(_S <= 0)
    {
      F.d = U_R[1];
      F.v(0) = U_R[1] * U_R[1] / d_R + p_R;
      F.v(1) = U_R[1] * U_R[2] / d_R;
      F.v(2) = U_R[1] * U_R[3] / d_R;
      F.p = (U_R[4] + p_R) * U_R[1] / d_R;

      flux_t buf = F;
      U_R *= S_R;
      U_R -= F;
      U_R *= _S;

      D *= S_R;
      D *= P_LR;

      F = (U_R + D);
      F /= (S_R - _S);
    }
  }

  f.f = F;
  f.f.v = (T.transpose()) * F.v;
  f.f *= f.geo.S;
}

#endif //! SOLVERS &&  HLLC