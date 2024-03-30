#include "linear_alg.h"

void GetRotationMatrix(const Vector3 &n, Matrix3 &T) {
  T = Matrix3::Zero();
  constexpr Type eps = 1e-12;

  if (LIKELY(fabs(n[2] * n[2] - 1) > eps)) {

    T(0, 0) = n[0];
    T(0, 1) = n[1];
    T(0, 2) = n[2];

    Type sqr = sqrt(1 - n[2] * n[2]);
    Type sqr_inv = 1. / sqr;

#ifdef DEBUG
    if (sqr < eps * eps - eps / 10)
      D_L;
#endif

    T(1, 0) = -n[1] * sqr_inv;
    T(1, 1) = n[0] * sqr_inv;

    T(2, 0) = -n[0] * n[2] * sqr_inv;
    T(2, 1) = -n[1] * n[2] * sqr_inv;
    T(2, 2) = sqr;
  } else if (n[2] > 0) // n_z == 1
  {
    T(0, 2) = 1;
    T(1, 1) = 1;
    T(2, 0) = -1;
  } else // n_z == -1
  {
    T(0, 2) = -1;
    T(1, 1) = -1;
    T(2, 0) = 1;
  }
}