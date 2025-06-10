#include "linear_alg.h"

void GetRotationMatrix(const Vector3 &n, Matrix3 &T)
{
  T = Matrix3::Zero();
  constexpr Type eps = 1e-12;

  if (LIKELY(fabs(n[2] * n[2] - 1) > eps))
  {

    T(0, 0) = n[0];
    T(0, 1) = n[1];
    T(0, 2) = n[2];

    Type sqr = sqrt(1 - n[2] * n[2]);
    Type sqr_inv = 1. / sqr;

#ifdef RRHD_DEBUG
    if (sqr < eps * eps - eps / 10)
      D_L;
#endif

    T(1, 0) = -n[1] * sqr_inv;
    T(1, 1) = n[0] * sqr_inv;

    T(2, 0) = -n[0] * n[2] * sqr_inv;
    T(2, 1) = -n[1] * n[2] * sqr_inv;
    T(2, 2) = sqr;
  }
  else if (n[2] > 0) // n_z == 1
  {
    T(0, 2) = 1;
    T(1, 1) = 1;
    T(2, 0) = -1;
  }
  else // n_z == -1
  {
    T(0, 2) = -1;
    T(1, 1) = -1;
    T(2, 0) = 1;
  }
}

void GetRotationMatrix2(const Vector3 &n, Matrix3 &T)
{
  T = Matrix3::Zero();

  auto Theta = ([](const Vector3 &v)
                { return acos(v[2]); });

  auto Phi = ([](const Vector3 &v)
              {
    constexpr Type PI = 3.141592653589793;
    Type p = atan2(v[1], v[0]);
    return (p < 0) ? p + 2 * PI : p; });

  Type phi = Phi(n);
  Type th = Theta(n);

  T(0, 0) = cos(th) * sin(phi);
  T(0, 1) = sin(th) * sin(phi);
  T(0, 2) = cos(th);

  T(1, 0) = -sin(th);
  T(1, 1) = cos(th);

  T(2, 0) = cos(th) * cos(phi);
  T(2, 1) = sin(th) * cos(phi);
  T(2, 2) = -sin(phi);
}

Vector3 Spherical2Cartesian(const Vector3 &rpt)
{
  Type x = rpt[0] * sin(rpt[1]) * cos(rpt[2]);
  Type y = rpt[0] * sin(rpt[1]) * sin(rpt[2]);
  Type z = rpt[0] * cos(rpt[1]);
  return Vector3(x, y, z);
}

Vector3 Cartesian2Spherical(const Vector3 &xyz)
{
  Vector3 x;
  Type r = xyz.norm();
  if (r < 1e-10)
  {
    return Vector3::Zero();
  }

  Type theta = acos(xyz[2] / r);
  Type phi = atan2(xyz[1], xyz[0]);

  if (phi > M_PI_2)
  {
    phi -= M_PI;
  }

  if (phi < -M_PI_2)
  {
    phi += M_PI;
  }

  //  return Vector3(r,theta,phi);
  return Vector3(theta, phi, 1);
}
