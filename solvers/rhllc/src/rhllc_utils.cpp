#if defined RHLLC && defined SOLVERS
#include "rhllc_utils.h"
#include "rhllc_flux.h"
#include "rhllc_flux_stab.h"

#include "gas_state.h"

#include "global_value.h"
#include "linear_alg.h"

#include <omp.h>

Type rhllc::max_signal_speed = 1;

Type rhllc::GetTimeStep(const hllc_value_t &hllc_set, const std::vector<elem_t> &cells) {

  const Type t = hllc_set.CFL * hllc_set.h_min / max_signal_speed;

  DIE_IF(t < 0);

  return t;
}

void rhllc::HllcPhysToConv(std::vector<elem_t> &cells) {

#pragma omp parallel for
  for (int i = 0; i < cells.size(); i++) {
    GetConvValueStab(cells[i].phys_val, cells[i].conv_val);
  }
}

void rhllc::HllcConvToPhys(std::vector<elem_t> &cells) {

#pragma omp parallel for
  for (int i = 0; i < cells.size(); i++) {
    GetPhysValueStab(cells[i].phys_val, cells[i].conv_val);
  }
}

void rhllc::BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound) {
  const int id_l = f.geo.id_l;
  const int id_r = f.geo.id_r;
  const elem_t &cell = cells[id_l];
  Matrix3 T;

  switch (id_r) // id соседа она же признак ГУ
  {
  case e_bound_free:
    bound.conv_val = cell.conv_val;
    bound.phys_val = cell.phys_val;
    break;
  case e_bound_inner_source:
#if GEOMETRY_TYPE == Sphere
    bound.phys_val.d = 0.1;
    bound.phys_val.v << 0, 0, 0;
    bound.phys_val.p = 1;

    GetConvValue(bound.phys_val, bound.conv_val);
#elif GEOMETRY_TYPE == Cone_JET
    bound.phys_val.d = 0.1;
    bound.phys_val.v << 0.99, 0, 0;
    bound.phys_val.p = 0.01;
    GetConvValue(bound.phys_val, bound.conv_val);
#else
    bound.conv_val = cell.conv_val;
    bound.phys_val = cell.phys_val;
#endif
    break;
  case e_bound_out_source:
#if GEOMETRY_TYPE == Cylinder
    bound.phys_val.d = 0.1;
    bound.phys_val.v << 0.99, 0, 0;
    bound.phys_val.p = 0.01;
    GetConvValue(bound.phys_val, bound.conv_val);
#elif GEOMETRY_TYPE == Cone

    bound.phys_val.d = kM_hydrogen * 1e15 / kDensity; //  0.1;
    bound.phys_val.p = GetPressure(bound.phys_val.d, 1e8);
    bound.phys_val.v = Vector3(1e3 / kVelocity, 0, 0);

    GetConvValue(bound.phys_val, bound.conv_val);
#else
    bound.conv_val = cell.conv_val;
    bound.phys_val = cell.phys_val;
#endif
    break;
  case e_bound_lock:

    bound.conv_val = cell.conv_val;
    bound.phys_val = cell.phys_val;

    GetRotationMatrix(f.geo.n, T);

    bound.conv_val.v = T * bound.conv_val.v;
    bound.phys_val.v = T * bound.phys_val.v;

    bound.conv_val.v[0] = -bound.conv_val.v[0];
    bound.phys_val.v[0] = -bound.phys_val.v[0];

    T = T.transpose();

    bound.conv_val.v = T * bound.conv_val.v;
    bound.phys_val.v = T * bound.phys_val.v;

    break;

  default:
    DIE_IF(id_r < 0); // Err bound in RHLLC_3d

    bound.conv_val = cells[id_r].conv_val;
    bound.phys_val = cells[id_r].phys_val;
    break;
  }
}

#if NUMBER_OF_MEASUREMENTS == 2
int ReBuildPhysicValue(const Vector4 &U, Vector4 &W) {
  Vector2 v(W(1), W(2));

  const Type vv = v.dot(v);
  const Type d = W(0);
  Type Gamma0 = 1. / sqrt(1 - vv);
  const Type h = 1 + gamma_g * W(3) / d;

  Type W0 = d * h * Gamma0 * Gamma0; // U[0] * Gamma0 * h;

  Vector2 m(U[1], U[2]);
  Type mm = m.dot(m);

  Type p = W(3);

  Type D = U[0];
  Type E = U[3];

  int cc = 0;

  Type err = 1;
  do {
    err = W0;

    Type fW = W0 - p - E;

    Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
    Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
    W0 -= (fW / dFdW);

    Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

    p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

    v[0] = m[0] / W0;
    v[1] = m[1] / W0;
    // v[2] = m[2] / W0;

    err -= W0;
    cc++;
  } while (fabs(err / W0) > 1e-14);

  if (p < 0 || D < 0 || std::isnan(p) || std::isnan(D)) {
    printf("Error (p = %lf, d= %lf)", p, D / Gamma0);
    return 1;
  }

  W(0) = D / Gamma0;
  W(1) = v(0);
  W(2) = v(1);
  W(3) = p;

  return 0;
}

int ReBuildPhysicValue(const std::vector<Vector4> &U, std::vector<Vector4> &W) {

  bool flag = false;
#pragma omp parallel default(none) shared(U, W, flag)
  {
    const int size = U.size();

#pragma omp for
    for (int num_cell = 0; num_cell < size; num_cell++) {
      if (!flag && ReBuildPhysicValue(U[num_cell], W[num_cell])) {
#pragma omp critical
        {
          flag = true;
          printf("Error cell= %d\n", num_cell);
          // MPI_RETURN(1);
        }
      }
    }
  }

  return flag;
}

int ReBuildConvValue(const std::vector<Vector4> &W, std::vector<Vector4> &U) {

  U.resize(size_grid);
  Vector4 cell;

  for (size_t i = 0; i < size_grid; i++) {
    const Type v = (W[i](1) * W[i](1) + W[i](2) * W[i](2));
    const Type d = W[i](0);
    const Type Gamma = 1. / sqrt(1 - v);
    const Type h = 1 + gamma_g * W[i](3) / d;
    const Type dhGG = d * h * Gamma * Gamma;

    cell[0] = Gamma * d;
    cell[1] = dhGG * W[i][1];
    cell[2] = dhGG * W[i][2];

    cell[3] = dhGG - W[i](3);

    U[i] = cell;
  }

  return 0;
}
#endif

#endif // RHLLC