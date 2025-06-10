#if defined RHLLC
#include "rhllc_utils.h"

#include "gas_state.h"
#include "linear_alg.h"
#include "rhllc_flux.h"

#include <omp.h>

Type rrhd::rhllc::max_signal_speed = 1;

using namespace rrhd;

void rhllc::HllcPhysToConv(grid_t &grid) {

#pragma omp parallel for
  for (int i = 0; i < grid.size; i++) {
    GetConvValue(grid.cells[i].phys_val, grid.cells[i].conv_val);
  }
}

void rhllc::HllcConvToPhys(grid_t &grid) {

  int myid = get_mpi_id();
  int start = grid.loc_shift;
  int end = start + grid.loc_size;

#pragma omp parallel default(none) firstprivate(start, end) shared(grid, glb_files)
  {
#pragma omp for
    for (int i = start; i < end; i++) {
      if (UNLIKELY(GetPhysValue(grid.cells[i].conv_val, grid.cells[i].phys_val))) {
        DIE_IF(PhysPressureFix(grid.cells[i].conv_val, grid.cells[i].phys_val));
      }
    }
  }
}

#ifdef RAD_RHD
#include "radRHD_utils.h"
void rhllc::AddRadFlux(grid_t &grid) {

  int myid = get_mpi_id();
  int start = grid.loc_shift;
  int end = start + grid.loc_size;

  constexpr Type ds = 1.0;
  Type tau = ds * _hllc_cfg.tau;
#pragma omp parallel default(none) firstprivate(start, end, tau) shared(grid)
  {

#pragma omp for
    for (int cell = start; cell < end; cell++) {

      Vector4 G;
      rad_rhd::GetRadSourceOpt(cell - start, grid.cells[cell], grid, G);
      G *= tau;

      flux_t &conv_val = grid.cells[cell].conv_val;
      conv_val.p += G[0];
      conv_val.v[0] += G[1];
      conv_val.v[1] += G[2];
      conv_val.v[2] += G[3];
    }
  }
}
#endif

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