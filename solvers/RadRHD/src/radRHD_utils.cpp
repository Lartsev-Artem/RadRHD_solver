#if defined RAD_RHD
#include "radRHD_utils.h"
#include "radRHD_config.h"

#include "gas_state.h"
#include "plunk.h"

Type rad_rhd::EddTensor(int i, int j, const Type U, const Vector3 &F) {
  /* --------------------------------------------------------
                  Compute the (i,j) component of the Eddington tensor
                  using the M1 closure, taking as input the state v.
                  (i,j) values must always be within FR1 and FR3.
     -------------------------------------------------------- */

  double ni, nj, f2, chi_2, edt;

  double flux_module = F.dot(F);

  if (flux_module < 1e-150 || U < 1e-150) {
    if (i == j)
      return 0.33333333333333;
    return 0.0;
  }

  flux_module = sqrt(flux_module);

  ni = F[i] / flux_module;
  nj = F[j] / flux_module;

  if (flux_module > U)
    return ni * nj;

  f2 = flux_module / U;
  f2 = f2 * f2;
  chi_2 = (1.5 + 2.0 * f2) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f2));
  edt = (3.0 * chi_2 - 0.5) * ni * nj;

  if (i == j)
    edt += 0.5 - chi_2;

  return edt;
}

void rad_rhd::GetRadSource(const int cell, const grid_t &grid, Vector4 &G) {

  /*-- Set opacities --*/
  const Type abs_op = grid.cells[cell].illum_val.absorp_coef;
  const Type scat_op = grid.cells[cell].illum_val.scat_coef;
  const Type tot_op = abs_op + scat_op;

  const Type rho = grid.cells[cell].phys_val.d;
  const Type prs = grid.cells[cell].phys_val.p;
  const Vector3 &v = grid.cells[cell].phys_val.v;

#ifdef USE_CUDA
  const Type Ur = grid.energy[cell];
  const Vector3 &Fr = grid.stream[cell];
  const Matrix3 &Tr = grid.impuls[cell];
#else
  const Type Ur = grid.cells[cell].illum_val.energy;
  const Vector3 &Fr = grid.cells[cell].illum_val.stream;
  const Matrix3 &Tr = grid.cells[cell].illum_val.impuls;
#endif

  /*-- Compute gamma, gamma^2, u and u^2 --*/
  Type gamma = v.dot(v);
  Type u2 = gamma;
  Type gamma2 = 1.0 / (1.0 - gamma);
  u2 *= gamma2;
  gamma = sqrt(gamma2);

  Vector3 u = gamma * v;

  /*-- Compute products involving the proper velocity --*/
  Type uuD = 0.;
  Type uF = u.dot(Fr);
  Vector3 uD;
  Type D;
  for (int i = 0; i < 3; i++) {
    uD[i] = 0.;
    for (int j = 0; j < 3; j++) {
#if ISOTROPIC_INTENSITY == 1
      D = EddTensor(j, i, Ur, Fr);
#else
      D = Tr(j, i);
#endif
      uD[i] += u[j] * D;
      uuD += u[i] * u[j] * D;
    }
  }

  /*-- Compute some useful quantities --*/
  Type rhogamma = rho * gamma;
  Type q = -abs_op * rho * B_Plank(GetTemperature(rho, prs));

  /*-- Compute source function in the Eulerian frame --*/
  G[0] = q * gamma + rhogamma * (abs_op - scat_op * (u2 + uuD)) * Ur - rho * (abs_op - scat_op * (u2 + gamma2)) * uF;

  for (int i = 0; i < 3; i++)
    G[i + 1] = q * u[i] - rho * (tot_op * uD[i] + scat_op * u[i] * (gamma2 + uuD)) * Ur + rhogamma * (tot_op * Fr[i] + scat_op * 2.0 * uF * u[i]);
}

#endif //! RAD_RHD

#if 0
void MyRadSourceFunction(int cell, grid_t &grid, Vector4 &G) {
  int i, j;
  static int comps = 4;
  double u2, gamma, gamma2, D, uD[4], uuD, uF;
  double B, rhogamma, q;

  /*-- Set opacities --*/
  double abs_op, scat_op, tot_op;
  abs_op = g_absorptionCoeff;
  scat_op = g_scatteringCoeff;
  tot_op = g_totalOpacity;

  const Vector3 &v = grid.cells[cell].phys_val.v;
  const Type rho = grid.cells[cell].phys_val.d;
  const Type prs = grid.cells[cell].phys_val.p;
  const Type Ur = grid.energy[cell];
  const Vector3 &Fr = grid.stream[cell];
  const Matrix3 &Tr = grid.impuls[cell];

  Matrix4 T;
  T(0, 0) = Ur;
  T(1, 0) = T(0, 1) = Fr[0];
  T(2, 0) = T(0, 2) = Fr[1];
  T(3, 0) = T(0, 3) = Fr[2];
  for (size_t j = 0; j < 3; j++) {
    for (size_t k = 0; k < 3; k++)
      T(j + 1, k + 1) = Tr(j, k);
  }
  gamma = v.dot(v);
  gamma2 = 1.0 / (1.0 - gamma);
  gamma = sqrt(gamma2);
  Vector4 u(gamma, v[0] * gamma, v[1] * gamma, v[2] * gamma);

  /*-- Compute products involving the proper velocity --*/
  uuD = 0.;
  for (i = 0; i < comps; i++) {
    uD[i] = 0.;
    for (j = 0; j < comps; j++) {
      uD[i] += u[j] * T(i, j);
      uuD += u[i] * u[j] * T(i, j);
    }
  }

  /*-- Compute some useful quantities --*/
  rhogamma = rho * gamma;
  B = Blackbody(GetTemperature(rho, prs));

  for (i = 0; i < comps; i++) {
    G[i] = -abs_op * rho * (B * u[i] + uD[i]) - scat_op * rho * (uD[i] + uuD * u[i]);
  }
}
#endif