#if defined SOLVERS && defined RHLLC && defined ILLUM
#include "global_value.h"
#include "rhllc_main.h"

#include "rhllc_calc.h"
#include "rhllc_flux_stab.h"
#include "rhllc_init.h"
#include "rhllc_utils.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

#include <chrono>
namespace tick = std::chrono;

double g_absorptionCoeff = 0.5; // Absorption coefficient (code_length^2/code_energy)
double g_scatteringCoeff = 0.5; // Scattering coefficient (code_length^2/code_energy)
double g_radiationConst = 1.0;  // Radiation constant (4*StefanBoltzmannConstant/c) (code_energy/(code_length^3*code_temperature^4))
double g_idealGasConst = 1.0;   //(code_temperature) g_idealGasConst*prs/rho determines the gas temperature.
double g_totalOpacity = 1.0;

#define RADIATION_MIN_ERAD 1e-16
double Blackbody(double temperature)
/*!
 * Return the blackbody intensity corresponding to the input temperature.
 *
 * \param [in]	temperature	Input temperature.
 *
 * \return	Return (4*PI)B(T)=radiationConst*T^4.
 *********************************************************************** */
{
  double T4 = temperature * temperature;
  T4 = T4 * T4;
  return g_radiationConst * T4;
}

double GetTemperature(double rho, double prs)
/*!
 * Return the (ideal) gas temperature corresponding to given gas pressure
 * and density.
 *
 * \param [in]	rho		Input mass density.
 * \param [in]	prs		Input (gas) pressure.
 *
 * \return	Return T = g_idealGasConst*prs/rho, where g_idealGasConst is the
 *			equilibrium ratio prs/(T*rho) for the considered ideal gas.
 *********************************************************************** */
{
  return g_idealGasConst * prs / rho;
}

double EddTensor(const Vector3 &F, const double U, int i, int j) {
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

void RadSourceFunction(int cell, grid_t &grid, Vector4 &G) {
  int i, j;
  static int comps = 4 - 1;
  double u[3], u2, gamma, gamma2, D, uD[3], uuD, uF;
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

  /*-- Compute gamma, gamma^2, u and u^2 --*/
  gamma = v.dot(v);
  u2 = gamma;
  gamma2 = 1.0 / (1.0 - gamma);
  u2 *= gamma2;
  gamma = sqrt(gamma2);
  for (i = 0; i < comps; i++)
    u[i] = gamma * v[i];

  /*-- Compute products involving the proper velocity --*/
  uuD = 0.;
  uF = 0.;
  for (i = 0; i < comps; i++) {
    uD[i] = 0.;
    uF += u[i] * Fr[i];
    for (j = 0; j < comps; j++) {
      D = EddTensor(Fr, Ur, j, i);
      uD[i] += u[j] * D;
      uuD += u[i] * u[j] * D;
    }
  }

  /*-- Compute some useful quantities --*/
  rhogamma = rho * gamma;
  B = Blackbody(GetTemperature(rho, prs));
  q = -abs_op * rho * B;

  /*-- Compute source function in the Eulerian frame --*/
  G[0] = q * gamma + rhogamma * (abs_op - scat_op * (u2 + uuD)) * Ur - rho * (abs_op - scat_op * (u2 + gamma2)) * uF;

  for (i = 0; i < comps; i++)
    G[i + 1] = q * u[i] - rho * (tot_op * uD[i] + scat_op * u[i] * (gamma2 + uuD)) * Ur + rhogamma * (tot_op * Fr[i] + scat_op * 2.0 * uF * u[i]);
}

#if 0
int RadIterToPrim(double *x, grid_t &rad_data)
/*!
 * Conver iterated to primitive fields if iterations are carried on MHD fields
 *
 * \param [in]      x       	vector of iterated fields
 * \param [in,out]  rad_data	pointer to rad_data structure that contains
 *                            the primitive fields and the current position
 * \return Returns 0 if success, and 1 if prs == NaN is found.
 *
 *********************************************************************** */
{
  int i, pos;
  static int comps = 4 - 1;
  double *v;
  double *u = rad_data->u;
  double u2, gamma, gamma2, gm1;
  double rhogamma, rhoh, vel2, vB, B2;
  double th, th2, wt;

  double gmmr = kGamma1 / (kGamma1 - 1.0);

  /*-- Get primitive variables at the current position --*/
  pos = rad_data->pos;
  v = rad_data->pv[pos];

  /*-- Set pressure --*/
  v[PRS] = x[0];

  /*-- Check p != NaN --*/
  if (v[PRS] != v[PRS]) {
    D_L;
    // print("! RadIterToPrim: NaN found while setting pressure, ");
    return 1;
  }

  /*-- Check p > 0 --*/
  if (v[PRS] < 0.0) {
    D_L;
    // WARNING(
    //     print("! RadIterToPrim: negative pressure (%8.2e) during implicit step, ", v[PRS]);
    v[PRS] = 1.e-20;
  }

  /*-- Set 4-velocity and Lorentz factor --*/
  u2 = 0.;
  for (i = 0; i < comps; i++) {
    u[i] = x[i + 1];
    u2 += u[i] * u[i];
  }
  gamma2 = u2 + 1.0;
  gamma = sqrt(gamma2);
  gm1 = 1.0 / gamma;

  /*-- Store gamma, gamma2 and u2 --*/
  rad_data->u2 = u2;
  rad_data->gamma = gamma;
  rad_data->gamma2 = gamma2;

  /*-- Set coordinate velocity --*/
  for (i = 0; i < comps; i++)
    v[VX1 + i] = gm1 * u[i];

  /*-- Set density --*/
  rhogamma = rad_data->cv[pos][RHO];
  v[RHO] = gm1 * rhogamma;

  /*-- Set enthalpy density --*/
  rhoh = v[RHO] + gmmr * v[PRS];

  /*-- Set radiation fields --*/
  v[ENR] = rad_data->Ttot[0] - rhoh * gamma2 + v[PRS];
  for (i = 0; i < comps; i++)
    v[FR1 + i] = rad_data->Ttot[i + 1] - rhoh * gamma * u[i];

  /*-- Check Er > 0 --*/
  if (v[ENR] < 0.0) {
    D_L;
    // print("! RadIterToPrim: negative radiation energy (%8.2e) during implicit step, ", v[ENR]);
    v[ENR] = RADIATION_MIN_ERAD;
  }

  return 0;
}

void RadNewtonMinusF(Rad_data *rad_data, double *x, double *mf)
/*!
 * If Newton's method is used, compute and store in mf the components
 * of -F, being F the function whose roots are being looked for. Otherwise,
 * if RADIATION_FIXEDPOINT_GAS, perform a fixed-point iteration.
 *
 * \param [in,out]  rad_data	pointer to rad_data structure
 * \param [in]      x       	vector of iterated fields
 * \param [in,out]  mf	      vector that stores the value of -F(x)
 *
 *********************************************************************** */
{
  int i, j, pos;
  static int comps = RADIATION_NEQS - 1;
  double *u, u2, gamma, gamma2;
  double D, uD[3], uuD, uF;
  double B, rhogamma, q, dt;

  double abs_op, scat_op, tot_op;

  double *v, *consv;

  pos = rad_data->pos;
  v = rad_data->pv[pos];
  consv = rad_data->cv[pos];

  u = rad_data->u;
  dt = rad_data->dt;

  /*-- Compute conserved fields --*/
  consv[ENR] = x[0];
  consv[ENG] = rad_data->Ttot[0] - x[0];
  for (i = 0; i < comps; i++) {
    consv[FR1 + i] = x[i + 1];
    consv[MX1 + i] = rad_data->Ttot[i + 1] - x[i + 1];
  }

  /*-- Convert conserved to primitive fields --*/
  ConsToPrim(rad_data->cv, rad_data->pv, pos, pos, rad_data->flag);

  /*-- Compute 4-velocity and derived quantities --*/
  gamma = consv[RHO] / v[RHO];
  gamma2 = gamma * gamma;
  u2 = gamma2 - 1.0;

  for (i = 0; i < comps; i++)
    u[i] = gamma * v[VX1 + i];

  /*-- Set opacities --*/
  abs_op = g_absorptionCoeff;
  scat_op = g_scatteringCoeff;
  tot_op = g_totalOpacity;

  /*-- Compute products involving the proper velocity --*/
  uuD = 0.;
  uF = 0.;
  for (i = 0; i < comps; i++) {
    uD[i] = 0.;
    uF += u[i] * v[FR1 + i];
    for (j = 0; j < comps; j++) {
      D = EddTensor(v, FR1 + j, FR1 + i);
      uD[i] += u[j] * D;
      uuD += u[i] * u[j] * D;
    }
  }

  /*-- Compute some useful quantities --*/
  rhogamma = consv[RHO];
  B = Blackbody(GetTemperature(v[RHO], v[PRS]));
  q = -abs_op * v[RHO] * B;

  /*-- Compute -F[x] or update iterated fields if
       RADIATION_FIXEDPOINT_GAS is used --*/
  mf[0] = q * gamma + rhogamma * (abs_op - scat_op * (u2 + uuD)) * v[ENR] - v[RHO] * (abs_op - scat_op * (u2 + gamma2)) * uF;

  consv[ENG] = rad_data->Rini[0] + dt * mf[0];
  consv[ENR] = rad_data->Ttot[0] - consv[ENG];

  for (i = 0; i < comps; i++) {
    mf[i + 1] = q * u[i] - v[RHO] * (tot_op * uD[i] + scat_op * u[i] * (gamma2 + uuD)) * v[ENR] + rhogamma * (tot_op * v[FR1 + i] + scat_op * 2.0 * uF * u[i]);

    consv[MX1 + i] = rad_data->Rini[i + 1] + dt * mf[i + 1];
    consv[FR1 + i] = rad_data->Ttot[i + 1] - consv[MX1 + i];
  }
}

void RadNewtonJacobian(double *x, double *mf, double **J, Rad_data *rad_data)
/*!
 * Compute and store in J the Jacobian matrix (dF^i/dx^j) of the system F(x) == 0,
 * calculated using difference quotients with small variations of the variables.
 * Store in mf the current value of -F(x).
 *
 * \param [in]      x       	vector of iterated fields
 * \param [in,out]  mf	      vector that stores the value of -F(x)
 * \param [in,out]  J	        Jacobian matrix
 * \param [in,out]  rad_data	pointer to rad_data structure, used in this case
 *                            to store the initial values of the radiation fields
 *                            as well as primitive fields
 *
 *********************************************************************** */
{
  int i, j;
  // * y, * mf2, * dvars ;
  double vel2, gvel;

  Vector4 dvars;

  /*- Define increments to compute derivatives -*/
  for (i = 0; i < 4; i++)
    dvars[i] = (fabs(x[i]) > 1e-20) ? -1e-4 * x[i] : 1e-20;

  /*- Store the current value of -F(x) -*/
  RadNewtonMinusF(rad_data, x, mf);

  /*- Update auxiliary variables for error computation -*/
  rad_data->exv = rad_data->pv[rad_data->pos][ENR];
  // rad_data->exv = rad_data->pv[rad_data->pos][PRS];

  /*- Compute and store the elements of the Jacobian -*/
  for (i = 0; i < 4; i++)
    y[i] = x[i];

  for (j = 0; j < 4; j++) {
    y[j] += dvars[j];

    RadNewtonMinusF(rad_data, y, mf2);

    for (i = 0; i < 4; i++)
      J[i][j] = (mf[i] - mf2[i]) / dvars[j];
    y[j] = x[j];
  }
}

double RadErr(double *x, double *dx, Rad_data *rad_data)
/*!
 * Compute the sum of the squares of the relative differences of the
 * several fields considered in the implicit step.
 *
 * \param [in]      x         Vector of primitive fields if RADIATION_FIXEDPOINT_RAD
 *                            and conserved fields if RADIATION_FIXEDPOINT_GAS.
 *                            Otherwise, it contains the iterated fields.
 * \param [in]      dx        Vector that stores the variations of
 *                            the iterated fields if Newton's method is used.
 * \param [in,out]  rad_data	Pointer to rad_data structure, used to store
 *                            fields from the previous iteration.
 *
 *********************************************************************** */
{
  int i;
  static int comps = 4 - 1;
  double err, er1, erf, mod;

  /*-- Compute sum of squared relative differences --*/
  er1 = fabs(x[ENR] - rad_data->Rprev[0]) / fabs(rad_data->Rprev[0]);
  err = er1 * er1;

  mod = 0;
  er1 = 0;
  for (i = 0; i < comps; i++) {
    mod += rad_data->Rprev[i + 1] * rad_data->Rprev[i + 1];
    erf = x[FR1 + i] - rad_data->Rprev[i + 1];
    er1 += erf * erf;
  }
  err += (mod > 1e-40) ? er1 / mod : er1 / 1e-40;

  /*-- Update previous step variables --*/
  rad_data->Rprev[0] = x[ENR];
  for (i = 0; i < comps; i++)
    rad_data->Rprev[i + 1] = x[FR1 + i];

  // /*-- Update previous additional variables if needed --*/
  // #if RADIATION_FULL_CONVERGENCE == YES
  //   er1 = fabs(x[PRS] - rad_data->exv_prev) / fabs(rad_data->exv_prev);
  //   err += er1 * er1;
  //   rad_data->exv_prev = x[PRS];
  // #endif

  return err;
}

void RadStep(int cell, double dt, grid_t &grid)
/*!
 * Perform the implicit step along a direction determined during RadStep3D.
 * Update both primitive and conserved quantities. Primitive and conserved
 * fields must be consistent before every call to this function.
 *
 * \param [in,out]	uprim	  array of primitive variables
 * \param [in,out]	ucons	  array of conservative variables
 * \param [in,out]	source	array of source terms
 * \param [in]  	  ibeg	  starting index of computation
 * \param [in] 		  iend	 	final index of computation
 * \param [in,out] 	flag    array of flags
 * \param [in]      dt      time step
 *
 *********************************************************************** */
{
  int i, j, m;
  double err, gamma;

  static int comps = 4 - 1;
  // static Rad_data rad_data;

  static double *x, *dx, *mf, **J;

  if (rad_data.Ttot == NULL) {
    rad_data.Ttot = ARRAY_1D(RADIATION_NEQS, double);
    rad_data.Rini = ARRAY_1D(RADIATION_NEQS, double);
    rad_data.Rprev = ARRAY_1D(RADIATION_NEQS, double);
    rad_data.u = ARRAY_1D(3, double);

    x = ARRAY_1D(RADIATION_NEQS, double);
    dx = ARRAY_1D(RADIATION_NEQS, double);
    mf = ARRAY_1D(RADIATION_NEQS, double);
    J = ARRAY_2D(RADIATION_NEQS, RADIATION_NEQS, double);
  }

  /* ----------------------------
       Main loop on positions
     ---------------------------- */

  /*-- Store current position --*/
  rad_data.pos = i;

  /*-- Set initial energies --*/

  // Store (Erad,Frad) in Rini and Rprev
  rad_data.Rini[0] = rad_data.Ttot[0] = consvar[ENR];
  rad_data.Ttot[0] += consvar[ENG];

  rad_data.Rprev[0] = consvar[ENR];

  /*-- Set initial extra variables if full convergence is imposed --*/
  rad_data.exv_prev = primvar[PRS];

  /*-- Set fluxes and total momentum --*/
  for (j = 0; j < comps; j++) {

    // Store (Erad,Frad) in Rini and Rprev
    rad_data.Rini[j + 1] = rad_data.Ttot[j + 1] = consvar[FR1 + j];
    rad_data.Ttot[j + 1] += consvar[MX1 + j];

    rad_data.Rprev[j + 1] = consvar[FR1 + j];
  }

  // /*-- Initial guess for the iterated fields --*/
  // #if RADIATION_IMPL == RADIATION_NEWTON_GAS || RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
  //     gamma = consvar[RHO] / primvar[RHO];
  //     x[0] = primvar[PRS];
  // #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
  //     rad_data.Rprev[0] = x[0]; // Store gas pressure in Rprev
  // #endif
  //     for (j = 0; j < comps; j++) {
  //       x[j + 1] = gamma * primvar[VX1 + j];
  // #if RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS
  //       rad_data.Rprev[j + 1] = x[j + 1]; // Store proper velocity in Rprev
  // #endif
  //     }
  // #elif RADIATION_IMPL == RADIATION_NEWTON_RAD
  //     x[0] = consvar[ENR];
  //     for (j = 0; j < comps; j++)
  //       x[j + 1] = consvar[FR1 + j];
  // #endif

  /* -----------------------------------------------------------------
      Implicit step until convergence or maximum number of iterations
     ----------------------------------------------------------------- */
  m = 0;
  err = 1.0;
  while (err > 1e-7 && m++ < 200) {

    /*********************************
                Update iterated fields
     *********************************/
#if RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD

    /*-- Set coefficients of the system C.(E,F^i)^{(m+1)} == b  --*/
    RadFPMatrices(&rad_data, dx, J);

    /*-- Update iterated fields and store them in x --*/
    if (GaussianSolve(J, dx, x, 4)) {

      D_LD;
    }

    /*-- Compute -F and the Jacobian --*/
    RadNewtonJacobian(x, mf, J, &rad_data);

    /*-- Solve the system J*dx == -F --*/
    if (GaussianSolve(J, mf, dx, RADIATION_NEQS)) {
      WARNING(Where(i, NULL);)
      for (j = 0; j < RADIATION_NEQS; j++)
        dx[j] = mf[j];
      //  QUIT_PLUTO(1);
    }

    /*-- Update iterated fields --*/
    for (j = 0; j < RADIATION_NEQS; j++)
      x[j] += dx[j];

      /*********************************
                     Error calculation
       *********************************/
#if RADIATION_IMPL == RADIATION_FIXEDPOINT_RAD

    /*-- Update conserved variables --*/
    consvar[ENR] = x[0];
    consvar[ENG] = rad_data.Ttot[0] - x[0];
    for (j = 0; j < comps; j++) {
      consvar[FR1 + j] = x[j + 1];
      consvar[MX1 + j] = rad_data.Ttot[j + 1] - x[j + 1];
    }

    /*-- Update primitive variables --*/
    ConsToPrim(ucons, uprim, i, i, flag);

    /*-- Compute relative differences --*/
    err = RadErr(primvar, NULL, &rad_data);

#elif RADIATION_IMPL == RADIATION_FIXEDPOINT_GAS

    /*-- Compute relative differences --*/
    err = RadErr(x, NULL, &rad_data);

#else

    /*-- Compute relative differences --*/
    err = RadErr(x, dx, &rad_data);

#endif

  } // End of iterations

  /*-- Check number of iterations --*/
  if (m > 200) {
    D_L;
  }

/*-- Compute and store source terms if needed --*/
#if RADIATION_IMEX_SSP2 == YES
  RadSourceFunction(primvar, source[i]);
#endif

#endif
} // End of loop on positions

#endif

int rhllc::RadRHDTest() {
  WRITE_LOG("Start RadRHDTest()\n");
  grid_t grid;

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), 0);

  DIE_IF(rhllc::Init(glb_files.hllc_init_value, grid.cells));

  grid_directions_t grid_direction;
  files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);
  cuda::interface::InitDevice(glb_files.base_address, grid_direction, grid);

#pragma omp parallel for
  for (int i = 0; i < grid.size; i++) {
    if (grid.cells[i].geo.center[0] < 0.5) {
      grid.energy[i] = 1;
      grid.stream[i] = Vector3(0.1, 0, 0);
      grid.impuls[i] = Matrix3::Zero();
      grid.impuls[i](0, 0) = grid.impuls[i](1, 1) = grid.impuls[i](2, 2) = grid.energy[i];
    } else {
      grid.energy[i] = 0;
      grid.stream[i] = Vector3(0, 0, 0);
      grid.impuls[i] = Matrix3::Zero();
      grid.impuls[i](0, 0) = grid.impuls[i](1, 1) = grid.impuls[i](2, 2) = grid.energy[i];
    }
    grid.cells[i].phys_val.d = 1;
    grid.cells[i].phys_val.p = 10;
    grid.cells[i].phys_val.v = Vector3::Zero();
  }

  HllcPhysToConv(grid.cells);

  Type t = 0.0;
  Type cur_timer = 0;
  int res_count = _solve_mode.start_point;

  _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);

  WRITE_LOG("tau = %lf\n", _hllc_cfg.tau);

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid); // начальное сохранение

  auto start_clock = tick::steady_clock::now();

  while (t < _hllc_cfg.T) {
    Hllc3dStab(_hllc_cfg.tau, grid);

    if (cur_timer >= _hllc_cfg.save_timer) {
      DIE_IF(files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid) != e_completion_success);

      WRITE_LOG("t= %lf, step= %d, time_step=%lf\n", t, res_count, (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
      cur_timer = 0;
    }

    constexpr Type alpha = 1;
    constexpr Type betta = 1;
    constexpr Type temperature4 = 1;

#pragma omp parallel for
    for (int cell = 0; cell < grid.size; cell++) {
#if 0      
      Matrix4 T;
      Vector4 G;

      const Vector3 &v = grid.cells[cell].phys_val.v;
      Type Gamma = 1. / sqrt(1 - v.dot(v));
      Vector4 u(Gamma, Gamma * v[0], Gamma * v[1], Gamma * v[2]);

      T(0, 0) = grid.energy[cell];
      T(1, 0) = T(0, 1) = grid.stream[cell][0];
      T(2, 0) = T(0, 2) = grid.stream[cell][1];
      T(3, 0) = T(0, 3) = grid.stream[cell][2];
      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++)
          T(j + 1, k + 1) = grid.impuls[cell](j, k);
      }

      Type d = grid.cells[cell].phys_val.d;
      for (size_t m = 0; m < 4; m++) {
        G[m] = 0;
        Type TT = 0;
        for (size_t a = 0; a < 4; a++) {
          for (size_t b = 0; b < 4; b++) {
            TT += T(a, b) * u[a] * u[b] * u[m];
          }
          G[m] += -alpha * d * (T(m, a) * u[a] + 4 /** kStefanBoltzmann * temperature4*/ * u[m]) - betta * d * (T(m, a) * u[a]);
        }
        G[m] += betta * d * (-TT);
      }

      // G << grid.energy[cell], grid.stream[cell][0], grid.stream[cell][1], grid.stream[cell][2];
#else
      Vector4 G;
      RadSourceFunction(cell, grid, G);
      /*
      u[ENG]   += dt*S[0] ;
         u[ENR]   -= dt*S[0] ;
         for (n=0; n<comps; n++ )  {
             u[MX1+n] += dt*S[n+1] ;
             u[FR1+n] -= dt*S[n+1] ;
         }
     */
#endif

      constexpr Type ds = 1;
      grid.cells[cell].conv_val.p += ds * _hllc_cfg.tau * G[0];
      grid.cells[cell].conv_val.v[0] += ds * _hllc_cfg.tau * G[1];
      grid.cells[cell].conv_val.v[1] += ds * _hllc_cfg.tau * G[2];
      grid.cells[cell].conv_val.v[2] += ds * _hllc_cfg.tau * G[3];

      GetPhysValueStab(grid.cells[cell].conv_val, grid.cells[cell].phys_val);
    }

    t += _hllc_cfg.tau;
    cur_timer += _hllc_cfg.tau;
    _hllc_cfg.tau = GetTimeStep(_hllc_cfg, grid.cells);
  }

  files_sys::bin::WriteSolution(glb_files.solve_address + std::to_string(res_count++), grid);

  WRITE_LOG("End RadRHDTest()\n");
  return e_completion_success;
}

#endif //! SOLVERS