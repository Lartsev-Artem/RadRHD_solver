/**
 * @file radRHD_utils.h
 * @brief Some functions used by the radiation module.
 *
 */
#if !defined RADRHD_UTILS_H && defined RAD_RHD
#define RADRHD_UTILS_H

#include "geo_types.h"
#include "global_value.h"

namespace rad_rhd {

/**
 * @brief Compute the (i,j) component of the Eddington tensor
          using the M1 closure, taking as input the state v.
          (i,j) values must always be within FR1 and FR3.
 *
 * @param[in] i component
 * @param[in] j component
 * @param[in] U radiation energy
 * @param[in] F radiation stream
 * @return  (i,j) component of the Eddington tensor
 */
Type EddTensor(int i, int j, const Type U, const Vector3 &F);

/**
 * @brief Compute 4-vector radiation source
 *
 * @param[in] cell num cell
 * @param[in] grid grid with computed Ur, Fr, Tr
 * @param[out] G radiation source
 */
void GetRadSource(const int cell, const grid_t &grid, Vector4 &G);

// #include "illum_rad_func.h"
#if 0
/**
 * @brief  Compute the blackbody intensity corresponding to the input temperature.
 *
 * @tparam T float-point-type
 * @param [in]	temperature	Input temperature.
 * @return the blackbody intensity
 */
template <typename T>
inline T Blackbody(T temperature) {
  constexpr T g_radiationConst = 1; /// \todo PLUTO здесь kStefanBoltzmann / (PI*C_Light)??
  // constexpr T g_radiationConst = kStefanBoltzmann / PI * (kTime * kTime * kTime / kMass);
  T T4 = temperature * temperature;
  return g_radiationConst * (T4 * T4);
}

/**
 * @brief Compute the (ideal) gas temperature corresponding to given gas pressure
 *
 * @tparam T float-point-type
 * \param [in]	rho Input mass density.
 * \param [in]	prs Input (gas) pressure.
 * \note размерная величина!
 * @return gas temperature
 */
template <typename T>
inline T GetTemperature(T rho, T prs) {
  constexpr T g_idealGasConst = 1;
  // constexpr T g_idealGasConst = (kM_hydrogen / k_boltzmann) * (kDist * kDist) / (kTime * kTime);
  return g_idealGasConst * prs / rho;
}
#endif
}; // namespace rad_rhd

#endif //! RADRHD_UTILS_H