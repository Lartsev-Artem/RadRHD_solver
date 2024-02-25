/**
 * @file gas_state.h
 * @brief Функции связи параметров газа через уравнение состояния идеального газа
 * @version 0.1
 * @date 2024-02-13
 *
 */
#ifndef GAS_STATE_H
#define GAS_STATE_H

#include "global_consts.h"

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
  // constexpr T g_idealGasConst = 1; //PLUTO
  constexpr T g_idealGasConst = (kM_hydrogen / k_boltzmann) * kPressure / kDensity; // * (kDist * kDist) / (kTime * kTime);
  return g_idealGasConst * prs / rho;
}

/**
 * @brief Compute the (ideal) gas pressure corresponding to given gas temperature
 *
 * @tparam T float-point-type
 * \param [in]	rho Input mass density.
 * \param [in]	tem Input (gas) temperature. (размерная)
 * \note безразмерная величина!
 * @return gas pressure
 */
template <typename T>
inline T GetPressure(T rho, T tem) {
  constexpr T g_idealGasConst = (kM_hydrogen / k_boltzmann);
  constexpr T g_idealGasConst_inv = kDensity / g_idealGasConst;
  return (g_idealGasConst_inv / kPressure) * tem * rho;
}

/**
 * @brief Compute the (ideal) gas density corresponding to given gas temperature
 *
 * @tparam T float-point-type
 * \param [in]	prs Input mass pressure.
 * \param [in]	tem Input (gas) temperature.
 * \note безразмерная величина!
 * @return gas density
 */
template <typename T>
inline T GetDensity(T prs, T tem) {
  constexpr T g_idealGasConst = (kM_hydrogen / k_boltzmann);
  return (((g_idealGasConst * kPressure / kDensity) * prs) / tem);
}
#endif //! GAS_STATE_H