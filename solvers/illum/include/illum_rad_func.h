#if 0 // !defined ILLUM_RADFUNC_H && defined ILLUM
#define ILLUM_RADFUNC_H

#include "geo_types.h"
#include "global_value.h"

/**
  @note
 * rho=m_h*n (плотность с концентрацией)
   phs = nkT (давление с концетрацией)
 */
namespace illum {

/**
 * @brief Приближение к интегралу планка
 *
 * @tparam T float-point-type
 * @param[in] x
 * @return значение интеграла
 * @note \c Четверушкин Б.Н. Моделирование задач РГД. с. 47
 */
template <typename T>
inline T f_plank(T x) {

#ifdef DEBUG
  if (x > 700) {
    WRITE_LOG_ERR("loss of accuracy exp(%lf)\n", x);
  }
#endif
  T x2 = x * x;
  if (x <= 2)
    return (x2 * x * ((1. / 3.) - x / 8. + x2 / 62.4));
  else
    return (PI * PI * PI * PI / 15) - exp(-x) * (x2 * x + 3 * x2 + 6 * x + 7.28);
}

/**
 * @brief Частотная составляющая константы стефана больцмана
 * @note SigmaSB(0, Infinity)=k_stefan_boltzmann
 *
 * @tparam T  float-point-type
 * @param [in]	temperature	Input temperature.
 * @param [in]  frq0 Input left frequency
 * @param [in]  frq1 Input right frequency
 * @return приближение к интегральной постоянной в диапазоне заданных частот
 */
template <typename T>
inline T SigmaSB(T temperature, T frq0, T frq1) {

  constexpr T a_rad = (2 * PI * k_boltzmann * k_boltzmann * k_boltzmann * k_boltzmann) / (kC_Light * kC_Light * kH_plank * kH_plank * kH_plank);
  constexpr T arg_c = kH_plank / k_boltzmann;
  const T A = arg_c / temperature;
  return a_rad * (f_plank(A * frq1) - f_plank(A * frq0));
}

/**
 * @brief Функция планка для частотного-группового приближения
 *
 * @warning предел частот для тела с температурой 1000К - начало ультрафиолетового спектра, т.к. (exp(-x0)-exp(-x1)) числено вырождается в ноль
 * также сильное дробление спектра приведет к неустойчивому поведению. Допустимый шаг 50нм
 * @tparam T float-point-type
 * @param [in]	temperature	Input temperature.
 * @param [in]  frq0 Input left frequency
 * @param [in]  frq1 Input right frequency
 * @return безразмерная интенсивность
 */
template <typename T>
inline T B_Plank(T temperature, T frq0, T frq1) {

  const T sigma = SigmaSB(temperature, frq0, frq1) / PI; //* (kTime * kTime * kTime / kMass); //обезразмеривание
  T T4 = temperature * temperature;
  return sigma * (T4 * T4);
}

/**
 * @brief Функция планка для абсолютно черного тела
 *
 * @tparam T float-point-type
 * @param [in]	temperature	Input temperature.
 * @return T безразмерная интенсивность
 */
template <typename T>
inline T B_Plank(T temperature) {
  constexpr T a_rad = kStefanBoltzmann / PI * (kTime * kTime * kTime / kMass);
  T T4 = temperature * temperature;
  return a_rad * (T4 * T4);
}

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
  // constexpr T g_idealGasConst = 1; //PLUTO
  // constexpr T g_idealGasConst = (kM_hydrogen / k_boltzmann); // * (kDist * kDist) / (kTime * kTime);
  constexpr T g_idealGasConst = (kM_hydrogen / k_boltzmann) * kPressure / kDensity; // * (kDist * kDist) / (kTime * kTime);
  return g_idealGasConst * prs / rho;
}
}; // namespace illum

#endif //! ILLUM_RADFUNC_H