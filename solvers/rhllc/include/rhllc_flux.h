/**
 * @file rhllc_flux.h
 * @brief Функции расчёта потоков
 *
 */
#if !defined RHLLC_FLUX_H && defined SOLVERS
#define RHLLC_FLUX_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {

/**
 * @brief Расчёт консервативных переменных через физические
 *
 * @details статья ost1098
 * @param[in] W физические переменные
 * @param[out] U консервативные переменные
 */
void GetConvValue(const flux_t &W, flux_t &U);

/**
 * @brief Восстановление физических переменных через консервативные
 *
 * @details алгоритм построен на основе метода Ньютона. из статьи ost1098
 * @param[in] U консервативные переменные
 * @param[out] W физические переменные
 */
void GetPhysValue(const flux_t &U, flux_t &W);

/**
 * @brief Восстановление физических переменных через консервативные
 *
 * @details алгоритм построен на основе метода Ньютона. из статьи ost1098
 * @param[in] U консервативные переменные
 * @param[out] W физические переменные
 * @return флаг ошибки (ограничение по давлению/плотности)
 * @note алгоритм содержит в себе ограничение на минимальное значения давления и плотности
 */
bool GetPhysValueSave(const flux_t &U, flux_t &W);

/**
 * @brief Расчёт потоков методом HLLC
 *
 * @details An HLLC Riemann solver for relativistic flows – I. Hydrodynamics A. Mignone and G. Bodo
    INAF Osservatorio Astronomico di Torino, 10025 Pino Torinese, Italy
    Accepted 2005 August 22. Received 2005 August 16; in original form 2005 May 16
    Mon. Not. R. Astron. Soc. 364, 126–136 (2005)

    https://github.com/PrincetonUniversity/Athena-Cversion/blob/master/src/rsolvers/hllc_sr.c
 *
 * @note комбинация кода из mignone 2005 и 2006. в части hllc
 * @param[in] conv_val_l консервативный ячейка слева
 * @param[in] conv_val_r консервативный ячейка справа
 * @param[in] phys_val_l физическая ячейка слева
 * @param[in] phys_val_r физическая ячейка справа
 * @param[out] f поток
 */
void GetFlux(const flux_t &conv_val_l, const flux_t &conv_val_r,
             const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f);

} // namespace rhllc

#endif //! RHLLC_FLUX_H