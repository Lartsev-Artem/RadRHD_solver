/**
 * @file rhllc_calc.h
 * @brief Функции расчета газодинамического состояние на текущий момент времени
 *
 */
#if !defined RHLLC_CALC_H && defined SOLVERS
#define RHLLC_CALC_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {

/**
 * @brief Газодинамический расчёт на текущий момент времени с шагом tau
 *
 * @param[in] tau шаг интегрирования
 * @param[inout] grid сетка
 */
void Hllc3d(const Type tau, grid_t &grid);

/**
 * @brief Газодинамический расчёт на текущий момент времени с шагом tau
 *
 * @details PLUTO архитектура
 * @param[in] tau шаг интегрирования
 * @param[inout] grid сетка
 */
void Hllc3dStab(const Type tau, grid_t &grid);

} // namespace rhllc
#endif //! RHLLC_CALC_H