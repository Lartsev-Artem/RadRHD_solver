/**
 * @file rhllc_3d.h
 * @brief Функции расчета газодинамического состояние на текущий момент времени
 *
 */
#if !defined RHLLC_3D_H && defined RHLLC
#define RHLLC_3D_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

/**
 * @brief Пространство имён модуля газовой динамики в релятивистской постановке
 *
 */
namespace rhllc
{

    /**
     * @brief Газодинамический расчёт на текущий момент времени с шагом tau
     *
     * @details PLUTO архитектура
     * @param[in] tau шаг интегрирования
     * @param[inout] grid сетка
     */
    void Hllc3d(const Type tau, grid_t &grid);

} // namespace rhllc
#endif //! RHLLC_3D_H