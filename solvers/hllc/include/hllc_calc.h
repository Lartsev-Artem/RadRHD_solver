#if !defined HLLC_CALC_H && defined SOLVERS && defined HLLC
#define HLLC_CALC_H

#include "solvers_struct.h"

/*! \addtogroup hllc Модуль расчета газовой динамики
    @{
*/

namespace hllc {
/**
 * @brief Газодинамический расчёт на текущий момент времени с шагом tau
 *
 * @param[in] tau шаг интегрирования
 * @param[inout] grid сетка
 */
void Hllc3d(const Type tau, grid_t &grid);
} // namespace hllc

#endif //! HLLC_CALC_H