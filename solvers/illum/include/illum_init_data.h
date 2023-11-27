/**
 * @file illum_init_data.h
 * @brief инициализация состояния излучения

 */
#if !defined ILLUM_INIT_DATA_H && defined SOLVERS && defined ILLUM
#define ILLUM_INIT_DATA_H

#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    @{
*/

namespace illum {

/**
 * @brief Функция устанавливает начальное состояние радиационных переменных
 *
 * @param[in] address_data адрес с бинарными данными
 * @param[inout] grid сетка
 * @return int ::e_type_completion
 */
int InitRadiationState(const std::string &address_data, grid_t &grid);

} // namespace illum
#endif //! ILLUM_INIT_DATA_H