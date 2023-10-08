/**
 * @file rhllc_utils.h
 * @brief Вспомогательные функции
 */
#if !defined RHLLC_UTILS_H && defined SOLVERS
#define RHLLC_UTILS_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {
/**
 * @brief Функция задаёт условия контакта ячеек
 * @details Устанавливает или граничные условия или состояние соседней ячейки
 * @param[in] f грань
 * @param[in] cells ячейки сетки
 * @param[out] bound соседняя ячейка к текущей
 */
void BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound);

/**
 * @brief Расчитывает текущий шаг интегрирования с учётом геометрии и настроек решателя
 *
 * @param[in] hllc_set настройки hllc решателя
 * @param[in] cells ячейки сетки
 * @return Type шаг интегрирования tau
 */
Type GetTimeStep(const hllc_value_t &hllc_set, const std::vector<elem_t> &cells);

/**
 * @brief Функция расчитывает консервативное состояние на всей сетки по физическому
 *
 * @note Физические переменные должно быть записаны
 * @param[inout] cells ячейки сетки
 */
void HllcPhysToConv(std::vector<elem_t> &cells);
} // namespace rhllc
#endif //! RHLLC_UTILS_H