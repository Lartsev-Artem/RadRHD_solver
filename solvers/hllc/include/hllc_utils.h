/**
 * @file hllc_utils.h
 * @brief Вспомогательные функции
 *
 */
#if !defined HLLC_UTILS_H && defined SOLVERS && defined HLLC
#define HLLC_UTILS_H

#include "solvers_struct.h"

/*! \addtogroup hllc Модуль расчета газовой динамики
    @{
*/

namespace hllc {

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

/**
 * @brief Функция расчитывает физическое состояние на всей сетки по консервативным
 *
 * @note Консервативные переменные должно быть записаны
 * @param[inout] cells ячейки сетки
 */
void HllcConvToPhys(std::vector<elem_t> &cells);

/**
 * @brief Функция задаёт условия контакта ячеек
 * @details Устанавливает или граничные условия или состояние соседней ячейки
 * @param[in] f грань
 * @param[in] cells ячейки сетки
 * @param[out] bound соседняя ячейка к текущей
 */
void BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound);

} // namespace hllc

#endif //! HLLC_UTILS_H