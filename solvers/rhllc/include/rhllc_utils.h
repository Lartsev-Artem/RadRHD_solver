/**
 * @file rhllc_utils.h
 * @brief Вспомогательные функции
 */
#if !defined RHLLC_UTILS_H && defined RHLLC
#define RHLLC_UTILS_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc
{

  extern Type max_signal_speed; ///< максимальная скорость сигнала

  /**
   * @brief Расчитывает текущий шаг интегрирования с учётом геометрии и настроек решателя
   *
   * @param[in] hllc_set настройки hllc решателя
   * @param[in] cells ячейки сетки
   * @return Type шаг интегрирования tau
   */
  inline Type GetTimeStep(const hllc_value_t &hllc_set)
  {
    const Type t = hllc_set.CFL * hllc_set.h_min / max_signal_speed;
    RRHD_ASSERT(t < 0);
    return t;
  }

  /**
   * @brief Функция расчитывает консервативное состояние на всей сетки по физическому
   *
   * @note Физические переменные должно быть записаны
   * @param[inout] grid сетка
   */
  void HllcPhysToConv(grid_t &grid);

  /**
   * @brief Функция расчитывает физическое состояние на всей сетки по консервативному
   *
   * @note консервативные переменные должно быть записаны
   * @param[inout] grid сетка
   */
  void HllcConvToPhys(grid_t &grid);

#ifdef RRHD_DEBUG
  /**
   * @brief Проверка условия физичности консервативных переменных
   *
   * @param U
   * @return true if bad
   */
  static inline bool CheckConvState(const flux_t &U)
  {
    if (U.p - sqrt(U.d * U.d + U.v.dot(U.v)) < 0)
    {
      return true;
    }
    return false;
  }
#endif

} // namespace rhllc
#endif //! RHLLC_UTILS_H