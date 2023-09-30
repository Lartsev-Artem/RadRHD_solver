#if !defined SCATTERING_H && defined SOLVERS && defined ILLUM
#define SCATTERING_H

#include "geo_types.h"
#include "solvers_struct.h"

namespace illum {
namespace scattering {

/**
 * @brief Функция расчёта индикатрисы рассеяния на электронах для двух направлений
 *
 * @details \f$ Gamma(\vec x, \vec\omega, \vec\omega^')= \dfrac34 \left(1+(\vec\omega, \vec\omega^')^2\right) \f$,  num_cell equals x
 * @param[in] direction направление 1
 * @param[in] direction2 направление 2
 * @return значение индикатрисы
 */
inline double Gamma(const Vector3 &direction, const Vector3 &direction2) {
  const Type dot = direction.dot(direction2);
  return (3. * (1 + dot * dot)) / 4.;
}

/**
 * @brief Функция рассчитывает интеграл рассеяния по указанному направлению я одной ячейке
 *
 * @details \f$ \int\limits_{\Omega}\Gamma(s, \vec\omega, \vec\omega^{'})I(s, \vec\omega^')d\vec\omega^' \f$
 * @param[in] num_cell_shift сдвиг данных для ячейки (num_cell*SIZE_CELLS, если значения на гранях)
 * @param[in] direction направление интегрирования
 * @param[in] size_one_direction размер данных излучения по одному направлению (обычно size_grid*SIZE_CELLS, т.к. излучения задано на гранях)
 * @param[in] illum_old массив излучения, расчитанный с предыдущей итерации
 * @param[in] grid_direction сфера направлений
 * @return значение интеграла
 */
Type GetIntScatteringByDirection(const int num_cell_shift, const Vector3 &direction, const int size_one_direction, const Type *illum_old,
                                 const grid_directions_t &grid_direction);

/**
 * @brief Функция расчёта интеграла рассеяния по всей сетке на cpu
 *
 * @param[in] grid_direction сфера направления
 * @param[inout] grid сетка
 */
void CalculateIntCPU(const grid_directions_t &grid_direction, grid_t grid);

} // namespace scattering
} // namespace illum

#endif //! SCATTERING_H