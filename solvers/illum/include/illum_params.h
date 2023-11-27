/**
 * @file illum_params.h
 * @brief Функции расчета параметров зависящих от излучения
 */
#if !defined ILLUM_PARAMS_H && defined SOLVERS && defined ILLUM
#define ILLUM_PARAMS_H

#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    @{
*/

namespace illum {

/**
 * @brief Функция расчёта энергии излучения на сетке
 *
 * @note не используйте эти функции совместно с Cuda
 * @param[in] grid_direction сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение и записано  в  ::elem_t
 */
void GetEnergy(const grid_directions_t &grid_direction, grid_t &grid);

/**
 * @brief Функция расчёта потока и дивергенции потока излучения на сетке
 *
 * @note не используйте эти функции совместно с Cuda
 * @param[in] grid_direction сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение и записано  в  ::elem_t
 */
void GetStream(const grid_directions_t &grid_direction, grid_t &grid);

/**
 * @brief Функция расчёта импульса и дивергенции импульса излучения на сетке
 *
 * @note не используйте эти функции совместно с Cuda
 * @param[in] grid_direction сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение и записано  в  ::elem_t
 */
void GetImpuls(const grid_directions_t &grid_direction, grid_t &grid);
} // namespace illum

#endif //! ILLUM_PARAMS_H