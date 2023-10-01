/**
 * @file illum_integrate.h
 * @brief Функций интегрирования по сфере направлений
 *
 */
#if !defined ILLUM_INTEGRATE_H && defined SOLVERS && defined ILLUM
#define ILLUM_INTEGRATE_H

#include "geo_types.h"

namespace illum {

/**
 * @brief Пространство имён для функций интегрирования по сфере направлений
 *
 */
namespace direction_integrator {

/**
 * @brief Функция интегрирования по направлениям, усреднённая по ячейке (плотность энергии излучения)
 *
 * @param[in] Illum  излучение в области, заданная на гранях
 * @param[in] grid_direction сфера направлений
 * @return значение интеграла (энергия на ячейке)
 */
Type IntegrateByCell(const std::vector<Type> &Illum, const grid_directions_t &grid_direction);

/**
 * @brief Функция интегрирования по направлениям с вектором направления на каждой грани (поток излучения)
 *
 * @param[in] Illum излучение в области, заданная на гранях
 * @param[in] grid_direction сфера направлений
 * @param[out] stream_face массив векторов длины ::CELL_SIZE(поток на грани
 */
void IntegrateByFace3(const std::vector<Type> &Illum, const grid_directions_t &grid_direction, Vector3 *stream_face);

/**
 * @brief @brief Функция интегрирования по направлениям с тензором направления на каждой грани (импульс излучения)
 *
 * @param[in] Illum излучение в области, заданная на гранях
 * @param[in] grid_direction сфера направлений
 * @param[out] impuls_face массив матриц длины ::CELL_SIZE(импульс на грани)
 */
void IntegrateByFace9(const std::vector<Type> &Illum, const grid_directions_t &grid_direction, Matrix3 *impuls_face);

} // namespace direction_integrator
} // namespace illum

#endif //! ILLUM_INTEGRATE_H