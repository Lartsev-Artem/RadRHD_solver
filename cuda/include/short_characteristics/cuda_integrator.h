/**
 * @file cuda_integrator.h
 * @brief Функции интегрирования по направлениям
 */
#if !defined CUDA_INTEGRATOR_H && defined USE_CUDA
#define CUDA_INTEGRATOR_H

#include "cuda_struct.h"

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    @{
*/

namespace cuda::device {

/**
 * @brief Пространство имён функций интегрирования по направлениям
 *
 */
namespace direction_integrator {

/**
 * @brief Функция расчёта индикатрисы рассеяния на электронах для двух направлений
 *
 * @param[in] direction направление 1
 * @param[in] direction2 направление 2
 * @return  значение индикатрисы
 */
__device__ Type Gamma(const Vector3 &direction, const Vector3 &direction2);

/**
 * @brief Функция интегрирования по направлениям, усреднённая по ячейке (плотность энергии излучения)
 *
 * @param[in] num_cell номер ячейки
 * @param[in] dir сфера направлений
 * @param[in] grid сетка
 * @return  значение интеграла (энергия на ячейке)
 */
__device__ Type IntegrateByCell(const int num_cell, const geo::grid_directions_device_t *dir, const geo::grid_device_t *grid);

/**
 * @brief Функция интегрирования по направлениям с вектором направления на каждой грани (поток излучения)
 *
 * @param[in] num_cell  номер ячейки
 * @param[in] dir_grid сфера направлений
 * @param[in] grid сетка
 * @param[out] Stream массив векторов длины ::CELL_SIZE(поток на грани
 */
__device__ void IntegrateByFaces3(const int num_cell, const geo::grid_directions_device_t *dir_grid, geo::grid_device_t *grid, Vector3 *Stream);

/**
 * @brief Функция интегрирования по направлениям с тензором направления на каждой грани (импульс излучения)
 *
 * @param num_cell номер ячейки
 * @param dir_grid сфера направлений
 * @param grid сетка
 * @param Impuls массив матриц длины ::CELL_SIZE(импульс на грани)
 */
__device__ void IntegrateByFaces9(const int num_cell, const geo::grid_directions_device_t *dir_grid, geo::grid_device_t *grid, Matrix3 *Impuls);

} // namespace direction_integrator
} // namespace cuda::device
#endif // CUDA_INTEGRATOR_H