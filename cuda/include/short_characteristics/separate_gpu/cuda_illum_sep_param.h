/**
 * @file cuda_illum_param.h
 * @brief Функции расчёта параметров зависящих от излучения

 */
#if !defined CUDA_ILLUM_SEP_PARAM_H && defined USE_CUDA
#define CUDA_ILLUM_SEP_PARAM_H

#include "cuda_struct.h"

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    @{
*/

namespace cuda {

namespace separate_device {

namespace device {
#ifdef ON_FULL_ILLUM_ARRAYS

/**
 * @brief Функция расчёта энергии излучения на сетке
 *
 * @param[in] dir сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение
 */
__device__ void MakeEnergy(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

/**
 * @brief Функция расчёта потока излучения на сетке
 *
 * @param[in] dir сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение
 */
__device__ void MakeStream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

/**
 * @brief  Функция расчёта импульса излучения на сетке
 *
 * @param[in] dir сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение
 */
__device__ void MakeImpuls(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

#endif

} // namespace device

namespace kernel {

/**
 * @brief Глобальная функция ядра для вычисления параметров излучения
 *
 * @param[in] dir сфера направлений
 * @param[inout] grid сетка
 * @return __global__ use cudaGetLastError
 */
__global__ void MakeIllumParam(const cuda::geo::grid_directions_device_t *dir, cuda::geo::grid_device_t *grid);

} // namespace kernel

} // namespace separate_device
} // namespace cuda
#endif //! CUDA_ILLUM_SEP_PARAM_H