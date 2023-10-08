/**
 * @file cuda_illum_param.h
 * @brief Функции расчёта параметров зависящих от излучения

 */
#if !defined CUDA_ILLUM_PARAM_H && defined USE_CUDA
#define CUDA_ILLUM_PARAM_H

#include "cuda_struct.h"

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    @{
*/

namespace cuda {

/**
 * @brief Пространство имён функций определённых только на видеокарте
 *
 */
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
#endif

/**
 * @brief Функция расчёта потока и дивергенции потока излучения на сетке
 *
 * @param[in] dir сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение
 */
__device__ void MakeDivStream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

/**
 * @brief  Функция расчёта импульса и дивергенции импульса излучения на сетке
 *
 * @param[in] dir сфера направлений
 * @param[inout] grid сетка
 * @warning предварительно должно быть вычислено излучение
 */
__device__ void MakeDivImpuls(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);
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
} // namespace cuda
#endif //! CUDA_ILLUM_PARAM_H