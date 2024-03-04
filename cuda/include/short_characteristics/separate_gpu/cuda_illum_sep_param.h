/**
 * @file cuda_illum_sep_param.h
 * @brief Функции расчёта параметров зависящих от излучения c разделением по видеокартам
 * \note отличается порядок хранения+ излучения хранится на ячейках.
 * Вычисление дивергенций в явном виде недоступны
 *
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
__global__ void MakeIllumParam(const cuda::geo::grid_directions_device_t *dir, cuda::geo::grid_device_t *grid, IdType size_params, IdType shift_params);

/**
 * @brief Устанавливает локальные сдвиги видеокарт, в случае симуляции нескольких видеокарт на одной
 *
 * @param grid
 * @param loc_size_gpu
 * @param shift_gpu
 * @param loc_size_params
 * @param shift_params
 * @return __global__
 */
__global__ void SetImDevice(geo::grid_device_t *grid,
                            const IdType loc_size_gpu,
                            const IdType shift_gpu,
                            const IdType loc_size_params,
                            const IdType shift_params);

} // namespace kernel

} // namespace separate_device
} // namespace cuda
#endif //! CUDA_ILLUM_SEP_PARAM_H