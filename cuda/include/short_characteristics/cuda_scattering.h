#if !defined CUDA_SCATTERING_H && defined USE_CUDA
#define CUDA_SCATTERING_H

#include "cuda_struct.h"

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    @{
*/

namespace cuda {

/**
 * @brief Пространство имён с глобальными функциями ядра
 *
 */
namespace kernel {

/**
 * @brief Расчёт интеграла рассеяния с одним потоком
 *
 * @param dir сфера направлений
 * @param grid сетка
 * @return __global__ use cudaGetLastError
 */
__global__ void GetS(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

__global__ void GetS_MPI(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

__global__ void GetS_MPI_Stream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid, const int start, const int end);

} // namespace kernel
} // namespace cuda
#endif //! CUDA_SCATTERING_H