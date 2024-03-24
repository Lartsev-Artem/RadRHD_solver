/**
 * @file cuda_scattering.h
 * @brief Функция расчета интеграла рассеяния
 */
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

__global__ void GetS_MPI_Stream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid, const IdType start, const IdType end);

/**
 * @brief Расчет интеграла рассеяния с mpi разделением по направлениям и gpu разделением по ячейкам
 *
 * @param[in] dir сетка направлений
 * @param[inout] grid пространственная сетка(определены только массивы значений и размеры)
 * @param[in] size_loc локальное число ячеек для текущего ядра
 * @param[in] start_dir начало локальных направлений
 * @param[in] end_dir конец локальных направлений
 */
__global__ void GetS_MPI_multi_device(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid,
                                      const IdType size_loc, const IdType start_dir, const IdType end_dir);

__global__ void Get_spectrum_multi_device(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid,
                                          const Type frq,
                                          const IdType end_dir);

__global__ void Get_full_spectrum_multi_device(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid,
                                               const IdType end_dir);

} // namespace kernel
} // namespace cuda
#endif //! CUDA_SCATTERING_H
