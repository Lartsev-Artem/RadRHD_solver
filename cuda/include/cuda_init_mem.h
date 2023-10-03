/**
 * @file cuda_init_mem.h
 * @brief Функции инициализации и очистки памяти на видеокарте
 *
 */
#if !defined CUDA_INIT_MEM_H && defined USE_CUDA
#define CUDA_INIT_MEM_H

#include "cuda_struct.h"
#include "geo_types.h"
#include "solvers_struct.h"

namespace cuda {

extern geo::device_host_ptr_t device_host_ptr; ///< связь хоста с массивами внутри видеокарты (через эту структуры идёт обращение к массивам на карте)

/**
 * @brief Функция инициализирует память и переносит структуру сферы направлений на видеокарту
 *
 * @param[in] grid_host сфера направлений
 * @param[out] grid_device указатель на сферу направлений на видеокарте
 */
void InitDirectionsOnDevice(const grid_directions_t &grid_host, geo::grid_directions_device_t *&grid_device);

/**
 * @brief Функция очищает память на видеокарте
 *
 * @param[in] grid_device указатель на сферу направлений на видеокарте
 */
void ClearDirectionsOnDevice(geo::grid_directions_device_t *&grid_device);

/**
 * @brief Функция инициализирует память под сетку на видеокарте
 *
 * @param[in] grid_host сетка на хосте
 * @param[in] size_dir число направлений
 * @param[in] start_dir номер начала направлений на узле (mpi)
 * @param[in] end_dir номер конца направлений на узле (mpi)
 * @param[out] grid_device указатель на структуру сетки на карте
 */
void InitGridOnDevice(const grid_t grid_host,
                      const int size_dir, const int start_dir, const int end_dir,
                      cuda::geo::grid_device_t *&grid_device);

/**
 * @brief  Функция очищает память на видеокарте
 *
 * @param[in] grid_device  указатель на структуру сетки на карте
 */
void ClearGridOnDevice(cuda::geo::grid_device_t *&grid_device);
} // namespace cuda
#endif //! CUDA_INIT_MEM_H