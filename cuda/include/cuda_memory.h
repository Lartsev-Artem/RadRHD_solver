/**
 * @file cuda_memory.h
 * @brief Функции безопасного обращения с памятью
 *
 */
#if !defined CUDA_MEMORY_H && defined USE_CUDA
#define CUDA_MEMORY_H

#include "cuda_def.h"
#include "dbgdef.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    @{
*/

namespace cuda {

/**
 * @brief Пространство имён безопасного обращения с памятью
 *
 */
namespace mem_protected {

#ifdef DEBUG
static double full_device_mem = 0;
static double full_host_mem = 0;
#endif

/**
 * @brief Выделение памяти на карте
 *
 * @tparam dataType тип данных
 * @tparam sizeT целочисленный тип
 * @param size размер в байтах
 * @param data указатель массив данных (передавать по адресу)
 */
template <typename dataType, typename sizeT>
void inline Malloc(sizeT size, dataType **data) {
#ifdef DEBUG
  double kb = (double)size / 1024.;
  full_device_mem += (kb / 1024.);
  WRITE_LOG("Malloc +%lf Kb. full=%lf Mb\n", kb, full_device_mem);
#endif
  if (CheckError(cudaMalloc((void **)data, size))) {
    EXIT_ERR("Error cudaMalloc\n");
  }
}

/**
 * @brief  Выделение памяти на хосте для быстрого копирования данных
 *
 * @tparam dataType тип данных
 * @tparam sizeT целочисленный тип
 * @param size размер в байтах
 * @param data указатель массив данных (передавать по адресу)
 */
template <typename dataType, typename sizeT>
void inline MallocHost(sizeT size, dataType **data) {
#ifdef DEBUG
  double kb = (double)size / 1024.;
  full_host_mem += (kb / 1024.);
  WRITE_LOG("MallocHost +%lf Kb, full=%lf Mb\n", kb, full_host_mem);
#endif
  if (CheckError(cudaMallocHost((void **)data, size))) {
    EXIT_ERR("Error MallocHost\n");
  }
}

/**
 * @brief Блокирующее копирование данных на карту
 *
 * @tparam distType тип данных на карте
 * @tparam srcType тип данных на хосте
 * @param dist массив на карте
 * @param src массив на хосте
 * @param size размер в байтах
 */
template <typename distType, typename srcType>
void inline CpyToDevice(distType *dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpy(dist, src, size, cudaMemcpyHostToDevice))) {
    EXIT_ERR("Error CpyToDevice\n");
  }
}

/**
 * @brief Блокирующее копирование данных на хост
 *
 * @tparam distType тип данных на хосте
 * @tparam srcType  тип данных на карте
 * @param dist массив на хосте
 * @param src массив на карте
 * @param size размер в байтах
 */
template <typename distType, typename srcType>
void inline CpyToHost(distType *dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpy(dist, src, size, cudaMemcpyDeviceToHost))) {
    EXIT_ERR("Error CpyToHost\n");
  }
}

/**
 * @brief Не блокирующее копирование данных на карту
 *
 * @tparam distType тип данных на карте
 * @tparam srcType тип данных на хосте
 * @param dist массив на карте
 * @param src массив на хосте
 * @param size размер в байтах
 */
template <typename distType, typename srcType>
void inline CpyToDeviceAsync(distType *dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyHostToDevice))) {
    EXIT_ERR("Error CpyToDeviceAsync\n");
  }
}

/**
 * @brief Не блокирующее копирование данных на хост
 *
 * @tparam distType тип данных на хосте
 * @tparam srcType  тип данных на карте
 * @param dist массив на хосте
 * @param src массив на карте
 * @param size размер в байтах
 */
template <typename distType, typename srcType>
void inline CpyToHostAsync(distType *dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyDeviceToHost))) {
    EXIT_ERR("Error CpyToHostAsync\n");
  }
}

/**
 * @brief Не блокирующее копирование данных на хост с разделением по потокам
 *
 * @tparam distType тип данных на хосте
 * @tparam srcType  тип данных на карте
 * @param dist массив на хосте
 * @param src массив на карте
 * @param size размер в байтах
 * @param st cuda поток
 */
template <typename distType, typename srcType>
void inline CpyToHostAsyncStream(distType *dist, const srcType *src, size_t size, cudaStream_t &st) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyDeviceToHost, st))) {
    EXIT_ERR("Error CpyToHostAsyncStream\n");
  }
}

/**
 * @brief Очистка памяти на карте
 *
 * @tparam dataType тип данных
 * @param data массив на карте
 */
template <typename dataType>
void inline FreeMem(dataType *&data) {
  if (!data) {
    return;
  }
  if (CheckError(cudaFree(data))) {
    EXIT_ERR("Error FreeMem\n");
  }
}

/**
 * @brief  Очистка памяти на хосте
 *
 * @tparam dataType тип данных
 * @param data массив на хосте
 */
template <typename dataType>
void inline FreeMemHost(dataType *&data) {
  if (CheckError(cudaFreeHost(data))) {
    EXIT_ERR("Error cudaFreeHost\n");
  }
}

} // namespace mem_protected
} // namespace cuda

#if 0 // устаревшие обращения
#define CUDA_MEMCPY_TO_DEVICE(dist, src, size) CUDA_CALL_FUNC(cudaMemcpy, dist, src, size, cudaMemcpyHostToDevice);
#define CUDA_MEMCPY_TO_HOST(dist, src, size) CUDA_CALL_FUNC(cudaMemcpy, dist, src, size, cudaMemcpyDeviceToHost);

#define CUDA_MEMCPY_TO_DEVICE_ASYNC(dist, src, size) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyHostToDevice);
#define CUDA_MEMCPY_TO_HOST_ASYNC(dist, src, size) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyDeviceToHost);
#define CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(dist, src, size, st) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyDeviceToHost, st);

#define CUDA_MALLOC(src, size) CUDA_CALL_FUNC(cudaMalloc, (void **)src, size);
#define CUDA_MALLOC_HOST(src, size) CUDA_CALL_FUNC(cudaMallocHost, (void **)src, size);

#define CUDA_FREE_MEMORY(val) CUDA_CALL_FUNC(cudaFree, val)
#define CUDA_FREE_HOST_MEMORY(val) CUDA_CALL_FUNC(cudaFreeHost, val)
#endif
#endif // CUDA_MEMORY_H