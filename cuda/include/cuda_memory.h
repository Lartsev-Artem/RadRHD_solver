#if !defined CUDA_MEMORY_H && defined USE_CUDA
#define CUDA_MEMORY_H

#include "dbgdef.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace cuda {
namespace mem_protected {

template <typename dataType, typename sizeT>
void inline Malloc(sizeT size, dataType **data) {
  if (CheckError(cudaMalloc((void **)data, size))) {
    EXIT_ERR("Error cudaMalloc\n");
  }
}

template <typename dataType>
void inline MallocHost(size_t size, dataType **data) {
  if (CheckError(cudaMallocHost((void **)data, size))) {
    EXIT_ERR("Error MallocHost\n");
  }
}

template <typename distType, typename srcType>
void inline CpyToDevice(distType *dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpy(dist, src, size, cudaMemcpyHostToDevice))) {
    EXIT_ERR("Error CpyToDevice\n");
  }
}

template <typename distType, typename srcType>
void inline CpyToHost(distType **dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpy(dist, src, size, cudaMemcpyDeviceToHost))) {
    EXIT_ERR("Error CpyToHost\n");
  }
}

template <typename distType, typename srcType>
void inline CpyToDeviceAsync(distType *dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyHostToDevice))) {
    EXIT_ERR("Error CpyToDeviceAsync\n");
  }
}

template <typename distType, typename srcType>
void inline CpyToHostAsync(distType **dist, const srcType *src, size_t size) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyDeviceToHost))) {
    EXIT_ERR("Error CpyToHostAsync\n");
  }
}

template <typename distType, typename srcType>
void inline CpyToHostAsyncStream(distType **dist, const srcType *src, size_t size, cudaStream_t &st) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyDeviceToHost, st))) {
    EXIT_ERR("Error CpyToHostAsyncStream\n");
  }
}

template <typename dataType>
void inline FreeMem(dataType *&data) {
  if (CheckError(cudaFree(data))) {
    EXIT_ERR("Error FreeMem\n");
  }
}
template <typename dataType>
void inline FreeMemHost(dataType *&data) {
  if (CheckError(cudaFreeHost(data))) {
    EXIT_ERR("Error cudaFreeHost\n");
  }
}

} // namespace mem_protected
} // namespace cuda

#define CUDA_MEMCPY_TO_DEVICE(dist, src, size) CUDA_CALL_FUNC(cudaMemcpy, dist, src, size, cudaMemcpyHostToDevice);
#define CUDA_MEMCPY_TO_HOST(dist, src, size) CUDA_CALL_FUNC(cudaMemcpy, dist, src, size, cudaMemcpyDeviceToHost);

#define CUDA_MEMCPY_TO_DEVICE_ASYNC(dist, src, size) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyHostToDevice);
#define CUDA_MEMCPY_TO_HOST_ASYNC(dist, src, size) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyDeviceToHost);
#define CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(dist, src, size, st) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyDeviceToHost, st);

#define CUDA_MALLOC(src, size) CUDA_CALL_FUNC(cudaMalloc, (void **)src, size);
#define CUDA_MALLOC_HOST(src, size) CUDA_CALL_FUNC(cudaMallocHost, (void **)src, size);

#define CUDA_FREE_MEMORY(val) CUDA_CALL_FUNC(cudaFree, val)
#define CUDA_FREE_HOST_MEMORY(val) CUDA_CALL_FUNC(cudaFreeHost, val)

#endif // CUDA_MEMORY_H