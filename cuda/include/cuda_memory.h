#if !defined CUDA_MEMORY_H && defined USE_CUDA
#define CUDA_MEMORY_H

#include "dbgdef.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace cuda {
namespace mem_protected {

template <typename distType, typename srcType, typename sizeT>
void inline CpyToDevice(distType *&dist, const srcType *src, sizeT size) {

  if (CheckError(cudaMemcpy(dist, src, size, cudaMemcpyHostToDevice))) {
    EXIT_ERR("Error CpyToDevice\n");
  }
}

template <typename distType, typename srcType, typename sizeT>
void inline CpyToHost(distType *&dist, const srcType *src, sizeT size) {

  if (CheckError(cudaMemcpy(dist, src, size, cudaMemcpyDeviceToHost))) {
    EXIT_ERR("Error CpyToHost\n")
  }
}

template <typename distType, typename srcType, typename sizeT>
void inline CpyToDeviceAsync(distType *&dist, const srcType *src, sizeT size) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyHostToDevice))) {
    EXIT_ERR("Error CpyToDeviceAsync\n")
  }
}

template <typename distType, typename srcType, typename sizeT>
void inline CpyToHostAsync(distType *&dist, const srcType *src, sizeT size) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyDeviceToHost))) {
    EXIT_ERR("Error CpyToHostAsync\n")
  }
}

template <typename distType, typename srcType, typename sizeT>
void inline CpyToHostAsyncStream(distType *&dist, const srcType *src, sizeT size, cudaStream_t &st) {

  if (CheckError(cudaMemcpyAsync(dist, src, size, cudaMemcpyDeviceToHost, st))) {
    EXIT_ERR("Error CpyToHostAsyncStream\n")
  }
}
template <typename srcType, typename sizeT>
void inline MallocHost(sizeT size, const srcType *&src) {
  if (CheckError(cudaMallocHost(dist, (void **)src, size))) {
    EXIT_ERR("Error MallocHost\n")
  }
}

template <typename srcType>
void inline FreeMem(srcType *&src) {
  if (CheckError(cudaFree(src))) {
    EXIT_ERR("Error FreeMem\n")
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