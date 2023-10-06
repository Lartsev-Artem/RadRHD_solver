#ifdef USE_CUDA
#include "cuda_def.h"
#include "cuda_interface.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

void cuda::interface::CudaWait() {
  CUDA_CALL_FUNC(cudaDeviceSynchronize);
}

#endif //! USE_CUDA