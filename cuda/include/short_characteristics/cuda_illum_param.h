#if !defined CUDA_ILLUM_PARAM_H && defined USE_CUDA
#define CUDA_ILLUM_PARAM_H

#include "cuda_struct.h"

namespace cuda {

namespace device {
#ifdef ON_FULL_ILLUM_ARRAYS
__device__ void MakeEnergy(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);
#endif

__device__ void MakeDivStream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

__device__ void MakeDivImpuls(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);
} // namespace device

namespace kernel {

__global__ void MakeIllumParam(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

} // namespace kernel
} // namespace cuda
#endif //! CUDA_ILLUM_PARAM_H