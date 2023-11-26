#ifdef USE_CUDA
#include "cuda_illum_sep_param.h"

#include "global_def.h"

namespace cuda_sep = cuda::separate_device;

__device__ void cuda_sep::device::MakeEnergy(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType N = grid->loc_size;
  const IdType M = dir->size;

  if (i >= N)
    return;
  if (i == 0) {
    printf("%ld %ld\n", N, M);
  }

  Type sum = 0;
  for (IdType k = 0; k < M; k++) {
    sum += grid->illum[M * i + k] * dir->directions[k].area;
  }

  grid->energy[i] = sum / dir->full_area; // direction_integrator::IntegrateByCell(shift + i, dir, grid);
}

///\todo наполнение функций
__device__ void cuda_sep::device::MakeStream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
}
__device__ void cuda_sep::device::MakeImpuls(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
}

/// \note если данные на ячейках, то дивергенции на прямую не посчитать (для rad_rhd они не нужны)
__global__ void cuda_sep::kernel::MakeIllumParam(const cuda::geo::grid_directions_device_t *dir, cuda::geo::grid_device_t *grid) {
  // эти функции можно объденить в одну. Тогда будет одно общее обращение в память к illum
  cuda_sep::device::MakeEnergy(dir, grid);
  cuda_sep::device::MakeStream(dir, grid);
  cuda_sep::device::MakeImpuls(dir, grid);
}

#endif //! USE_CUDA