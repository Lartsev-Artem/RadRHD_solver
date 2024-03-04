#ifdef USE_CUDA
#include "cuda_illum_sep_param.h"
#ifdef ON_FULL_ILLUM_ARRAYS

#include "global_def.h"

namespace cuda_sep = cuda::separate_device;

__device__ void cuda_sep::device::MakeEnergy(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType N = grid->loc_size_params;
  const IdType M = dir->size;

  if (i >= N)
    return;

  IdType Illum_shift = grid->shift_params; // откуда начинаются данные

  Type sum = 0;
  for (IdType k = 0; k < M; k++) {
    sum += grid->illum[M * (Illum_shift + i) + k] * dir->directions[k].area;
  }

  grid->energy[i] = sum / dir->full_area; // direction_integrator::IntegrateByCell(shift + i, dir, grid);
}

__device__ void cuda_sep::device::MakeStream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType N = grid->loc_size_params;
  const IdType M = dir->size;

  if (i >= N)
    return;

  IdType Illum_shift = grid->shift_params; // откуда начинаются данные

  Vector3 stream = Vector3::Zero();
  const Type *I = grid->illum + M * (i + Illum_shift);
  for (IdType k = 0; k < M; k++) {
    stream += dir->directions[k].dir * (I[k] * dir->directions[k].area);
  }

  grid->stream[i] = stream / dir->full_area;
}

__device__ void cuda_sep::device::MakeImpuls(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const IdType cell = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType N = grid->loc_size_params; //сколько хотим посчитать
  const IdType M = dir->size;

  if (cell >= N)
    return;

  IdType Illum_shift = grid->shift_params; // откуда начинаются данные

  Matrix3 impuls = Matrix3::Zero();
  const Type *illum = grid->illum + M * (cell + Illum_shift);

  for (IdType d = 0; d < M; d++) {

    const Type Ia = illum[d] * dir->directions[d].area;
    for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++) {
        impuls(i, k) += dir->directions[d].dir[i] * dir->directions[d].dir[k] * Ia;
      }
  }
  grid->impuls[cell] = impuls / dir->full_area;
}

/// \note если данные на ячейках, то дивергенции на прямую не посчитать (для rad_rhd они не нужны)
__global__ void cuda_sep::kernel::MakeIllumParam(const cuda::geo::grid_directions_device_t *dir, cuda::geo::grid_device_t *grid, IdType size_params, IdType shift_params) {
  // эти функции можно объденить в одну. Тогда будет одно общее обращение в память к illum

#ifdef SINGLE_GPU
  grid->loc_size_params = size_params;
  grid->shift_params = shift_params;
#endif

  cuda_sep::device::MakeEnergy(dir, grid);
  cuda_sep::device::MakeStream(dir, grid);
  cuda_sep::device::MakeImpuls(dir, grid);
}

#endif
#endif //! USE_CUDA