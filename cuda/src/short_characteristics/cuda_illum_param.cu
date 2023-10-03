#ifdef USE_CUDA
#include "cuda_illum_param.h"
#include "cuda_integrator.h"

#include "global_def.h"

#ifdef ON_FULL_ILLUM_ARRAYS
#define CUDA_CONVERT_FACE_TO_CELL(val, size, src) \
  for (int k = 0; k < size; k++) {                \
    val[k] = 0;                                   \
    for (int f = 0; f < CELL_SIZE; f++)           \
      val[k] += src[f][k];                        \
    val[k] /= CELL_SIZE;                          \
  }

__device__ void cuda::kernel::MakeEnergy(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const int N = grid->loc_size;
  const int shift = grid->shift;
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= N)
    return;

  grid->energy[i] = direction_integrator::IntegrateByCell(shift + i, dir, grid);
}
#endif // ON_FULL_ILLUM_ARRAYS

__device__ void cuda::kernel::MakeDivStream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const int M = dir->size;
  const int N = grid->loc_size;
  const int shift = grid->shift;

  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= N)
    return;

  Vector3 Stream[CELL_SIZE];
  direction_integrator::IntegrateByFaces3(i + shift, dir, grid, Stream);

#ifdef ON_FULL_ILLUM_ARRAYS
  CUDA_CONVERT_FACE_TO_CELL(grid->stream[i], 3, Stream);
#endif

  grid->divstream[i] = 0;
  int pos = (i + shift) * CELL_SIZE;
  Type div = 0;
  for (int f = 0; f < CELL_SIZE; f++) {
    Type sum = 0;
    for (int k = 0; k < 3; k++) {
      sum += Stream[f][k] * grid->normals[pos + f][k];
    }
    div += sum * grid->areas[pos + f];
  }

  grid->divstream[i] = div / grid->volume[i + shift];
  return;
}

__device__ void cuda::kernel::MakeDivImpuls(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const int M = dir->size;
  const int N = grid->loc_size;
  const int shift = grid->shift;

  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= N)
    return;

  Matrix3 impuls[CELL_SIZE];
  direction_integrator::IntegrateByFaces9(i + shift, dir, grid, impuls);

#ifdef ON_FULL_ILLUM_ARRAYS
  CUDA_CONVERT_FACE_TO_CELL(grid->impuls[i], 9, impuls);
#endif

  Vector3 div = Vector3::Zero();

  for (int j = 0; j < CELL_SIZE; j++) {
    int pos = (i + shift) * CELL_SIZE + j;
    for (int h = 0; h < 3; h++) {
      Type sum = 0;
      for (int k = 0; k < 3; k++) {
        sum += impuls[j][h * 3 + k] * grid->normals[pos][k];
      }

      div[h] += sum * grid->areas[pos];
    }
  }

  grid->divimpuls[i] = div / grid->volume[i + shift];
  return;
}

#undef CUDA_CONVERT_FACE_TO_CELL

__global__ void cuda::MakeIllumParam(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  // эти функции можно объденить в одну. Тогда будет одно общее обращение в память к illum
  kernel::MakeEnergy(dir, grid);
  kernel::MakeDivStream(dir, grid);
  kernel::MakeDivImpuls(dir, grid);
}

#endif //! USE_CUDA