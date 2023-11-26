#ifdef USE_CUDA

#include "cuda_integrator.h"
#include "cuda_scattering.h"
#include "global_def.h"

__global__ void cuda::kernel::GetS(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const IdType M = dir->size;
  const IdType N = grid->size;

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType k = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= N || k >= M)
    return;

  const Vector3 &cur_dir = dir->directions[k].dir;
  const Type *Illum = grid->illum;
  const geo::direction_device_t *all_dir = dir->directions;

  Type scatter = 0;
  for (IdType num_direction = 0; num_direction < M; num_direction++) {
    IdType pos = CELL_SIZE * (num_direction * N + i);
    Type I = (Illum[pos] + Illum[pos + 1] + Illum[pos + 2] + Illum[pos + 3]) / 4.;
    scatter += device::direction_integrator::Gamma(all_dir[num_direction].dir, cur_dir) * I * all_dir[num_direction].area;
  }

  grid->int_scattering[k * N + i] = scatter / dir->full_area;
}

__global__ void cuda::kernel::GetS_MPI(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const IdType N = grid->size;

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType k = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= N || k >= grid->local_scattering_size)
    return;

  const IdType M = dir->size;

  const Vector3 &cur_dir = dir->directions[grid->local_scattering_disp + k].dir;
  const Type *Illum = grid->illum;
  const geo::direction_device_t *all_dir = dir->directions;

  Type scatter = 0;
  for (IdType num_direction = 0; num_direction < M; num_direction++) {
    IdType pos = CELL_SIZE * (num_direction * N + i);
    Type I = (Illum[pos] + Illum[pos + 1] + Illum[pos + 2] + Illum[pos + 3]) / 4.;
    scatter += device::direction_integrator::Gamma(all_dir[num_direction].dir, cur_dir) * I * all_dir[num_direction].area;
  }

  grid->int_scattering[k * N + i] = scatter / dir->full_area;
}

__global__ void cuda::kernel::GetS_MPI_Stream(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid, const IdType start, const IdType end) {
  const IdType N = grid->size;

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType k = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= N || k >= end || k < start)
    return;

  const IdType M = dir->size;

  const Vector3 &cur_dir = dir->directions[grid->local_scattering_disp + k].dir;
  const Type *Illum = grid->illum;
  const geo::direction_device_t *all_dir = dir->directions;

  Type scatter = 0;
  for (IdType num_direction = 0; num_direction < M; num_direction++) {
    IdType pos = CELL_SIZE * (num_direction * N + i);
    Type I = (Illum[pos] + Illum[pos + 1] + Illum[pos + 2] + Illum[pos + 3]) / 4.;
    scatter += device::direction_integrator::Gamma(all_dir[num_direction].dir, cur_dir) * I * all_dir[num_direction].area;
  }

  grid->int_scattering[k * N + i] = scatter / dir->full_area;
}

__global__ void cuda::kernel::GetS_MPI_multi_device(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid,
                                                    const IdType size_loc, const IdType start_dir, const IdType end_dir) {
  const IdType N = grid->size;

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType k = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= size_loc || k >= end_dir || k < start_dir)
    return;

  const IdType M = dir->size;

  const Vector3 &cur_dir = dir->directions[grid->local_scattering_disp + k].dir;
  const Type *Illum = grid->illum;
  const geo::direction_device_t *all_dir = dir->directions;

  Type scatter = 0;
  for (IdType num_direction = 0; num_direction < M; num_direction++) {
    Type I = Illum[i * M + num_direction];
    scatter += device::direction_integrator::Gamma(all_dir[num_direction].dir, cur_dir) * I * all_dir[num_direction].area;
  }

  grid->int_scattering[k * N + i] = scatter / dir->full_area;
}

#endif //! USE_CUDA