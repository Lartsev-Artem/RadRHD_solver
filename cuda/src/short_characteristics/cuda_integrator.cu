#ifdef USE_CUDA
//***********************************************************************//
//*********************Functions from device*****************************//
//***********************************************************************//
#include "cuda_integrator.h"
#include "global_def.h"

namespace c_dir = cuda::device::direction_integrator;

__device__ Type c_dir::Gamma(const Vector3 &direction, const Vector3 &direction2) {
  Type sum = direction.dot(direction2);
  return (3. * (1 + sum * sum)) / 4.;
}

__device__ Type c_dir::IntegrateByCell(const int num_cell, const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const int M = dir->size;
  const int N = grid->size;

  Type res = 0;
  for (int i = 0; i < M; i++) {
    int pos = CELL_SIZE * (N * i + num_cell);

    Type I = 0;
    for (int k = 0; k < CELL_SIZE; k++) {
      I += grid->illum[pos + k];
    }
    I /= CELL_SIZE;

    res += I * dir->directions[i].area;
  }

  return res / dir->full_area;
}

__device__ void c_dir::IntegrateByFaces3(const int num_cell, const geo::grid_directions_device_t *dir_grid, geo::grid_device_t *grid, Vector3 *Stream) {

  const int M = dir_grid->size;
  const int N = grid->size;

  for (int h = 0; h < CELL_SIZE; h++) {
    Stream[h] = Vector3::Zero();
  }

  for (int i = 0; i < M; i++) {
    int pos = CELL_SIZE * (N * i + num_cell);
    for (int f = 0; f < CELL_SIZE; f++) {
      Stream[f] += dir_grid->directions[i].dir * (grid->illum[pos + f] * dir_grid->directions[i].area);
    }
  }

  for (int h = 0; h < CELL_SIZE; h++)
    Stream[h] /= dir_grid->full_area;

  return;
}

__device__ void c_dir::IntegrateByFaces9(const int num_cell, const geo::grid_directions_device_t *dir_grid, geo::grid_device_t *grid, Matrix3 *Impuls) {
  const int M = dir_grid->size;
  const int N = grid->size;

  for (int h = 0; h < CELL_SIZE; h++) {
    Impuls[h] = Matrix3::Zero();
  }

  for (int dir = 0; dir < M; dir++) {
    int pos = CELL_SIZE * (N * dir + num_cell);
    for (int f = 0; f < CELL_SIZE; f++) {
      for (int i = 0; i < 3; i++)
        for (int k = 0; k < 3; k++) {
          Type I = grid->illum[pos + f];
          Impuls[f][i * 3 + k] += dir_grid->directions[dir].dir[i] * dir_grid->directions[dir].dir[k] * (I * dir_grid->directions[dir].area);
        }
    }
  }

  for (int h = 0; h < CELL_SIZE; h++)
    Impuls[h] /= dir_grid->full_area;

  return;
}

#endif // USE_CUDA