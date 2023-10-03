#if !defined CUDA_INTEGRATOR_H && defined USE_CUDA
#define CUDA_INTEGRATOR_H

#include "cuda_struct.h"
namespace cuda::kernel {
namespace direction_integrator {

__device__ Type Gamma(const Vector3 &direction, const Vector3 &direction2);

__device__ Type IntegrateByCell(const int num_cell, const geo::grid_directions_device_t *dir, geo::grid_device_t *grid);

__device__ void IntegrateByFaces3(const int num_cell, const geo::grid_directions_device_t *dir_grid, geo::grid_device_t *grid, Vector3 *Stream);

__device__ void IntegrateByFaces9(const int num_cell, const geo::grid_directions_device_t *dir_grid, geo::grid_device_t *grid, Matrix3 *Impuls);

} // namespace direction_integrator
} // namespace cuda::kernel
#endif // CUDA_INTEGRATOR_H