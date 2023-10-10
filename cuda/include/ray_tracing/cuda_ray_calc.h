#if !defined CUDA_RAY_CALC_H && defined USE_CUDA
#define CUDA_RAY_CALC_H

#include "cuda_ray_struct.h"

namespace cuda::ray_tracing {

__device__ Type RayIntersectsTriangle(const Ray &ray, const Face &triangle, Vector3 &intersection);

__global__ void RayTracing(const int M, const Ray *rays, const int N, const Face *triangles, int *intersections);

} // namespace cuda::ray_tracing
#endif //! CUDA_RAY_CALC_H