#if !defined CUDA_RAY_INTERFACE_H && defined USE_CUDA
#define CUDA_RAY_INTERFACE_H

#include "geo_types.h"

namespace cuda::ray_tracing {

enum e_ray_intersect_code {
  e_ray_intersect_none = -1,
  e_ray_intersect_disk = -2,
  e_ray_intersect_sphere = -3,
  e_ray_intersect_rosh = -4,
};

namespace interface {

int StartTracing(const std::vector<Ray_t> &rays_host, std::vector<IntId> &intersections);

void InitDevice(const std::vector<FaceCell> &faces_host, const std::vector<Ray_t> &rays_host);

void ClearDevice();

} // namespace interface
} // namespace cuda::ray_tracing
#endif //! CUDA_RAY_INTERFACE_H