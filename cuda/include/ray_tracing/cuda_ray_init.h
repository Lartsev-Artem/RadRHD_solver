#if !defined CUDA_RAY_INIT_H && defined USE_CUDA
#define CUDA_RAY_INIT_H

#include "cuda_ray_struct.h"
#include "geo_types.h"

namespace cuda::ray_tracing {

int InitMemory(const std::vector<FaceCell> &faces_host, Face *&faces_device,
               const std::vector<Ray_t> &rays_host, Ray *&rays_device,
               int *&intersection_device);

int ClearMemory(Face *&faces_device, Ray *&rays_device, int *&intersection_device);

} // namespace cuda::ray_tracing
#endif //! CUDA_RAY_INIT_H