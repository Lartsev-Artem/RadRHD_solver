#if defined USE_CUDA

#include "cuda_memory.h"
#include "cuda_ray_init.h"

namespace mem = cuda::mem_protected;

int cuda::ray_tracing::InitMemory(const std::vector<FaceCell> &faces_host, Face *&faces_device,
                                  const std::vector<Ray_t> &rays_host, Ray *&rays_device,
                                  int *&intersection_device) {

  mem::Malloc(faces_host.size() * sizeof(FaceCell), &faces_device);
  mem::CpyToDevice(faces_device, faces_host.data(), faces_host.size() * sizeof(FaceCell));

  mem::Malloc(rays_host.size() * sizeof(Ray_t), &rays_device);       //матрица заданного размера
  mem::Malloc(rays_host.size() * sizeof(int), &intersection_device); // 1 луч - 1 пересечение

  return e_completion_success;
}

int cuda::ray_tracing::ClearMemory(Face *&faces_device, Ray *&rays_device, int *&intersection_device) {
  mem::FreeMem(faces_device);
  mem::FreeMem(rays_device);
  mem::FreeMem(intersection_device);
  return e_completion_success;
}

#endif //! USE_CUDA