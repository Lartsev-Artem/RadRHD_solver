#ifdef USE_CUDA
#include "cuda_def.h"

#include "cuda_memory.h"
#include "cuda_ray_calc.h"
#include "cuda_ray_init.h"
#include "cuda_ray_interface.h"

namespace mem = cuda::mem_protected;
namespace ray = cuda::ray_tracing;

static int size_surface = 0;
static ray::Face *faces_device;
static ray::Ray *rays_device;
static int *intersection_device;

int ray::interface::StartTracing(const std::vector<Ray_t> &rays_host, std::vector<int> &intersections) {

  int M = rays_host.size();

  CUDA_TREADS_1D(threads);
  CUDA_BLOCKS_1D(blocks, M);

  mem::CpyToDevice(rays_device, rays_host.data(), M * sizeof(Ray_t));

  cuda::ray_tracing::RayTracing<<<blocks, threads>>>(M, rays_device, size_surface, faces_device, intersection_device);

  CUDA_CALL_FUNC(cudaGetLastError);

  mem::CpyToHost(intersections.data(), intersection_device, M * sizeof(int));

  return e_completion_success;
}

int ray::interface::FindInnerIntersection(const std::vector<Ray_t> &rays_host, std::vector<int> &intersections) {

  int M = rays_host.size();

  CUDA_TREADS_1D(threads);
  CUDA_BLOCKS_1D(blocks, M);

  mem::CpyToDevice(rays_device, rays_host.data(), M * sizeof(Ray_t));

  cuda::ray_tracing::InnerRayTracing<<<blocks, threads>>>(M, rays_device, size_surface, faces_device, intersection_device);

  CUDA_CALL_FUNC(cudaGetLastError);

  mem::CpyToHost(intersections.data(), intersection_device, M * sizeof(int));

  return e_completion_success;
}

void ray::interface::InitDevice(const std::vector<FaceCell> &faces_host, const std::vector<Ray_t> &rays_host) {

  CUDA_CALL_FUNC(cudaSetDevice, 0);
  size_surface = faces_host.size();
  ray::InitMemory(faces_host, faces_device, rays_host, rays_device, intersection_device);
}

void ray::interface::ClearDevice() {
  ray::ClearMemory(faces_device, rays_device, intersection_device);
  size_surface = 0;
  WRITE_LOG("Free device arrays\n");
}

#endif