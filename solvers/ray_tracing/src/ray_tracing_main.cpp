#include "ray_tracing_main.h"
#include "cuda_ray_interface.h"
#include "global_value.h"
#include "ray_tracing_build_plane.h"
#include "ray_tracing_const.h"

#include "reader_bin.h"
#include "writer_bin.h"

/// \todo Отрисовка самой плоскости и её перевод в vtk формат
int ray_tracing::RunRayTracing() {

#ifndef USE_CUDA
  EXIT_ERR("Need cuda support\n");
#else

  std::vector<FaceCell> faces;
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE, faces))
    return e_completion_fail;

  std::vector<Ray_t> rays(k_pixels_width * k_pixels_height);
  std::vector<IntId> intersections(k_pixels_width * k_pixels_height, -1);

  cuda::ray_tracing::interface::InitDevice(faces, rays);

  for (int i = 0; i < k_number_of_frame; i++) {
    MakePlane(i, rays);
    cuda::ray_tracing::interface::StartTracing(rays, intersections);
    files_sys::bin::WriteSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(i) + ".bin", intersections);
  }

  cuda::ray_tracing::interface::ClearDevice();

#endif // USE_CUDA
  return e_completion_success;
}