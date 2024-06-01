#ifdef USE_CUDA
#include "cuda_ray_interface.h"
#include "global_value.h"
#include "ray_tracing_build_plane.h"
#include "ray_tracing_calc_illum.h"
#include "ray_tracing_const.h"
#include "ray_tracing_main.h"

#include "reader_bin.h"
#include "writer_bin.h"

static uint8_t it = 0;

int ray_tracing::FindIntersections() {

  ///\note Этот луч задает центр проекции

#ifdef GRB_TASK
  Vector3 start_ray;
  Vector3 end_ray;
  if (it == 0) {
    start_ray = Vector3(1.1, 0, 0);
    end_ray = Vector3(0.95, 0, 0);
  } else {
    start_ray = Vector3(0.5, 0, 1);
    end_ray = Vector3(0.5, 0, 0);
  }
#else
  const Vector3 start_ray = Vector3(0.5, 0, 1);
  const Vector3 end_ray(0.5, 0, 0);
#endif

  int myid = get_mpi_id();
  int np = get_mpi_np();

  if (get_mpi_id() != 0) {
    return e_completion_success;
  }

  std::vector<FaceCell> faces; ///< поверхность сетки
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE, faces))
    return e_completion_fail;

  std::vector<Ray_t> rays(k_pixels_width * k_pixels_height);
  std::vector<IntId> intersections(k_pixels_width * k_pixels_height, -1);

  cuda::ray_tracing::interface::InitDevice(faces, rays);

  Ray_t center_ray(start_ray, (end_ray - start_ray).normalized()); ///< центр плоскости
  WRITE_LOG("Ray_center: [%lf, %lf, %lf] -> [%lf %lf %lf]\n", start_ray[0], start_ray[1], start_ray[2], center_ray.direction[0], center_ray.direction[1], center_ray.direction[2]);

  Timer time;
  MakeRays(center_ray, rays);
  cuda::ray_tracing::interface::StartTracingGrid(rays, intersections);

  WRITE_LOG("Tracing end %lf\n", time.get_delta_time_sec());

  files_sys::bin::WriteSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(it++) + ".bin", intersections);
  cuda::ray_tracing::interface::ClearDevice();

  MPI_BARRIER(MPI_COMM_WORLD);
  return e_completion_success;
}
#endif ///! USE_CUDA