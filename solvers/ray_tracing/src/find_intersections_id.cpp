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

int ray_tracing::FindIntersections()
{

  ///\note Этот луч задает центр проекции

#if TASK_TYPE == GRB_TASK
  Vector3 start_ray;
  Vector3 end_ray;
  if (it == 0)
  {
    start_ray = Vector3(1.1, 0, 0);
    end_ray = Vector3(0.95, 0, 0);
  }
  else
  {
    start_ray = Vector3(0.5, 0, 1);
    end_ray = Vector3(0.5, 0, 0);
  }
#else
  const Vector3 start_ray = Vector3(0.5, 0, 1);
  const Vector3 end_ray(0.5, 0, 0);
#endif

  int myid = get_mpi_id();
  int np = get_mpi_np();

  if (get_mpi_id() != 0)
  {
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

int ray_tracing::FindObserverIntersections()
{

  int myid = get_mpi_id();
  int np = get_mpi_np();

  if (get_mpi_id() != 0)
  {
    return e_completion_success;
  }

  std::vector<FaceCell> faces; ///< поверхность сетки
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE, faces))
    return e_completion_fail;

  ParamTraceProjection plane_cfg(PlaneParams(1, 1, 1, 1), Vector3::Zero(), Vector3::Zero());
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_PLANE_CFG, (uint8_t *)(&plane_cfg)))
    return e_completion_fail;

  std::vector<Ray_t> rays;
  PlaneParams &plane_params = plane_cfg.params2D;
  Vector3 &plane_orig = plane_cfg.plane_orig;
  Vector3 &observer = plane_cfg.observer;
  MakeRays(plane_params, plane_orig, observer, rays);

  std::vector<IntId> intersections(rays.size(), -1);
  cuda::ray_tracing::interface::InitDevice(faces, rays);

  WRITE_LOG(R"(Tracing start:
            PixX= %d,PixY= %d, width= %lf, height= %lf,
            plane_orig: %lf %lf %lf,
            observer: %lf %lf %lf
            )",
            plane_params._pixels_width,
            plane_params._pixels_height,
            plane_params._width,
            plane_params._height,
            plane_orig[0], plane_orig[1], plane_orig[2],
            observer[0], observer[1], observer[2]);

  Timer time;
  cuda::ray_tracing::interface::StartTracingGrid(rays, intersections);

  WRITE_LOG("Tracing end %lf\n", time.get_delta_time_sec());

  files_sys::bin::WriteSimple(glb_files.trace_address + F_RAY_TRACE + ".bin", intersections);
  cuda::ray_tracing::interface::ClearDevice();

  MPI_BARRIER(MPI_COMM_WORLD);
  return e_completion_success;
}
#endif ///! USE_CUDA