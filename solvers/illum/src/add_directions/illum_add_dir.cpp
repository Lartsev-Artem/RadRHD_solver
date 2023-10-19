#if defined ILLUM
#include "illum_add_dir.h"
#include "global_value.h"
#include "ray_tracing_const.h"

#include "reader_bin.h"
#include "writer_bin.h"
#include "writer_txt.h"

#ifdef USE_CUDA
#include "cuda_ray_interface.h"
#else
#include "intersections.h"
#include <omp.h>
#endif

int illum::additional_direction::MakeDirectionReInterpolation(const std::string &base_address, const std::string &out_address,
                                                              const int add_directions, const grid_directions_t &grid_dir) {

  std::vector<Face> direction_surface;
  if (files_sys::bin::ReadSimple(base_address + F_SURFACE_SPHERE_DIRECTION, direction_surface))
    return e_completion_fail;

  DIE_IF(grid_dir.size != direction_surface.size()); //разные сферы направлений

  std::vector<Ray_t> rays(add_directions);
  std::vector<int> dir_intersection_id(add_directions);

  for (size_t i = 0; i < add_directions; i++) {

    Type start_angle = 2 * PI / add_directions * i; // текущий поворот картинной плоскости
    // выполняется поворот относительно центра масс на угол start_angle относительно нулевого положения
    const Vector3 start_ray = Vector3(cos(start_angle), sin(start_angle), ray_tracing::k_height_above_center) + ray_tracing::k_center_of_mass;
    const Vector3 end_ray(ray_tracing::k_center_of_mass);

    rays[i].direction = -(end_ray - start_ray).normalized();
    rays[i].orig = Vector3::Zero(); //сфера направлений задана в нуле
  }

#ifndef USE_CUDA

#pragma omp parallel for
  for (int dir = 0; dir < add_directions; dir++) {

    int loc_id = e_ray_intersect_none;
    Type loc_dist = 1e10;

    //по всей сетке
    Vector3 p;
    for (int i = 0; i < grid_dir.size; i++) {

      Type d = intersection::RayIntersectsTriangle(rays[dir].orig, rays[dir].direction, direction_surface[i], p);

      //новое пересечение есть и оно ближе имеющегося
      if (d > 0 && d < loc_dist) {
        loc_id = i;
        loc_dist = d;
      }
    }

    dir_intersection_id[dir] = loc_id;
  }
#else
  std::vector<FaceCell> faces(direction_surface.size());
  for (size_t i = 0; i < direction_surface.size(); i++) {
    faces[i].face = direction_surface[i];
    faces[i].face_id = i;
  }
  direction_surface.clear();

  cuda::ray_tracing::interface::InitDevice(faces, rays);
  cuda::ray_tracing::interface::TracingGrid(rays, dir_intersection_id);
  cuda::ray_tracing::interface::ClearDevice();
#endif //! USE_CUDA

  grid_directions_t grid_plane(add_directions);
  for (size_t i = 0; i < add_directions; i++) {
    grid_plane.directions[i].dir = rays[i].direction;
    grid_plane.directions[i].area = 1;
  }
  grid_plane.full_area = add_directions;

  files_sys::txt::WriteSphereDirectionCartesian(out_address + F_ADDITIONAL_DIRECTION_GRID, grid_plane);

  return files_sys::bin::WriteSimple(out_address + F_DIRECTION_INTERPOLATION, dir_intersection_id); //номера ячеек переинтерполяции для сетки направлений
}

void illum::additional_direction::SaveInterpolationScattering(const std::string &address_add_dir, const grid_t &grid) {
  std::vector<int> interpolation_direction;
  if (files_sys::bin::ReadSimple(address_add_dir + F_DIRECTION_INTERPOLATION, interpolation_direction) == e_completion_success) {
    for (size_t i = 0; i < interpolation_direction.size(); i++) {
      int dir = interpolation_direction[i];
      files_sys::bin::WriteSimple(address_add_dir + F_SCATTERING + std::to_string(i) + ".bin", grid.size, &grid.scattering[dir * grid.size]);
    }
  }
}
#endif //! ILLUM